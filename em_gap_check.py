"""接触部位の EM を切り出し、2 つの細胞の「間に何があるか」を数える。

  python em_gap_check.py --pairs results/mb_em_targets.csv

骨格ベースの「接触 1 µm 以内」では、隙間にグリアや別のニューロンが入っているかどうかは
判定できない。そこで 2 段階で実体を確かめる。

  1. メッシュ（実際の膜）どうしの最短距離を測り、最接近点の座標を出す
  2. その座標でセグメンテーションを切り出し、A のボクセルと B のボクセルが
     直接隣接しているか、隣接していなければ間にどの物体があるかを数える

hemibrain のセグメンテーションは全ボクセルにいずれかの物体 ID が割り当てられているので、
A と B のボクセルが隣接していれば膜が直接向かい合っている。

入力 CSV の列: t1, t2, pre, post （pre/post は hemibrain の bodyId）
出力: results/em_gap/<t1>_<t2>.txt と、まとめの results/em_gap/summary.csv
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

from common.runner import ROOT

VOX = 8.0  # nm / voxel（hemibrain は 8 nm 等方）
SEG = "precomputed://gs://neuroglancer-janelia-flyem-hemibrain/v1.2/segmentation"


def meshes(body_ids: list[int], lod: int = 2) -> dict[int, np.ndarray]:
    from navis.interfaces import neuprint as nvn
    from neuprint import NeuronCriteria as NC

    from common.loaders import _retry, neuprint_client
    client = neuprint_client("hemibrain:v1.2.1")
    out = {}
    for b in body_ids:
        m = _retry(nvn.fetch_mesh_neuron, NC(bodyId=int(b)), client=client, lod=lod)
        m = m[0] if hasattr(m, "__len__") else m
        out[int(b)] = np.asarray(m.vertices, float) * VOX
    return out


def closest_spots(va: np.ndarray, vb: np.ndarray, k: int = 3, apart_nm: float = 3000.0):
    """最接近点を、互いに apart_nm 以上離れた代表点として k 個返す。"""
    t = cKDTree(vb)
    d, j = t.query(va, k=1)
    picks = []
    for i in np.argsort(d)[:600]:
        p = (va[i] + vb[j[i]]) / 2
        if all(np.linalg.norm(p - q) > apart_nm for _, q in picks):
            picks.append((float(d[i]), p))
        if len(picks) == k:
            break
    return picks, d


def gap_at(vol, a_id: int, b_id: int, pos_vox: np.ndarray, radius=(160, 160, 80)):
    """その座標でセグメンテーションを切り出し、直接隣接と間にある物体を数える。"""
    sl = tuple(slice(int(pos_vox[i]) - radius[i], int(pos_vox[i]) + radius[i]) for i in range(3))
    seg = np.asarray(vol[sl])[..., 0]
    A, B = seg == a_id, seg == b_id
    if not A.any() or not B.any():
        return None
    st = ndimage.generate_binary_structure(3, 1)
    touch = int((ndimage.binary_dilation(A, st) & B).sum())
    dist, idx = ndimage.distance_transform_edt(~B, return_indices=True)
    surf = A & ndimage.binary_dilation(~A, st)
    pts = np.argwhere(surf & (dist <= 12))  # 96 nm 以内まで
    if not len(pts):
        return {"touch": touch, "n": 0, "direct": 0, "between": {}}
    pts = pts[:: max(1, len(pts) // 4000)]
    from collections import Counter
    between, direct = Counter(), 0
    for q in pts:
        tgt = idx[(slice(None), *q)]
        n = int(max(abs(tgt - q)))
        if n <= 1:
            direct += 1
            continue
        got = {int(seg[tuple(np.round(q + (tgt - q) * s / n).astype(int))]) for s in range(1, n)}
        got -= {a_id, b_id, 0}
        if got:
            between.update(got)
        else:
            direct += 1
    return {"touch": touch, "n": len(pts), "direct": direct, "between": dict(between.most_common(6))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True, help="列 t1,t2,pre,post の CSV")
    ap.add_argument("--spots", type=int, default=2, help="1 組あたり調べる最接近点の数")
    args = ap.parse_args()

    from cloudvolume import CloudVolume
    vol = CloudVolume(SEG, mip=0, use_https=True, progress=False, fill_missing=True)

    pairs = pd.read_csv(args.pairs)
    ids = sorted({int(b) for b in list(pairs.pre) + list(pairs.post)})
    print(f"メッシュを取得: {len(ids)} 細胞")
    V = meshes(ids)
    for b, v in V.items():
        print(f"  {b} 頂点 {len(v)}")

    out_dir = ROOT / "results" / "em_gap"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in pairs.itertuples(index=False):
        a, b = int(r.pre), int(r.post)
        picks, d = closest_spots(V[a], V[b], k=args.spots)
        near = {f"<{t}nm": int((d < t).sum()) for t in (50, 100, 200, 500)}
        print(f"\n########## {r.t1} → {r.t2}  (body {a} → {b})")
        print(f"  膜の距離: " + "  ".join(f"{k} {v}" for k, v in near.items()))
        for gap_nm, p in picks:
            pv = (p / VOX).round().astype(int)
            g = gap_at(vol, a, b, pv)
            if g is None:
                print(f"  voxel {list(pv)}: 切り出し範囲にどちらかが無い")
                continue
            pct = g["direct"] * 100 // max(g["n"], 1)
            print(f"  voxel {pv[0]},{pv[1]},{pv[2]}  隙間 {gap_nm:.0f}nm  "
                  f"直接隣接 {g['touch']} voxel  間に何も無い表面点 {pct}% ({g['direct']}/{g['n']})")
            if g["between"]:
                print("    間にあった物体: " + ", ".join(f"{k}×{v}" for k, v in g["between"].items()))
            rows.append({"t1": r.t1, "t2": r.t2, "pre": a, "post": b,
                         "x": int(pv[0]), "y": int(pv[1]), "z": int(pv[2]),
                         "gap_nm": round(gap_nm), "touch_voxels": g["touch"],
                         "direct_pct": pct, **near,
                         "between": ";".join(str(k) for k in g["between"])})
    d = pd.DataFrame(rows)
    f = out_dir / "summary.csv"
    if f.exists():  # 前回までの結果を消さないよう、同じ組だけ入れ替えて残す
        old = pd.read_csv(f)
        keep = old[~old.set_index(["t1", "t2"]).index.isin(d.set_index(["t1", "t2"]).index)]
        d = pd.concat([keep, d], ignore_index=True)
    d.to_csv(f, index=False, encoding="utf-8")
    print(f"\n-> {f}")


if __name__ == "__main__":
    main()
