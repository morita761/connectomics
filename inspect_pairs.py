"""候補のタイプの組について、「どの脳領域（ROI）で接触しているか」と「その領域で相手と結合しているか」を調べる。

  python inspect_pairs.py --data data/hemibrain_mb_all.pkl --kind output --top 40
  python inspect_pairs.py --data data/hemibrain_mb_all.pkl --pair PAM07 MBON01

各タイプの組の代表ペア（type_pairs.py の ex_a, ex_b）について:
  - 接触点（上流側の出力部位が下流側に 1 µm 以内）を ROI ごとに数える（接触長 µm）
  - 上流側のシナプス前部位が、どの ROI にあるか（上流側がそこで出力しているか）
  - 下流側のシナプス後部位が、どの ROI にあるか（下流側がそこで入力を受けているか）
ROI の形は neuPrint から取得し data/roi_meshes/ に保存する（2 回目以降は再利用）。

出力: results/<データ名>/inspect_<kind>.csv と、標準出力の一覧。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import trimesh
from scipy.spatial import cKDTree

from common.config import Thresholds
from common.geometry import cable_points
from common.model import Dataset
from common.runner import ROOT, save_merged

UM = 1000.0
ROIS = ["CA(R)", "PED(R)", "a1(R)", "a2(R)", "a3(R)", "a'1(R)", "a'2(R)", "a'3(R)", "b1(R)", "b2(R)",
        "b'1(R)", "b'2(R)", "g1(R)", "g2(R)", "g3(R)", "g4(R)", "g5(R)", "LH(R)", "SMP(R)", "SIP(R)",
        "CRE(R)", "SLP(R)", "SCL(R)", "ICL(R)", "PLP(R)", "AL(R)",
        # キノコ体以外の領域（中心複合体・視覚の投射先など）
        "EB", "FB", "PB", "NO", "BU(R)", "LAL(R)", "AOTU(R)", "PVLP(R)", "AVLP(R)", "WED(R)", "SPS(R)", "IB",
        # 前運動・聴覚・下行ニューロン周辺（LAL・PS・WED の解析用）
        "VES(R)", "GOR(R)", "SAD", "CAN(R)", "IPS(R)", "EPA(R)", "GNG", "ATL(R)", "FLA(R)", "PRW"]


def load_meshes(dataset="hemibrain:v1.2.1", voxel_nm=8.0) -> dict[str, trimesh.Trimesh]:
    d = ROOT / "data" / "roi_meshes"
    d.mkdir(parents=True, exist_ok=True)
    client = None
    meshes = {}
    for r in ROIS:
        f = d / (r.replace("'", "p").replace("(", "_").replace(")", "") + ".obj")
        if not f.exists():
            if client is None:
                from common.loaders import neuprint_client
                client = neuprint_client(dataset)
            f.write_bytes(client.fetch_roi_mesh(r))
        m = trimesh.load(f, force="mesh")
        m.apply_scale(voxel_nm)  # ボクセル → nm
        meshes[r] = m
    return meshes


def roi_of(points: np.ndarray, meshes: dict) -> np.ndarray:
    """点ごとの ROI 名（キノコ体の区画を優先。どれにも入らなければ 'other'）。"""
    out = np.full(len(points), "other", dtype=object)
    if len(points) == 0:
        return out
    for r, m in meshes.items():  # ROIS の順 = 細かい区画が先
        todo = out == "other"
        if not todo.any():
            break
        lo, hi = m.bounds
        box = todo & np.all((points >= lo) & (points <= hi), axis=1)
        idx_box = np.flatnonzero(box)
        for k in range(0, len(idx_box), 2000):  # メモリを抑えるため少しずつ
            sub = idx_box[k:k + 2000]
            out[sub[m.contains(points[sub])]] = r
    return out


def summarize_counts(labels, weights=None, top=4) -> str:
    s = pd.Series(weights if weights is not None else np.ones(len(labels)), index=labels).groupby(level=0).sum()
    s = s.sort_values(ascending=False)
    tot = s.sum()
    if tot == 0:
        return "-"
    return ", ".join(f"{k} {v / tot:.0%}" for k, v in s.head(top).items())


def roi_synapse_counts(body_ids, dataset="hemibrain:v1.2.1") -> pd.DataFrame:
    """neuPrint の「ニューロン × ROI ごとのシナプス数」（pre = 出力、post = 入力）。"""
    from neuprint import NeuronCriteria as NC, fetch_neurons
    from common.loaders import neuprint_client
    _, roi = fetch_neurons(NC(bodyId=[int(b) for b in body_ids]), client=neuprint_client(dataset))
    return roi.set_index(["bodyId", "roi"])[["pre", "post"]]


def _roi_summary(counts: pd.DataFrame, body: int, col: str, top=4) -> str:
    if body not in counts.index.get_level_values(0):
        return "-"
    s = counts.loc[body][col]
    s = s[s.index.isin(ROIS)].sort_values(ascending=False)
    tot = s.sum()
    return ", ".join(f"{k} {v / tot:.0%}" for k, v in s.head(top).items() if v > 0) if tot else "-"


def inspect_pair(ds: Dataset, a: int, b: int, th: Thresholds, meshes, counts: pd.DataFrame) -> dict:
    na, nb = ds.neurons[a], ds.neurons[b]
    ca = cable_points(na, th.cable_spacing_um * UM, th.default_radius_nm)
    cb = cable_points(nb, th.cable_spacing_um * UM, th.default_radius_nm)
    out = np.zeros(len(ca["points"]), bool)
    if na.n_pre:
        d, _ = cKDTree(na.pre).query(ca["points"], distance_upper_bound=th.output_site_radius_um * UM)
        out = np.isfinite(d)
    d, _ = cKDTree(cb["points"]).query(ca["points"], distance_upper_bound=th.skeleton_contact_um * UM)
    hit = out & np.isfinite(d)
    pts, w = ca["points"][hit], ca["weight"][hit] / UM
    lab = roi_of(pts, meshes)
    main_roi = pd.Series(w, index=lab).groupby(level=0).sum().idxmax() if len(w) else None
    # 主な接触領域での、上流の出力数・下流の入力数（その領域で両者とも「働いている」か）
    get = lambda body, col: int(counts[col].get((body, main_roi), 0)) if main_roi else 0
    return {
        "contact_um_g": float(w.sum()),
        "contact_rois": summarize_counts(lab, w),
        "main_roi": main_roi,
        "pre_out_in_main": get(a, "pre"),
        "post_in_in_main": get(b, "post"),
        "pre_output_rois": _roi_summary(counts, a, "pre"),
        "post_input_rois": _roi_summary(counts, b, "post"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--kind", default="output", choices=["output", "directed", "undirected"])
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--all", action="store_true", help="feedforward・記憶系で絞らない")
    ap.add_argument("--pairs", default=None, help="調べる組を書いた CSV（列 t1, t2）")
    ap.add_argument("--ranked", action="store_true", help="rank_candidates.py の candidates_ranked.csv の上位 --top 件を調べる")
    ap.add_argument("--any-region", action="store_true", help="記憶系（キノコ体）でなくてもよい（上流→下流だけで絞る）")
    ap.add_argument("--pair", nargs=2, default=None, metavar=("TYPE1", "TYPE2"))
    args = ap.parse_args()

    th = Thresholds()
    ds = Dataset.load(args.data)
    res = ROOT / "results" / Path(args.data).stem
    tp = pd.read_csv(res / f"type_pairs_{args.kind}.csv")
    if args.pairs:  # 組を書いた CSV（列 t1, t2）。shortlist.csv などを渡せる
        want = set(zip(*pd.read_csv(args.pairs)[["t1", "t2"]].values.T))
        tp = tp[[(a, b) in want for a, b in zip(tp.t1, tp.t2)]]
    elif args.ranked:  # rank_candidates.py の並び順で上から
        tp = pd.read_csv(res / "candidates_ranked.csv").head(args.top)
    elif args.pair:
        tp = tp[(tp.t1 == args.pair[0]) & (tp.t2 == args.pair[1])]
    else:
        sel = tp.flagged
        if not args.all and "feedforward" in tp:
            sel &= tp.feedforward
            if not args.any_region:
                sel &= tp.memory
        tp = tp[sel].head(args.top)
    tp = tp.dropna(subset=["ex_a"])
    meshes = load_meshes()
    counts = roi_synapse_counts(set(tp.ex_a.astype(int)) | set(tp.ex_b.astype(int)))
    rows = []
    for r in tp.itertuples(index=False):
        info = inspect_pair(ds, int(r.ex_a), int(r.ex_b), th, meshes, counts)
        rows.append({"t1": r.t1, "t2": r.t2, "g1": r.g1, "g2": r.g2, "contact_total_um": r.contact_total_um,
                     "synapses": r.synapses, "expected": r.expected, "ex_a": int(r.ex_a), "ex_b": int(r.ex_b),
                     "ex_weight": r.ex_weight, **info})
        print(f"{r.t1:>14s} → {r.t2:<12s} 接触 {r.contact_total_um:7.1f} µm  syn {r.synapses:4d} / 期待 {r.expected:6.1f} | "
              f"接触: {info['contact_rois']} | 主領域で 上流出力 {info['pre_out_in_main']} / 下流入力 {info['post_in_in_main']}")
    out = pd.DataFrame(rows)
    f = res / f"inspect_{args.kind}.csv"
    save_merged(out, f, ["t1", "t2"])
    print(f"-> {f}")


if __name__ == "__main__":
    main()
