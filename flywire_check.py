"""FlyWire（別個体・全脳）の公開データで、候補のタイプの組の結合を数える。

必要なファイル（どちらも公開データ。作業用フォルダに置く）:
  proofread_connections_783.feather          Zenodo 10676866（852 MB、神経どうしの結合表）
  Supplemental_file1_neuron_annotations.tsv  GitHub flyconnectome/flywire_annotations（32 MB、細胞タイプ）

  # 1) 1 回だけ: 大きいファイルをタイプ単位に集計して小さい表に落とす
  python flywire_check.py --build --connections <...>.feather --annotations <...>.tsv
  # 2) 以降は候補の照合だけ（大きいファイルは不要）
  python flywire_check.py --data data/hemibrain_cx.pkl

タイプ名の対応は注釈の hemibrain_type 列（無い行は cell_type）を使う。
hemibrain_type が "AVLP386,AVLP388" のように複数あるときは、どちらにも数える
（「シナプスが無い」と言い切るうえでは、多めに数える方が安全なため）。

注意: Schlegel et al. 2024 によると、hemibrain のタイプの約 3 分の 1 は FlyWire で確実に同定できない。
      照合できなかった組は total_syn に加算せず、detail に理由を書く。

出力: data/flywire_type_edges.csv.gz（集計表）と results/<データ名>/cross_flywire_783.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common.runner import ROOT

EDGES = ROOT / "data" / "flywire_type_edges.csv.gz"
DATASET = "flywire:783"
TAG = "flywire_783"


def type_map(annotations: str) -> pd.DataFrame:
    """root_id → (タイプ, 半球) の表。複数タイプの行は行を分ける。"""
    a = pd.read_csv(annotations, sep="\t", low_memory=False,
                    usecols=["root_id", "side", "cell_type", "hemibrain_type"])
    a["type"] = a.hemibrain_type.fillna(a.cell_type)
    a = a.dropna(subset=["type"])
    a["type"] = a.type.astype(str).str.split(",")
    a = a.explode("type")
    a["type"] = a.type.str.strip()
    a["side"] = a.side.map({"left": "L", "right": "R"}).fillna("?")
    return a[["root_id", "type", "side"]].drop_duplicates()


def build(connections: str, annotations: str) -> pd.DataFrame:
    import pyarrow.feather as ft
    tm = type_map(annotations)
    print(f"注釈 {len(tm)} 行、タイプ {tm.type.nunique()} 種")
    cols = ["pre_pt_root_id", "post_pt_root_id", "neuropil", "syn_count"]
    c = ft.read_table(connections, columns=cols, memory_map=True).to_pandas()
    print(f"結合表 {len(c)} 行")
    pre = tm.rename(columns={"root_id": "pre_pt_root_id", "type": "t1", "side": "s1"})
    post = tm.rename(columns={"root_id": "post_pt_root_id", "type": "t2", "side": "s2"})
    c = c.merge(pre, on="pre_pt_root_id", how="inner").merge(post, on="post_pt_root_id", how="inner")
    print(f"タイプが付いた行 {len(c)}")
    g = (c.groupby(["t1", "t2", "s1", "s2", "neuropil"], observed=True)
         .syn_count.sum().reset_index())
    EDGES.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(EDGES, index=False, compression="gzip")
    print(f"{len(g)} 行 -> {EDGES}")
    return g


def load_edges() -> pd.DataFrame:
    if not EDGES.exists():
        raise SystemExit(f"{EDGES} がありません。先に --build を実行してください")
    return pd.read_csv(EDGES)


def check(data: str, top: int, min_visible: int, min_expected: float,
          pairs: str | None = None, extra_data: list[str] | None = None,
          tag: str | None = None) -> pd.DataFrame:
    from common.model import Dataset
    from validate_cross import reverse_in_source

    res = ROOT / "results" / (tag or Path(data).stem)
    res.mkdir(parents=True, exist_ok=True)
    # 組を外から指定できる（LexA を含む組だけ、など）。複数の領域にまたがってよい
    cand = pd.read_csv(pairs or (ROOT / "results" / Path(data).stem / "candidates_ranked.csv"))
    cand = cand[(cand.synapses == 0) & (cand.expected >= min_expected)]
    cand = cand.sort_values("expected", ascending=False).head(top)
    # 逆向きのシナプスは元データから数える。渡された全データを合わせる
    rev = {}
    for d in [data, *(extra_data or [])]:
        rev.update(reverse_in_source(Dataset.load(d)))

    e = load_edges()
    # 半球は「細胞体の左右」ではなく「シナプスがある領域の左右」で分ける。
    # MBON01 のように細胞体と樹状突起が反対側にある細胞があり、細胞体の左右で揃えると本物の結合を見落とす。
    e["nside"] = e.neuropil.astype(str).str.extract(r"_([LR])$")[0].fillna("other")
    fwd = e.groupby(["t1", "t2", "nside"], observed=True).syn_count.sum()
    total = e.groupby(["t1", "t2"], observed=True).syn_count.sum()  # 判定にはこちら（左右をまたぐ分も含む）
    # 「もし結合があれば見えるはずか」= そのタイプに出力・入力が見えているか
    out_w = e.groupby("t1").syn_count.sum()
    in_w = e.groupby("t2").syn_count.sum()
    types = set(e.t1) | set(e.t2)

    rows = []
    for r in cand.itertuples(index=False):
        src_fwd, src_rev = int(r.synapses), int(rev.get((r.t2, r.t1), 0))
        base = {"t1": r.t1, "t2": r.t2, "dataset": DATASET,
                "contact_um_src": r.contact_total_um, "expected_src": r.expected,
                "syn_src_fwd": src_fwd, "syn_src_rev": src_rev,
                "comp_type": r.comp_type, "comp_rate": r.comp_rate,
                "availability": getattr(r, "availability", np.nan)}
        missing = [t for t in (r.t1, r.t2) if t not in types]
        thin = []
        if not missing:
            if out_w.get(r.t1, 0) < min_visible:
                thin.append(f"{r.t1} の出力が見えない")
            if in_w.get(r.t2, 0) < min_visible:
                thin.append(f"{r.t2} の入力が見えない")
        if missing or thin:
            why = ("（" + "・".join(missing) + " が FlyWire に無い）") if missing else "（" + "・".join(thin) + "）"
            rows.append({**base, "n_sides": 0, "syn_cross_fwd": -1, "syn_cross_rev": -1,
                         "cross_total": -1, "detail": why, "total_syn": src_fwd + src_rev})
            continue
        f_tot = int(total.get((r.t1, r.t2), 0))
        b_tot = int(total.get((r.t2, r.t1), 0))
        det = [f"{s}:{int(fwd.get((r.t1, r.t2, s), 0))}/{int(fwd.get((r.t2, r.t1, s), 0))}"
               for s in ("R", "L", "other") if fwd.get((r.t1, r.t2, s), 0) or fwd.get((r.t2, r.t1, s), 0)]
        rows.append({**base, "n_sides": 2, "syn_cross_fwd": f_tot, "syn_cross_rev": b_tot,
                     "cross_total": f_tot + b_tot,
                     "detail": " ".join(det) if det else "合計 0/0",
                     "total_syn": src_fwd + src_rev + f_tot + b_tot})
        print(f"{r.t1:>14s} → {r.t2:<14s} 元 {src_fwd}/{src_rev} | FlyWire 合計 {f_tot}/{b_tot}"
              + ("  ★" if rows[-1]["total_syn"] == 0 else ""))
    out = pd.DataFrame(rows)
    f = res / f"cross_{TAG}.csv"
    out.to_csv(f, index=False, encoding="utf-8")
    ok = out[(out.total_syn == 0) & (out.n_sides > 0)]
    print(f"{len(out)} 組（照合できたのは {int((out.n_sides > 0).sum())} 組）、合計 0 の組 {len(ok)} 組 -> {f}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true", help="大きいファイルからタイプ単位の集計表を作る")
    ap.add_argument("--connections", default=None, help="proofread_connections_783.feather のパス")
    ap.add_argument("--annotations", default=None, help="Supplemental_file1_neuron_annotations.tsv のパス")
    ap.add_argument("--data", nargs="*", default=[], help="元データ（hemibrain）の pickle")
    ap.add_argument("--top", type=int, default=150)
    ap.add_argument("--min-visible", type=int, default=100)
    ap.add_argument("--pairs", default=None, help="調べる組の CSV（candidates_ranked と同じ列）")
    ap.add_argument("--tag", default=None, help="出力先 results/<tag>/")
    ap.add_argument("--min-expected", type=float, default=8.0,
                    help="元データでの期待シナプス数の下限（validate_cross.py と揃える）")
    args = ap.parse_args()

    if args.build:
        if not (args.connections and args.annotations):
            ap.error("--build には --connections と --annotations が必要です")
        build(args.connections, args.annotations)
    if args.pairs:
        if not args.data:
            ap.error("--pairs には --data（逆向きを数える元データ）が必要です")
        print(f"\n----- {args.pairs}")
        check(args.data[0], args.top, args.min_visible, args.min_expected,
              pairs=args.pairs, extra_data=args.data[1:], tag=args.tag)
    else:
        for d in args.data:
            print(f"\n----- {d}")
            check(d, args.top, args.min_visible, args.min_expected)


if __name__ == "__main__":
    main()
