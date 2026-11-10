"""複数のアプローチで一致した「近いのに結合しない」ペアを探す。

  python compare_methods.py --results results/hemibrain_mb --data data/hemibrain_mb.pkl
  python compare_methods.py --results results/hemibrain_mb --top 300 --min-methods 3

一致の基準は 2 つ:
  1. 上位一致: 各アプローチのスコア上位 --top 件（無向ペア単位）に、--min-methods 個以上のアプローチで入る
  2. flagged 一致: 各アプローチの閾値判定 (flagged) が --min-methods 個以上で True
     （D, E, F は FDR 補正が厳しくほとんど flagged にならないので、1 の方が実用的）

有向のアプローチ (C–F) は、A→B と B→A のうちスコアの高い方をそのペアの値とする。
F はパッチ単位なので、そのペアで最もスコアの高いパッチの値を使う。

出力 (results/<データ名>/):
  consensus_pairs.csv   一致したペアと、各アプローチでの順位・接触量・シナプス数
  consensus_summary.md  アプローチ間の重なりの表と、一致ペアの一覧
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common.model import Dataset

SHORT = {"approach_a_skeleton_kdtree": "A", "approach_b_surface_contact": "B",
         "approach_c_axon_dendrite": "C", "approach_d_peters_null_model": "D",
         "approach_e_celltype_consistency": "E", "approach_f_contact_patches": "F",
         "approach_g_output_sites": "G"}


def load_pairs(csv: Path) -> pd.DataFrame:
    df = pd.read_csv(csv)
    a, b = ("pre", "post") if "pre" in df else ("a", "b")
    lo, hi = np.minimum(df[a], df[b]), np.maximum(df[a], df[b])
    df = df.assign(n1=lo.astype("int64"), n2=hi.astype("int64"))
    df = df.sort_values("score", ascending=False).drop_duplicates(["n1", "n2"])
    df["rank"] = np.arange(1, len(df) + 1)
    if a == "pre":
        df["direction"] = df[a].astype(str) + "→" + df[b].astype(str)
    return df.set_index(["n1", "n2"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="run_all.py の出力フォルダ")
    ap.add_argument("--data", default=None, help="Dataset pickle（細胞タイプとシナプス数を付ける）")
    ap.add_argument("--top", type=int, default=200)
    ap.add_argument("--top-pct", type=float, default=None,
                    help="件数の代わりに各アプローチの上位何 %% を「上位」とするか（アプローチごとにペア数が違うため）")
    ap.add_argument("--min-methods", type=int, default=3)
    args = ap.parse_args()

    res = Path(args.results)
    tables = {}
    for name, s in SHORT.items():
        f = res / f"{name}.csv"
        if f.exists():
            tables[s] = load_pairs(f)
    methods = list(tables)
    print("methods:", methods)

    idx = sorted(set().union(*[set(t.index) for t in tables.values()]))
    out = pd.DataFrame(index=pd.MultiIndex.from_tuples(idx, names=["n1", "n2"]))
    for s, t in tables.items():
        out[f"rank_{s}"] = t["rank"]
        out[f"flag_{s}"] = t["flagged"].astype("boolean")
    cut = {s: (int(np.ceil(len(t) * args.top_pct / 100)) if args.top_pct else args.top)
           for s, t in tables.items()}
    top_label = f"{args.top_pct}%" if args.top_pct else f"{args.top} 件"
    out["n_top"] = sum((out[f"rank_{s}"] <= cut[s]).fillna(False).astype(int) for s in methods)
    out["top_methods"] = ["".join(s for s in methods if r[f"rank_{s}"] <= cut[s])
                          for _, r in out[[f"rank_{s}" for s in methods]].fillna(np.inf).iterrows()]
    out["n_flagged"] = sum(out[f"flag_{s}"].fillna(False).astype(int) for s in methods)
    out["mean_rank_pct"] = np.nanmean(np.column_stack(
        [out[f"rank_{s}"] / len(tables[s]) for s in methods]), axis=1) * 100

    # 補足情報
    for col, src in (("contact_um", "A"), ("contact_area_um2", "B"), ("expected", "D"), ("p_lower", "E")):
        if src in tables and col in tables[src]:
            out[f"{col}_{src}"] = tables[src][col]
    if "C" in tables:
        out["direction_C"] = tables["C"]["direction"]
    if args.data:
        ds = Dataset.load(args.data)
        w = ds.directed_weights()
        n1, n2 = out.index.get_level_values(0), out.index.get_level_values(1)
        out["type1"] = [ds.neurons[i].meta.get("type") for i in n1]
        out["type2"] = [ds.neurons[i].meta.get("type") for i in n2]
        out["w_12"] = [w.get((i, j), 0) for i, j in zip(n1, n2)]
        out["w_21"] = [w.get((j, i), 0) for i, j in zip(n1, n2)]

    cons = out[(out.n_top >= args.min_methods) | (out.n_flagged >= args.min_methods)]
    cons = cons.sort_values(["n_top", "mean_rank_pct"], ascending=[False, True])
    cons.reset_index().to_csv(res / "consensus_pairs.csv", index=False)

    # 上位 K 件どうしの重なり
    tops = {s: set(t.index[t["rank"] <= cut[s]]) for s, t in tables.items()}
    ov = pd.DataFrame([[len(tops[x] & tops[y]) for y in methods] for x in methods], index=methods, columns=methods)
    flags = {s: set(t.index[t["flagged"].astype(bool)]) for s, t in tables.items()}
    fv = pd.DataFrame([[len(flags[x] & flags[y]) for y in methods] for x in methods], index=methods, columns=methods)

    lines = [f"# 複数アプローチで一致したペア ({res.name})", "",
             f"上位 {top_label}どうしの重なり（対角 = 各アプローチの件数）:", "", ov.to_string(), "",
             "flagged どうしの重なり:", "", fv.to_string(), "",
             f"上位 {top_label}に入ったアプローチ数ごとのペア数:", "",
             out.n_top.value_counts().sort_index(ascending=False).to_string(), "",
             f"一致ペア（{args.min_methods} 個以上で上位 or flagged）: {len(cons)} 組", "",
             cons.head(30).to_string(), ""]
    text = "\n".join(lines)
    (res / "consensus_summary.md").write_text(text, encoding="utf-8")
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(text)


if __name__ == "__main__":
    main()
