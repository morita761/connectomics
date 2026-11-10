"""まだ探索していない細胞タイプ・領域を洗い出す。

  python coverage.py

これまでの探索は領域ごとにデータを取ってきた（results/hemibrain_<領域>/）。
hemibrain 全体のうち、どの細胞型がまだ 1 度も候補の母集団に入っていないかを出す。

やること:
  1. results/*/type_pairs_output.csv に出てくる型を「探索済み」として集める
  2. neuPrint から hemibrain の Traced なニューロンを全部取り、型ごとの
     主な領域（シナプスが最も多い ROI）と細胞数を出す
  3. 探索済みと突き合わせ、領域ごとに「未探索の型」を数える

出力: results/coverage.csv（型ごと）と results/coverage_roi.csv（領域ごと）
"""
from __future__ import annotations

import argparse
import json
from collections import Counter

import pandas as pd

from common.runner import ROOT

RESULTS = ROOT / "results"
# neuPrint の roiInfo には親領域と子領域の両方が入っている。親は数えない。
SUPER_ROIS = {"SNP", "VLNP", "VMNP", "INP", "LX", "PENP", "OL", "MB", "CX",
              "MB(+ACA)", "VX", "SLP+SIP", "GNG+", ""}


def screened_types() -> tuple[set[str], dict[str, set[str]]]:
    """これまでの探索に出てきた型と、どの領域の探索で出てきたか。"""
    seen: set[str] = set()
    where: dict[str, set[str]] = {}
    for f in sorted(RESULTS.glob("*/type_pairs_output.csv")):
        region = f.parent.name
        try:
            d = pd.read_csv(f, usecols=["t1", "t2"])
        except Exception:
            continue
        ts = set(d.t1.dropna()) | set(d.t2.dropna())
        seen |= ts
        for t in ts:
            where.setdefault(t, set()).add(region)
    return seen, where


def all_types(dataset: str) -> pd.DataFrame:
    """hemibrain の Traced なニューロンを型ごとにまとめ、主な領域を決める。"""
    from neuprint import NeuronCriteria as NC, fetch_neurons

    from common.loaders import _retry, neuprint_client
    client = neuprint_client(dataset)
    n, _ = _retry(fetch_neurons, NC(status="Traced"), client=client)
    rows = []
    for r in n.itertuples(index=False):
        if not isinstance(r.type, str) or not r.type:
            continue
        ri = r.roiInfo if isinstance(r.roiInfo, dict) else json.loads(r.roiInfo or "{}")
        # 親領域（SNP = SMP+SIP+SLP など）は細かい領域と重複して数えられているので除く
        best, best_n = "", -1
        for k, v in ri.items():
            if k.split("(")[0] in SUPER_ROIS:
                continue
            tot = int(v.get("pre", 0)) + int(v.get("post", 0))
            if tot > best_n:
                best, best_n = k, tot
        rows.append({"type": r.type, "bodyId": r.bodyId, "roi": best,
                     "pre": int(r.pre or 0), "post": int(r.post or 0)})
    d = pd.DataFrame(rows)
    g = d.groupby("type").agg(n_cells=("bodyId", "size"), pre=("pre", "sum"),
                              post=("post", "sum"),
                              roi=("roi", lambda s: Counter(s).most_common(1)[0][0]))
    return g.reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="hemibrain:v1.2.1")
    ap.add_argument("--min-cells", type=int, default=1)
    args = ap.parse_args()

    seen, where = screened_types()
    print(f"これまでの探索に出てきた型: {len(seen)}")

    d = all_types(args.dataset)
    print(f"{args.dataset} の Traced な型: {len(d)}（細胞 {int(d.n_cells.sum())}）")
    d["screened"] = d.type.isin(seen)
    d["regions"] = d.type.map(lambda t: " ".join(sorted(where.get(t, ()))))
    d = d.sort_values(["screened", "n_cells"], ascending=[True, False])
    f = RESULTS / "coverage.csv"
    d.to_csv(f, index=False, encoding="utf-8")

    n_un = int((~d.screened).sum())
    print(f"未探索の型: {n_un} / {len(d)}（細胞 {int(d[~d.screened].n_cells.sum())}）\n")

    roi = d.groupby("roi").agg(型=("type", "size"), 未探索=("screened", lambda s: int((~s).sum()),),
                               細胞=("n_cells", "sum"),
                               未探索の細胞=("n_cells", lambda s: 0)).reset_index()
    un = d[~d.screened].groupby("roi").n_cells.sum()
    roi["未探索の細胞"] = roi.roi.map(un).fillna(0).astype(int)
    roi["未探索率"] = (roi.未探索 / roi.型 * 100).round(0)
    roi = roi.sort_values("未探索", ascending=False)
    roi.to_csv(RESULTS / "coverage_roi.csv", index=False, encoding="utf-8")

    print("領域ごとの未探索（型の数が多い順）")
    print(roi.head(25).to_string(index=False))
    print(f"\n-> {f} / {RESULTS / 'coverage_roi.csv'}")


if __name__ == "__main__":
    main()
