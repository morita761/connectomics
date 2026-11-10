"""実データ (Dataset pickle) に全アプローチを適用し、候補ペアを 1 つの表にまとめる。

  python run_all.py --data data/hemibrain_mb.pkl
  python run_all.py --data data/hemibrain_mb.pkl --set skeleton_contact_um=1.5 weak_synapse_max=2

出力:
  results/<dataset>/<approach>.csv   各アプローチの全近接ペア
  results/<dataset>/consensus.csv    無向ペアごとに各アプローチの順位を並べ、平均順位でソート
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

from common.config import Thresholds, parse_overrides
from common.model import Dataset
from common.runner import ROOT, list_approaches, load_approach


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--approaches", nargs="*", default=None)
    ap.add_argument("--compartments", default="auto")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()

    th = Thresholds().update(parse_overrides(args.set))
    ds = Dataset.load(args.data)
    out = ROOT / "results" / Path(args.data).stem  # ds.name は ":" を含むことがある（Windows で不可）
    out.mkdir(parents=True, exist_ok=True)

    ranks = []
    for name in args.approaches or list_approaches():
        t0 = time.perf_counter()
        df = load_approach(name).find_candidates(ds, th, compartments=args.compartments)
        df.to_csv(out / f"{name}.csv", index=False)
        a, b = ("pre", "post") if "pre" in df else ("a", "b")
        key = df[[a, b]].min(axis=1).astype(str) + "-" + df[[a, b]].max(axis=1).astype(str)
        best = df.assign(key=key).groupby("key").score.max()
        ranks.append(best.rank(ascending=False).rename(name))
        print(f"{name:32s} {len(df):6d} pairs  {int(df.flagged.sum()):5d} flagged  "
              f"{time.perf_counter() - t0:.1f}s")

    cons = pd.concat(ranks, axis=1)
    cons["mean_rank"] = cons.mean(axis=1)
    cons = cons.sort_values("mean_rank")
    cons.to_csv(out / "consensus.csv")
    print(f"\nconsensus top 20 (lower rank = more suspicious):\n{cons.head(20)}")


if __name__ == "__main__":
    main()
