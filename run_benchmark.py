"""全アプローチを正解つき合成データで走らせて比較する。

  python run_benchmark.py                       # 既定: 80 ニューロン, seed 0..2
  python run_benchmark.py --n-neurons 150 --seeds 0 1 2 3 4
  python run_benchmark.py --compartments given  # 正解の軸索/樹状突起を使う（C, D の上限性能）
  python run_benchmark.py --set skeleton_contact_um=2.0

approach_*/run.py に find_candidates(ds, th, **kw) を持つフォルダを追加すれば、自動で比較対象に入る。
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from common.config import Thresholds, parse_overrides
from common.evaluate import evaluate, evaluate_patches
from common.runner import ROOT, list_approaches, load_approach
from common.synthetic import make_dataset


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-neurons", type=int, default=80)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--compartments", default="navis",
                    choices=["synapse_vote", "navis", "given", "auto"],
                    help="C, D で使う軸索/樹状突起の推定法（実データと同じ条件にするなら synapse_vote か navis）")
    ap.add_argument("--approaches", nargs="*", default=None)
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()

    th = Thresholds().update(parse_overrides(args.set))
    names = args.approaches or list_approaches()
    res_dir = ROOT / "results"
    res_dir.mkdir(exist_ok=True)

    rows = []
    for seed in args.seeds:
        ds = make_dataset(n_neurons=args.n_neurons, seed=seed, th=th)
        ds.save(ROOT / "data" / f"{ds.name}.pkl")
        if args.compartments != "given":
            for n in ds.neurons.values():
                n.compartment = None  # 正解ラベルを隠す（実データと同じ条件）
        n_pt = len(ds.meta.get("patch_truth", []))
        print(f"\n=== {ds.name}: {len(ds.edges)} edges, {len(ds.truth)} planted pairs, "
              f"{n_pt} avoided patches ===")
        for name in names:
            mod = load_approach(name)
            t0 = time.perf_counter()
            df = mod.find_candidates(ds, th, compartments=args.compartments)
            dt = time.perf_counter() - t0
            df.to_csv(res_dir / f"{name}__{ds.name}.csv", index=False)
            m = evaluate(df, ds.truth)
            m.update(evaluate_patches(df, ds.meta.get("patch_truth")))
            m.update(approach=name, seed=seed, seconds=dt, n_pairs=len(df))
            rows.append(m)
            extra = f"  patch AP={m['AP_patch']:.3f}" if "AP_patch" in m else ""
            print(f"  {name:32s} AP={m['AP_undirected']:.3f}  P@K={m['P@K_undirected']:.3f}  "
                  f"flag P/R={m['precision_flagged']:.2f}/{m['recall_flagged']:.2f}  {dt:.1f}s{extra}")

    res = pd.DataFrame(rows)
    res.to_csv(res_dir / "benchmark_runs.csv", index=False)
    cols = ["AP_undirected", "P@K_undirected", "AP_directed", "P@K_directed",
            "precision_flagged", "recall_flagged", "flagged", "seconds"]
    cols += [c for c in ("AP_patch", "P@K_patch", "precision_flagged_patch", "recall_flagged_patch")
             if c in res]
    summary = res.groupby("approach")[cols].mean().sort_values("AP_undirected", ascending=False)
    md = summary.round(3).to_markdown() if _has_tabulate() else summary.round(3).to_string()
    text = (f"# Benchmark\n\nn_neurons={args.n_neurons}, seeds={args.seeds}, "
            f"compartments={args.compartments}, overrides={args.set}\n\n{md}\n")
    (res_dir / "benchmark.md").write_text(text, encoding="utf-8")
    print("\n" + text)


def _has_tabulate() -> bool:
    try:
        import tabulate  # noqa: F401
        return True
    except ImportError:
        return False


if __name__ == "__main__":
    main()
