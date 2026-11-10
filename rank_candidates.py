"""上流→下流の「接触のわりに結合しない」タイプの組を、比較対象と合わせて並べ直す。

  python rank_candidates.py --data data/hemibrain_lh.pkl

各候補（type_pairs_output.csv の flagged かつ feedforward）について:
  - comp_type / comp_rate: 同じ下流タイプに対して、同じグループの他の上流タイプのうち、
    接触 --min-comp-contact µm 以上で最も密に結合しているものと、その密度（シナプス/µm）
    → 「下流側がそもそも入力を受けない場所」ではなく、「他の上流からは受けるのに、この上流からは受けない」ことの裏付け
  - self_rate: この組の密度（シナプス/µm）
  - ratio: comp_rate / self_rate（self_rate は下限 0.002 で割る）
  - 両側のドライバー（driver_matches.csv があれば）
並べ方: 期待値、比較対象との差、接触しているペア数を合わせた score_rank。
出力: results/<データ名>/candidates_ranked.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common.runner import ROOT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--min-comp-contact", type=float, default=20.0)
    ap.add_argument("--min-pairs", type=int, default=1, help="接触しているニューロンのペア数の下限")
    ap.add_argument("--top", type=int, default=40)
    args = ap.parse_args()

    res = ROOT / "results" / Path(args.data).stem
    tp = pd.read_csv(res / "type_pairs_output.csv")
    tp["rate"] = tp.synapses / tp.contact_total_um
    cand = tp[tp.flagged & tp.feedforward & (tp.n_pairs_contact >= args.min_pairs)].copy()

    pool = tp[tp.contact_total_um >= args.min_comp_contact]
    best = (pool.sort_values("rate", ascending=False)
            .groupby(["t2", "g1"])[["t1", "rate", "contact_total_um"]].apply(lambda d: d.head(3).values.tolist()))

    comp_t, comp_r, comp_c = [], [], []
    for r in cand.itertuples(index=False):
        rows = [x for x in best.get((r.t2, r.g1), []) if x[0] != r.t1]
        if rows:
            comp_t.append(rows[0][0]); comp_r.append(rows[0][1]); comp_c.append(rows[0][2])
        else:
            comp_t.append(None); comp_r.append(np.nan); comp_c.append(np.nan)
    cand["comp_type"], cand["comp_rate"], cand["comp_contact_um"] = comp_t, comp_r, comp_c
    cand["self_rate"] = cand.rate
    cand["ratio"] = cand.comp_rate / cand.self_rate.clip(lower=0.002)  # self が 0 でも割れるように下限 0.002
    cand["score_rank"] = (np.log10(cand.expected) + np.log10(cand.ratio.fillna(1).clip(lower=1))
                          + 0.3 * np.log2(cand.n_pairs_contact))

    dm = res / "driver_matches.csv"
    if dm.exists():
        d = pd.read_csv(dm)[["t1", "t2", "availability", "up_gal4", "down_lexa", "down_gal4", "up_lexa"]]
        cand = cand.merge(d, on=["t1", "t2"], how="left")
    cand = cand.sort_values("score_rank", ascending=False)
    cols = ["t1", "t2", "g1", "g2", "n_cells_1", "n_cells_2", "n_pairs_contact", "contact_total_um", "synapses",
            "expected", "self_rate", "comp_type", "comp_rate", "ratio", "score_rank", "ex_a", "ex_b", "ex_weight"]
    cols += [c for c in ("availability", "up_gal4", "down_lexa", "down_gal4", "up_lexa") if c in cand]
    cand[cols].to_csv(res / "candidates_ranked.csv", index=False, encoding="utf-8")
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 30):
        print(cand[cols[:15]].head(args.top).round(3).to_string(index=False))
    print(f"{len(cand)} 組 -> {res / 'candidates_ranked.csv'}")


if __name__ == "__main__":
    main()
