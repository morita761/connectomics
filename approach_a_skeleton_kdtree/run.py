"""アプローチ A: スケルトン点群 + KD-tree による近接検出（最も手軽なベースライン）

考え方:
  - 全ニューロンのスケルトンを 0.5 µm 間隔の点に分割し、1 本の KD-tree に入れる
  - 中心線どうしが skeleton_contact_um (1 µm) 以内にある点ペアを全列挙
  - 細胞ペアごとに「相手の近くを通るケーブル長」= 接触ケーブル長 [µm] を集計
  - 接触が min_contact_cable_um 以上 かつ 両方向合計シナプス数 <= weak_synapse_max を候補とする

長所: 速い・スケルトンだけで済む・数千ニューロンまで現実的。
短所: 太さを無視（半径はチューブ表面の B で考慮）、軸索/樹状突起の区別なし（→ C）、
      接触量に対して何シナプスが「普通」かの基準がない（→ D）。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from common.config import Thresholds
from common.geometry import cable_points, close_point_pairs, contacts, stack
from common.model import Dataset, attach_weights

UM = 1000.0


def find_candidates(ds: Dataset, th: Thresholds, **_) -> pd.DataFrame:
    cps = {i: cable_points(n, th.cable_spacing_um * UM, th.default_radius_nm)
           for i, n in ds.neurons.items()}
    cloud = stack(cps)
    df = contacts(cloud, th.skeleton_contact_um * UM, None, "a", "b")
    df["contact_um"] = (df.contact_a + df.contact_b) / 2 / UM
    df = attach_weights(df, ds, "a", "b", directed=False)
    df["score"] = df.contact_um / (1 + df.weight)
    df["flagged"] = (df.contact_um >= th.min_contact_cable_um) & (df.weight <= th.weak_synapse_max)
    cols = ["a", "b", "contact_um", "weight", "weight_ab", "weight_ba", "score", "flagged",
            "n_point_pairs", "cx", "cy", "cz"]
    return df[cols].sort_values("score", ascending=False).reset_index(drop=True)


if __name__ == "__main__":
    from common.runner import cli
    cli("approach_a_skeleton_kdtree", find_candidates)
