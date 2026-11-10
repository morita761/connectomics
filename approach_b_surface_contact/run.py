"""アプローチ B: 膜表面どうしの接触面積（メッシュ間距離ベース）

考え方:
  - メッシュがあればその頂点（面積の重み付き）、なければスケルトン半径からチューブ表面を作る
  - 表面点どうしが surface_contact_um (0.2 µm) 以内 = 膜が接している
  - 細胞ペアごとの接触面積 [µm²] を計算し、面積が大きいのに結合が弱いペアを候補とする

A との違い: 中心線ではなく実際の膜の位置を見るので、太い突起どうしの「見かけの近接」や
  細い突起どうしの見落としが減る。シナプスは膜接触面に作られるので、接触面積は
  「シナプスを作る機会」の指標としてケーブル長より直接的。
短所: 点数が多く重い（メッシュは頂点数が多いので、ROI を絞るか簡略化メッシュを使う）。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from common.config import Thresholds
from common.geometry import cable_points, contacts, stack, surface_points
from common.model import Dataset, attach_weights

UM = 1000.0


def find_candidates(ds: Dataset, th: Thresholds, **_) -> pd.DataFrame:
    surf = {}
    for i, n in ds.neurons.items():
        cp = cable_points(n, th.cable_spacing_um * UM, th.default_radius_nm)
        surf[i] = surface_points(n, cp, th.surface_points_per_ring)
    cloud = stack(surf)
    df = contacts(cloud, th.surface_contact_um * UM, None, "a", "b")
    df["contact_area_um2"] = (df.contact_a + df.contact_b) / 2 / UM**2
    df = attach_weights(df, ds, "a", "b", directed=False)
    df["score"] = df.contact_area_um2 / (1 + df.weight)
    df["flagged"] = (df.contact_area_um2 >= th.min_contact_area_um2) & (df.weight <= th.weak_synapse_max)
    cols = ["a", "b", "contact_area_um2", "weight", "weight_ab", "weight_ba", "score", "flagged",
            "n_point_pairs", "cx", "cy", "cz"]
    return df[cols].sort_values("score", ascending=False).reset_index(drop=True)


if __name__ == "__main__":
    from common.runner import cli
    cli("approach_b_surface_contact", find_candidates)
