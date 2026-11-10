"""アプローチ C: 軸索 → 樹状突起の局所近傍判定（有向）

考え方:
  - シナプスはふつう「A の軸索（出力部位）」と「B の樹状突起（入力部位）」の間にできる。
    軸索どうし・樹状突起どうしが近いだけではシナプスの機会にならない。
  - 各ニューロンを軸索/樹状突起に分け（navis の synapse flow centrality、なければ
    近傍シナプスの pre/post 多数決）、A の軸索 × B の樹状突起 の近接だけを数える。
  - 有向ペア A→B ごとに接触ケーブル長を出し、A→B のシナプス数が弱いものを候補とする。

A との違い: 方向性がある。B→A は強いが A→B はゼロ、という非対称なケースも拾える。
注意: ハエのニューロンは軸索/樹状突起が混在するものも多い（局所介在ニューロンなど）。
      そういう細胞は分割が不安定なので、結果を見るときは compartment 推定の質に注意。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from common.compartments import label_cable
from common.config import Thresholds
from common.geometry import cable_points, close_point_pairs, contacts, stack
from common.model import Dataset, attach_weights

UM = 1000.0


def labeled_cloud(ds: Dataset, th: Thresholds, compartments: str = "auto") -> dict:
    """全ニューロンのケーブル点群（軸索/樹状突起ラベル付き）。"""
    cps = {}
    for i, n in ds.neurons.items():
        cp = cable_points(n, th.cable_spacing_um * UM, th.default_radius_nm)
        cp["compartment"] = label_cable(n, cp, th, compartments)
        cps[i] = cp
    return stack(cps)


def axon_dendrite_pairs(ds: Dataset, th: Thresholds, compartments: str = "auto"):
    """軸索の点 i と樹状突起の点 j の近接ペア（点単位）と、その点群を返す。F からも使う。"""
    cloud = labeled_cloud(ds, th, compartments)
    i, j = close_point_pairs(cloud["points"], cloud["label"], th.skeleton_contact_um * UM,
                             roles=cloud["compartment"])
    return cloud, i, j


def axon_dendrite_contacts(ds: Dataset, th: Thresholds, compartments: str = "auto") -> pd.DataFrame:
    """有向ペア (pre, post) ごとの 軸索–樹状突起 接触ケーブル長 [µm]。D, E からも使う。"""
    cloud = labeled_cloud(ds, th, compartments)
    df = contacts(cloud, th.skeleton_contact_um * UM, cloud["compartment"], "pre", "post")
    df = df.rename(columns={"contact_a": "axon_contact", "contact_b": "dendrite_contact"})
    df["axon_contact_um"] = df.pop("axon_contact") / UM
    df["dendrite_contact_um"] = df.pop("dendrite_contact") / UM
    df["contact_um"] = (df.axon_contact_um + df.dendrite_contact_um) / 2
    return df


def find_candidates(ds: Dataset, th: Thresholds, compartments: str = "auto", **_) -> pd.DataFrame:
    df = axon_dendrite_contacts(ds, th, compartments)
    df = attach_weights(df, ds, "pre", "post", directed=True)
    df["score"] = df.contact_um / (1 + df.weight)
    df["flagged"] = (df.contact_um >= th.min_contact_cable_um) & (df.weight <= th.weak_synapse_max)
    cols = ["pre", "post", "contact_um", "axon_contact_um", "dendrite_contact_um", "weight",
            "weight_reverse", "score", "flagged", "n_point_pairs", "cx", "cy", "cz"]
    return df[cols].sort_values("score", ascending=False).reset_index(drop=True)


def _args(ap):
    ap.add_argument("--compartments", default="auto", choices=["auto", "given", "navis", "synapse_vote"])


if __name__ == "__main__":
    from common.runner import cli
    cli("approach_c_axon_dendrite", find_candidates, _args)
