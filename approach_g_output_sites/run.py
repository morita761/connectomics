"""アプローチ G: シナプス前の「出力部位」から相手のどこへでも（有向）

考え方:
  - C は「A の軸索 × B の樹状突起」だけを見るが、キノコ体では
      DAN → KC の軸索（軸索–軸索）、KC → KC（軸索–軸索）、APL → KC など
    樹状突起以外が受け手になる結合が多い。また navis の軸索/樹状突起の分割は、入力も多い
    KC の軸索などで不安定。
  - そこで分割を使わず、A の「出力部位」= A 自身の T-bar（シナプス前部位）から
    output_site_radius_um 以内のケーブル、と定義する。シナプス A→B は必ずここから出る。
  - 有向ペア A→B ごとに「A の出力部位が B のどこかに skeleton_contact_um 以内で接する長さ」を数え、
    A→B のシナプス数と比べる。

長所: 受け手の区画を仮定しない。上流（出力側）→ 下流（受け手）の向きがはっきりする。
短所: T-bar の少ない細胞（再構築の途中など）は出力部位が短くなる。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from common.config import Thresholds
from common.geometry import cable_points, contacts, stack
from common.model import AXON, UNKNOWN, Dataset, attach_weights

UM = 1000.0


def output_cloud(ds: Dataset, th: Thresholds) -> dict:
    cps = {}
    for i, n in ds.neurons.items():
        cp = cable_points(n, th.cable_spacing_um * UM, th.default_radius_nm)
        role = np.full(len(cp["points"]), UNKNOWN, dtype=np.int8)
        if n.n_pre:
            d, _ = cKDTree(n.pre).query(cp["points"], distance_upper_bound=th.output_site_radius_um * UM)
            role[np.isfinite(d)] = AXON  # ここでは AXON = 出力部位 の意味で使う
        cp["compartment"] = role
        cps[i] = cp
    return stack(cps)


def output_contacts(ds: Dataset, th: Thresholds) -> pd.DataFrame:
    cloud = output_cloud(ds, th)
    df = contacts(cloud, th.skeleton_contact_um * UM, cloud["compartment"], "pre", "post", post_any=True)
    df["output_contact_um"] = df.pop("contact_a") / UM   # 上流側の出力部位のうち、相手に接する長さ
    df["target_contact_um"] = df.pop("contact_b") / UM   # 下流側のうち、相手の出力部位に接する長さ
    df["contact_um"] = (df.output_contact_um + df.target_contact_um) / 2
    return df


def find_candidates(ds: Dataset, th: Thresholds, **_) -> pd.DataFrame:
    df = output_contacts(ds, th)
    df = attach_weights(df, ds, "pre", "post", directed=True)
    df["score"] = df.contact_um / (1 + df.weight)
    df["flagged"] = (df.contact_um >= th.min_contact_cable_um) & (df.weight <= th.weak_synapse_max)
    cols = ["pre", "post", "contact_um", "output_contact_um", "target_contact_um", "weight",
            "weight_reverse", "score", "flagged", "n_point_pairs", "cx", "cy", "cz"]
    return df[cols].sort_values("score", ascending=False).reset_index(drop=True)


if __name__ == "__main__":
    from common.runner import cli
    cli("approach_g_output_sites", find_candidates)
