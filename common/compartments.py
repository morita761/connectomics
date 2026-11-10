"""軸索 / 樹状突起の推定（アプローチ C, D で使用）。

method:
  "given"         : Neuron.compartment をそのまま使う（合成データの正解など）
  "navis"         : navis.split_axon_dendrite（synapse flow centrality, Schneider-Mizell et al. 2016）
  "synapse_vote"  : 各点の近くのシナプスが pre / post のどちらが多いかで多数決（依存なしの簡易版）
  "auto"          : given → navis → synapse_vote の順に使えるものを使う
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .config import Thresholds
from .model import AXON, DENDRITE, UNKNOWN, Neuron

MAX_NAVIS_SYNAPSES = 15_000
MAX_NAVIS_NODES = 40_000


def synapse_vote(n: Neuron, points: np.ndarray, th: Thresholds) -> np.ndarray:
    pre = n.pre if n.pre is not None else np.empty((0, 3))
    post = n.post if n.post is not None else np.empty((0, 3))
    syn = np.vstack([pre, post])
    if len(syn) == 0:
        return np.zeros(len(points), dtype=np.int8)
    is_pre = np.r_[np.ones(len(pre)), np.zeros(len(post))]
    k = min(th.compartment_k, len(syn))
    d, nn = cKDTree(syn).query(points, k=k, distance_upper_bound=th.compartment_radius_um * 1000)
    d, nn = d.reshape(len(points), -1), nn.reshape(len(points), -1)
    valid = np.isfinite(d)
    votes = np.where(valid, is_pre[np.minimum(nn, len(syn) - 1)], 0.0)
    cnt = valid.sum(1)
    frac = np.divide(votes.sum(1), cnt, out=np.full(len(points), 0.5), where=cnt > 0)
    comp = np.full(len(points), UNKNOWN, dtype=np.int8)
    comp[(cnt > 0) & (frac > 0.5)] = AXON
    comp[(cnt > 0) & (frac < 0.5)] = DENDRITE
    # シナプスが近くにない点は、最も近いラベル付きの点のラベルを引き継ぐ
    lab = comp != UNKNOWN
    if lab.any() and (~lab).any():
        _, nn2 = cKDTree(points[lab]).query(points[~lab])
        comp[~lab] = comp[lab][nn2]
    return comp


def navis_split(n: Neuron) -> np.ndarray | None:
    """ノード単位のラベルを返す（失敗したら None）。"""
    try:
        import navis
        import pandas as pd
    except ImportError:
        return None
    if n.n_pre == 0 or n.n_post == 0:
        return None
    # 巨大なニューロン（APL・MBON など、シナプス数万）は navis の flow centrality が
    # 100 GB 級のメモリを使うことがあるので、分割しない（呼び出し側で synapse_vote に切り替わる）
    if n.n_pre + n.n_post > MAX_NAVIS_SYNAPSES or len(n.nodes) > MAX_NAVIS_NODES:
        return None
    nodes = pd.DataFrame({"node_id": np.arange(len(n.nodes)), "parent_id": n.parents,
                          "x": n.nodes[:, 0], "y": n.nodes[:, 1], "z": n.nodes[:, 2]})
    tn = navis.TreeNeuron(nodes, units="nm")
    syn = np.vstack([n.pre, n.post])
    _, nid = cKDTree(n.nodes).query(syn)
    tn.connectors = pd.DataFrame({"connector_id": np.arange(len(syn)), "node_id": nid,
                                  "x": syn[:, 0], "y": syn[:, 1], "z": syn[:, 2],
                                  "type": ["pre"] * n.n_pre + ["post"] * n.n_post})
    try:
        out = navis.split_axon_dendrite(tn, label_only=True)
    except Exception:
        return None
    lab = out.nodes.set_index("node_id").compartment.astype(object).reindex(np.arange(len(n.nodes)))
    return np.select([lab == "axon", lab == "dendrite"], [AXON, DENDRITE], UNKNOWN).astype(np.int8)


def label_cable(n: Neuron, cp: dict, th: Thresholds, method: str = "auto") -> np.ndarray:
    """cable_points() の各点に軸索/樹状突起ラベルを付ける。"""
    if method in ("given", "auto") and n.compartment is not None and np.any(n.compartment):
        return cp["compartment"]
    if method == "given":
        raise ValueError(f"neuron {n.id} has no given compartments")
    if method in ("navis", "auto"):
        node_lab = navis_split(n)
        if node_lab is not None and np.any(node_lab):
            _, nn = cKDTree(n.nodes).query(cp["points"])
            return node_lab[nn]
        # pre か post のどちらかしかない細胞などは navis で分割できない → 多数決にフォールバック
    return synapse_vote(n, cp["points"], th)
