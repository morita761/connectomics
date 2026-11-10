"""アプローチ F: 接触パッチごとの解析（ペアの中で「シナプスを避ける部位」を探す）

考え方:
  - A–E はペア単位で「近いのに結合が弱い」を探す。しかし、全体としては結合しているペアでも、
    特定の部位（ある枝どうしの接触）だけシナプスがない、ということがありうる。
  - A の軸索 × B の樹状突起 の接触点を、空間的につながったまとまり =「接触パッチ」に分ける
    （中点どうしが patch_link_um 以内なら同じパッチ、単連結クラスタリング）。
  - 各シナプスを最寄りのパッチに割り当て、パッチごとのシナプス数を数える。
  - そのパッチの期待シナプス数 = 同じペアの「他のパッチ」のシナプス密度 × パッチの接触長
    （他のパッチが小さいと推定が不安定なので、全体平均の密度に向けて縮小推定する）
  - Poisson で P(X <= 観測値 | 期待値) を計算し、FDR 補正。

出力は 1 行 = 1 パッチ。ペア単位で比べるときは、そのペアの最もスコアの高いパッチを使う。
必要なもの: シナプスごとの (pre, post, 位置)。Dataset.synapses が無い場合は、A の pre 部位のうち
  B の post 部位が近くにあるものを A→B のシナプスとみなす（近似）。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

from common.config import Thresholds
from common.geometry import contact_patches
from common.model import Dataset
from common.runner import load_approach

UM = 1000.0
axon_dendrite_pairs = load_approach("approach_c_axon_dendrite").axon_dendrite_pairs
bh_qvalues = load_approach("approach_d_peters_null_model").bh_qvalues


def pair_synapses(ds: Dataset, pre: int, post: int, tol_nm: float = 1000.0) -> np.ndarray:
    if ds.synapses is not None:
        s = ds.synapses
        return s.loc[(s.pre == pre) & (s.post == post), ["x", "y", "z"]].to_numpy(float)
    a, b = ds.neurons[pre], ds.neurons[post]
    if a.n_pre == 0 or b.n_post == 0:
        return np.empty((0, 3))
    d, _ = cKDTree(b.post).query(a.pre, distance_upper_bound=tol_nm)
    return a.pre[np.isfinite(d)]


def find_candidates(ds: Dataset, th: Thresholds, compartments: str = "auto", **_) -> pd.DataFrame:
    cloud, i, j = axon_dendrite_pairs(ds, th, compartments)
    patch = contact_patches(i, j, cloud, th.patch_link_um * UM)
    w, P, lab = cloud["weight"], cloud["points"], cloud["label"]
    pp = pd.DataFrame({"pre": lab[i], "post": lab[j], "i": i, "j": j, "patch": patch})

    # パッチごとの接触長と重心
    la = pp.drop_duplicates(["patch", "i"]).assign(v=lambda d: w[d.i]).groupby("patch").v.sum()
    lb = pp.drop_duplicates(["patch", "j"]).assign(v=lambda d: w[d.j]).groupby("patch").v.sum()
    mid = (P[i] + P[j]) / 2
    pp[["mx", "my", "mz"]] = mid
    pt = pp.groupby("patch").agg(pre=("pre", "first"), post=("post", "first"),
                                 cx=("mx", "mean"), cy=("my", "mean"), cz=("mz", "mean"))
    pt["contact_um"] = (la + lb) / 2 / UM
    pt["n_syn"] = 0

    # シナプスを最寄りのパッチに割り当てる
    r = th.patch_synapse_radius_um * UM
    for (a, b), grp in pp.groupby(["pre", "post"]):
        syn = pair_synapses(ds, a, b)
        if len(syn) == 0:
            continue
        d, k = cKDTree(grp[["mx", "my", "mz"]].to_numpy()).query(syn, distance_upper_bound=r)
        hit = np.isfinite(d)
        if hit.any():
            cnt = pd.Series(grp.patch.to_numpy()[k[hit]]).value_counts()
            pt.loc[cnt.index, "n_syn"] += cnt.to_numpy()

    pt = pt.reset_index()
    pair_L = pt.groupby(["pre", "post"]).contact_um.transform("sum")
    pair_S = pt.groupby(["pre", "post"]).n_syn.transform("sum")
    wd = ds.directed_weights()
    pt["pair_weight"] = [wd.get((int(a), int(b)), 0) for a, b in zip(pt.pre, pt.post)]
    pt["pair_contact_um"] = pair_L
    pt["n_patches_in_pair"] = pt.groupby(["pre", "post"]).patch.transform("size")

    # 期待値: 同じペアの他のパッチの密度（全体平均へ縮小）
    global_rate = pt.n_syn.sum() / pt.contact_um.sum()
    a0 = th.patch_prior_synapses
    other_S = pair_S - pt.n_syn
    other_L = pair_L - pt.contact_um
    rate = (other_S + a0) / (other_L + a0 / global_rate)
    pt["rate_other_per_um"] = rate
    pt["expected"] = rate * pt.contact_um
    pt["p_lower"] = stats.poisson.cdf(pt.n_syn, pt.expected)

    big = pt.contact_um >= th.min_patch_contact_um
    pt["q"] = 1.0
    pt.loc[big, "q"] = bh_qvalues(pt.loc[big, "p_lower"].to_numpy())
    pt["score"] = np.where(big, -np.log10(np.maximum(pt.p_lower, 1e-300)), 0.0)
    pt["flagged"] = (big & (pt.q < th.fdr_alpha) & (pt.expected >= th.min_patch_expected)
                     & (pt.n_syn == 0))
    pt = pt.rename(columns={"patch": "patch_id"})
    cols = ["pre", "post", "patch_id", "contact_um", "n_syn", "expected", "rate_other_per_um",
            "pair_weight", "pair_contact_um", "n_patches_in_pair", "p_lower", "q", "score", "flagged",
            "cx", "cy", "cz"]
    return pt[cols].sort_values("score", ascending=False).reset_index(drop=True)


def _args(ap):
    ap.add_argument("--compartments", default="auto", choices=["auto", "given", "navis", "synapse_vote"])


if __name__ == "__main__":
    from common.runner import cli
    cli("approach_f_contact_patches", find_candidates, _args)
