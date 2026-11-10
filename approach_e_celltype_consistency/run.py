"""アプローチ E: 細胞タイプの一貫性（同じタイプの組の他のペアと比べる）

考え方:
  - ハエのコネクトームでは、同じ細胞タイプのニューロンは配線もよく似ている
    （タイプ T1 → T2 の結合は、どの個体のペアでも同程度にある）。
  - そこで、ペア A→B の「期待シナプス数」を、同じタイプの組 (type(A), type(B)) に属する
    **A→B 以外のペア**から求める（leave-one-out）:
        rate_T = (タイプの組のシナプス総数 − w_AB) / (タイプの組の接触長総数 − L_AB)
        期待値 = rate_T × L_AB
  - 負の二項分布で P(W <= w_AB) を計算し、FDR 補正。

D との違い: D は全ペア共通の 1 本の回帰式で期待値を出すので、「もともとほとんど結合しない
  タイプの組」の自然な 0 まで怪しく見えてしまう。E はタイプの組ごとの結合のしやすさを基準にするので、
  「同じタイプの仲間はみんな結合しているのに、このペアだけ結合しない」という個別の逸脱を拾える。
必要なもの: Neuron.meta["type"]（hemibrain は neuPrint の type、FlyWire はアノテーション）。
  タイプが不明 or 同じタイプの組のペアが少ない場合は、全体の rate（D と同じ考え方）にフォールバック。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from common.config import Thresholds
from common.model import Dataset, attach_weights
from common.runner import load_approach

axon_dendrite_contacts = load_approach("approach_c_axon_dendrite").axon_dendrite_contacts
_d = load_approach("approach_d_peters_null_model")
nb_alpha, lower_tail_p, bh_qvalues = _d.nb_alpha, _d.lower_tail_p, _d.bh_qvalues


def find_candidates(ds: Dataset, th: Thresholds, compartments: str = "auto", type_key: str = "type",
                    **_) -> pd.DataFrame:
    df = axon_dendrite_contacts(ds, th, compartments)
    df = attach_weights(df, ds, "pre", "post", directed=True)
    df = df[df.contact_um > 0].reset_index(drop=True)

    # type_key="instance" にすると左右（_L/_R）なども区別したグループで比べる
    ntype = {i: (n.meta.get(type_key) or n.meta.get("type") or f"untyped_{i}") for i, n in ds.neurons.items()}
    df["type_pre"] = df.pre.map(ntype)
    df["type_post"] = df.post.map(ntype)

    g = df.groupby(["type_pre", "type_post"])
    W = g.weight.transform("sum")
    L = g.contact_um.transform("sum")
    N = g.weight.transform("size")
    strong = (df.weight > th.weak_synapse_max).astype(int)
    S = strong.groupby([df.type_pre, df.type_post]).transform("sum")

    w, l = df.weight.to_numpy(float), df.contact_um.to_numpy(float)
    global_rate = w.sum() / l.sum()
    n_other = (N - 1).to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        rate_loo = (W.to_numpy() - w) / (L.to_numpy() - l)
    use_type = (n_other >= th.min_type_pairs) & np.isfinite(rate_loo)
    rate = np.where(use_type, rate_loo, global_rate)
    mu = np.maximum(rate * l, 1e-9)

    # 過分散はタイプ単位の期待値が使えるペアで推定（外れ値を除いて 2 回）
    fit = use_type.copy() if use_type.sum() > 10 else np.ones(len(df), bool)
    alpha = nb_alpha(w[fit], mu[fit])
    p = lower_tail_p(w, mu, alpha)
    fit &= p >= th.robust_refit_p
    alpha = nb_alpha(w[fit], mu[fit])
    p = lower_tail_p(w, mu, alpha)

    df["n_type_pairs_other"] = n_other
    df["type_strong_frac_other"] = np.where(n_other > 0, (S.to_numpy() - strong.to_numpy()) / np.maximum(n_other, 1), np.nan)
    df["rate_source"] = np.where(use_type, "type", "global")
    df["expected"] = mu
    df["p_lower"] = p
    df["q"] = bh_qvalues(p)
    df["score"] = -np.log10(np.maximum(p, 1e-300))
    df["flagged"] = ((df.q < th.fdr_alpha) & (df.expected >= th.min_expected_synapses)
                     & (df.weight <= th.weak_synapse_max))
    print(f"  E: {use_type.sum()}/{len(df)} pairs use type-level rate, NB alpha={alpha:.3f}")
    cols = ["pre", "post", "type_pre", "type_post", "contact_um", "weight", "weight_reverse",
            "expected", "rate_source", "n_type_pairs_other", "type_strong_frac_other",
            "p_lower", "q", "score", "flagged", "cx", "cy", "cz"]
    return df[cols].sort_values("score", ascending=False).reset_index(drop=True)


def _args(ap):
    ap.add_argument("--compartments", default="auto", choices=["auto", "given", "navis", "synapse_vote"])
    ap.add_argument("--type-key", dest="type_key", default="type", choices=["type", "instance"],
                    help="グループ分けに使う名前（instance は左右なども区別する）")


if __name__ == "__main__":
    from common.runner import cli
    cli("approach_e_celltype_consistency", find_candidates, _args)
