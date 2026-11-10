"""アプローチ D: Peters' rule に基づく結合確率モデルとの比較（統計的な「結合の欠落」）

考え方:
  - Peters' rule: 軸索と樹状突起が出会う量（接触ケーブル長）に比例してシナプスができる、という帰無仮説。
  - データ全体から「接触量 → 期待シナプス数」の関係を回帰で学習する:
        log E[w_AB] = b0 + b1·log(接触長_AB) + b2·log(1+A の出力シナプス総数) + b3·log(1+B の入力シナプス総数)
    （出力/入力総数を入れるのは、もともとシナプスを多く作る細胞かどうかを補正するため）
  - 過分散を負の二項分布で扱い、各ペアについて p = P(W <= 観測値 | 期待値) を計算。
  - p が小さい（= 接触量のわりに不自然に結合が少ない）ペアを、BH 法で FDR 補正して候補とする。
  - 欠落ペア自身が回帰を歪めないよう、外れ値を除いて再フィット（ロバスト化）。

A–C との違い: 「接触が大きい」だけでなく「その接触量なら普通は何シナプスあるはずか」を基準にするので、
  シナプスを作りにくい細胞どうしの自然な 0 と、本当に避け合っている 0 を区別しやすい。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy import stats

from common.config import Thresholds
from common.model import Dataset, attach_weights
from common.runner import load_approach

# アプローチ C の 軸索–樹状突起 接触計算を再利用
axon_dendrite_contacts = load_approach("approach_c_axon_dendrite").axon_dendrite_contacts


def poisson_glm(X: np.ndarray, y: np.ndarray, iters: int = 50) -> np.ndarray:
    """IRLS による Poisson 回帰（statsmodels 不要）。"""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 1e-3))
    for _ in range(iters):
        eta = np.clip(X @ beta, -20, 20)
        mu = np.exp(eta)
        z = eta + (y - mu) / mu
        XtW = X.T * mu
        new = np.linalg.solve(XtW @ X + 1e-8 * np.eye(X.shape[1]), XtW @ z)
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    return beta


def nb_alpha(y: np.ndarray, mu: np.ndarray) -> float:
    """NB2 の過分散 α（モーメント法）: Var = μ + α μ²"""
    return float(max(0.0, np.sum((y - mu) ** 2 - mu) / np.sum(mu ** 2)))


def lower_tail_p(y: np.ndarray, mu: np.ndarray, alpha: float) -> np.ndarray:
    if alpha < 1e-6:
        return stats.poisson.cdf(y, mu)
    size = 1.0 / alpha
    return stats.nbinom.cdf(y, size, size / (size + mu))


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1)
    return out


def find_candidates(ds: Dataset, th: Thresholds, compartments: str = "auto", **_) -> pd.DataFrame:
    df = axon_dendrite_contacts(ds, th, compartments)
    df = attach_weights(df, ds, "pre", "post", directed=True)
    df = df[df.contact_um > 0].reset_index(drop=True)

    n_out = {i: n.n_pre for i, n in ds.neurons.items()}
    n_in = {i: n.n_post for i, n in ds.neurons.items()}
    X = np.column_stack([
        np.ones(len(df)),
        np.log(df.contact_um),
        np.log1p(df.pre.map(n_out).to_numpy(float)),
        np.log1p(df.post.map(n_in).to_numpy(float)),
    ])
    y = df.weight.to_numpy(float)

    fit = np.ones(len(df), bool)
    for _ in range(3):  # ロバスト再フィット
        beta = poisson_glm(X[fit], y[fit])
        mu = np.exp(np.clip(X @ beta, -20, 20))
        alpha = nb_alpha(y[fit], mu[fit])
        p = lower_tail_p(y, mu, alpha)
        new_fit = p >= th.robust_refit_p
        if (new_fit == fit).all():
            break
        fit = new_fit

    df["expected"] = mu
    df["p_lower"] = p
    df["q"] = bh_qvalues(p)
    df["log2_obs_over_exp"] = np.log2((y + 0.5) / (mu + 0.5))
    df["score"] = -np.log10(np.maximum(p, 1e-300))
    df["flagged"] = ((df.q < th.fdr_alpha) & (df.expected >= th.min_expected_synapses)
                     & (df.weight <= th.weak_synapse_max))
    df.attrs["model"] = dict(beta=beta.tolist(), alpha=alpha,
                             terms=["intercept", "log contact_um", "log1p n_out(pre)", "log1p n_in(post)"])
    print(f"  model: beta={np.round(beta, 3).tolist()}, NB alpha={alpha:.3f}, "
          f"fit on {fit.sum()}/{len(fit)} pairs")
    cols = ["pre", "post", "contact_um", "weight", "weight_reverse", "expected", "p_lower", "q",
            "log2_obs_over_exp", "score", "flagged", "cx", "cy", "cz"]
    return df[cols].sort_values("score", ascending=False).reset_index(drop=True)


def _args(ap):
    ap.add_argument("--compartments", default="auto", choices=["auto", "given", "navis", "synapse_vote"])


if __name__ == "__main__":
    from common.runner import cli
    cli("approach_d_peters_null_model", find_candidates, _args)
