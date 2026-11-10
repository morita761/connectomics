"""合成データの正解ペアに対する評価指標。"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _pair_cols(df: pd.DataFrame) -> tuple[str, str, bool]:
    if "pre" in df.columns:
        return "pre", "post", True
    return "a", "b", False


def _ap(ranked_hits: np.ndarray, n_truth: int) -> float:
    if n_truth == 0:
        return float("nan")
    hits = np.cumsum(ranked_hits)
    prec = hits / np.arange(1, len(ranked_hits) + 1)
    return float((prec * ranked_hits).sum() / n_truth)


def evaluate(result: pd.DataFrame, truth: pd.DataFrame) -> dict:
    a, b, directed = _pair_cols(result)
    out = {}

    # 無向で評価（全アプローチ共通の物差し）
    key = [tuple(sorted((int(x), int(y)))) for x, y in zip(result[a], result[b])]
    r = result.assign(key=key)
    u = r.groupby("key").agg(score=("score", "max"), flagged=("flagged", "max"))
    u = u.sort_values("score", ascending=False)
    tu = {tuple(sorted((int(x), int(y)))) for x, y in zip(truth.pre, truth.post)}
    hits = np.array([k in tu for k in u.index], dtype=float)
    K = len(tu)
    out["AP_undirected"] = _ap(hits, K)
    out["P@K_undirected"] = float(hits[:K].mean()) if K else float("nan")
    fl = u.flagged.to_numpy(bool)
    out["flagged"] = int(fl.sum())
    out["precision_flagged"] = float(hits[fl].mean()) if fl.any() else 0.0
    out["recall_flagged"] = float(hits[fl].sum() / K) if K else float("nan")

    # 有向で評価（方向を持つアプローチのみ）
    if directed:
        td = {(int(x), int(y)) for x, y in zip(truth.pre, truth.post)}
        d = result.sort_values("score", ascending=False).drop_duplicates([a, b])  # F は 1 ペア複数行
        hd = np.array([(int(x), int(y)) in td for x, y in zip(d[a], d[b])], dtype=float)
        out["AP_directed"] = _ap(hd, len(td))
        out["P@K_directed"] = float(hd[:len(td)].mean())
    else:
        out["AP_directed"] = float("nan")
        out["P@K_directed"] = float("nan")
    return out


def evaluate_patches(result: pd.DataFrame, patch_truth: pd.DataFrame, tol_um: float = 3.0) -> dict:
    """パッチ単位の評価（アプローチ F 用）。同じ有向ペアで重心が tol_um 以内なら一致とみなす。"""
    if patch_truth is None or len(patch_truth) == 0 or "patch_id" not in result:
        return {}
    truth = {}
    for t in patch_truth.itertuples():
        truth.setdefault((int(t.pre), int(t.post)), []).append(np.array([t.cx, t.cy, t.cz]))
    used = set()
    r = result.sort_values("score", ascending=False)
    hits = np.zeros(len(r))
    for k, row in enumerate(r.itertuples()):
        for m, c in enumerate(truth.get((int(row.pre), int(row.post)), [])):
            key = (int(row.pre), int(row.post), m)
            if key not in used and np.linalg.norm(c - [row.cx, row.cy, row.cz]) <= tol_um * 1000:
                used.add(key)
                hits[k] = 1
                break
    K = len(patch_truth)
    fl = r.flagged.to_numpy(bool)
    return {"AP_patch": _ap(hits, K), "P@K_patch": float(hits[:K].mean()),
            "precision_flagged_patch": float(hits[fl].mean()) if fl.any() else 0.0,
            "recall_flagged_patch": float(hits[fl].sum() / K)}
