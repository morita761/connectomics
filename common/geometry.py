"""幾何処理: スケルトン → 点群 / 表面点群、近接ペアの高速検出。"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .model import AXON, DENDRITE, UNKNOWN, Neuron

CHUNK_POINTS = 50_000  # 近傍探索を何点ずつ行うか（メモリと速度のトレードオフ）


# ----------------------------------------------------------------------
# スケルトン → ケーブル点群
# ----------------------------------------------------------------------
def cable_points(n: Neuron, spacing_nm: float, default_radius_nm: float = 150.0) -> dict:
    """各エッジ（子→親）を spacing_nm 間隔の点に分割する。

    各点はそのエッジのケーブル長を等分した重み (nm) を持つので、
    「ある条件を満たす点の重みの和」=「その条件を満たすケーブル長」になる。
    """
    child = np.flatnonzero(n.parents >= 0)
    par = n.parents[child]
    a, b = n.nodes[child].astype(float), n.nodes[par].astype(float)
    seg = b - a
    L = np.linalg.norm(seg, axis=1)
    k = np.maximum(1, np.ceil(L / spacing_nm).astype(int))
    idx = np.repeat(np.arange(len(child)), k)
    starts = np.cumsum(k) - k
    j = np.arange(k.sum()) - np.repeat(starts, k)
    f = (j + 0.5) / k[idx]

    pts = a[idx] + seg[idx] * f[:, None]
    w = (L / k)[idx]

    if n.radius is not None and np.any(n.radius > 0):
        r = np.where(n.radius > 0, n.radius, default_radius_nm).astype(float)
        rad = r[child][idx] * (1 - f) + r[par][idx] * f
    else:
        rad = np.full(len(pts), default_radius_nm)

    comp = (n.compartment[child][idx] if n.compartment is not None
            else np.zeros(len(pts), dtype=np.int8))

    with np.errstate(invalid="ignore", divide="ignore"):
        t = seg / L[:, None]
    t = np.nan_to_num(t)[idx]
    return dict(points=pts, weight=w, radius=rad, compartment=comp.astype(np.int8), tangent=t)


# ----------------------------------------------------------------------
# 表面点群（メッシュがあればメッシュ、なければ半径からチューブを作る）
# ----------------------------------------------------------------------
def _orthonormal(t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = np.where(np.abs(t[:, [0]]) < 0.9, [[1.0, 0, 0]], [[0, 1.0, 0]])
    u = np.cross(t, h)
    nu = np.linalg.norm(u, axis=1, keepdims=True)
    u = np.where(nu > 0, u / np.where(nu > 0, nu, 1), [[0, 0, 1.0]])
    v = np.cross(t, u)
    return u, v


def tube_surface(cp: dict, per_ring: int) -> dict:
    """ケーブル点ごとに半径 r の円周上へ per_ring 点を置く。重みは表面積 [nm²]。"""
    u, v = _orthonormal(cp["tangent"])
    th = 2 * np.pi * np.arange(per_ring) / per_ring
    c, s = np.cos(th), np.sin(th)
    r = cp["radius"][:, None, None]
    off = r * (c[None, :, None] * u[:, None, :] + s[None, :, None] * v[:, None, :])
    pts = (cp["points"][:, None, :] + off).reshape(-1, 3)
    area = np.repeat(2 * np.pi * cp["radius"] * cp["weight"] / per_ring, per_ring)
    comp = np.repeat(cp["compartment"], per_ring)
    return dict(points=pts, weight=area, compartment=comp)


def mesh_surface(n: Neuron, cp: dict | None = None) -> dict:
    """メッシュ頂点 + 頂点面積（隣接面の面積の 1/3 ずつ）。"""
    V = n.mesh_vertices.astype(float)
    F = n.mesh_faces.astype(int)
    tri = V[F]
    fa = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    area = np.zeros(len(V))
    for col in range(3):
        np.add.at(area, F[:, col], fa / 3)
    if cp is not None and np.any(cp["compartment"] != UNKNOWN):
        _, nn = cKDTree(cp["points"]).query(V)
        comp = cp["compartment"][nn]
    else:
        comp = np.zeros(len(V), dtype=np.int8)
    return dict(points=V, weight=area, compartment=comp)


def surface_points(n: Neuron, cp: dict, per_ring: int) -> dict:
    if n.mesh_vertices is not None and n.mesh_faces is not None:
        return mesh_surface(n, cp)
    return tube_surface(cp, per_ring)


# ----------------------------------------------------------------------
# 近接ペアの検出
# ----------------------------------------------------------------------
def stack(clouds: dict[int, dict], keys=("points", "weight", "compartment")) -> dict:
    """{neuron_id: 点群} を 1 本の配列にまとめ、label（細胞 ID）列を付ける。"""
    ids = list(clouds)
    out = {k: np.concatenate([clouds[i][k] for i in ids]) for k in keys}
    out["label"] = np.concatenate([np.full(len(clouds[i]["points"]), i, dtype=np.int64) for i in ids])
    return out


def close_point_pairs(points: np.ndarray, labels: np.ndarray, radius_nm: float,
                      roles: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """距離 radius_nm 以内にある「別の細胞に属する」点のペア (i, j) を返す。

    roles=None : 無向。label[i] < label[j] になるよう並べ替える。
    roles あり : 有向。i が軸索 (AXON)、j が樹状突起 (DENDRITE) のペアだけ返す。
    """
    out_i, out_j = [], []
    for i, j in _iter_pair_chunks(points, labels, radius_nm, roles):
        out_i.append(i)
        out_j.append(j)
    if not out_i:
        return np.empty(0, int), np.empty(0, int)
    return np.concatenate(out_i), np.concatenate(out_j)


def _iter_pair_chunks(points, labels, radius_nm, roles=None, chunk=CHUNK_POINTS, post_any=False):
    """点を chunk 個ずつに分けて近傍探索し、条件を満たす点ペアを少しずつ返す。

    全点ペアを一度にメモリに載せないため（KC の軸索束のような密な領域では
    1 点あたり数百の近傍があり、一括だと数億ペアになる）。
    各点ペアはちょうど 1 回だけ返る:
      無向: label[i] < label[j] の向きのときだけ
      有向: i が軸索・j が樹状突起の向きのときだけ
      有向 (post_any=True): i が軸索（出力部位）で、j は相手のどの部位でもよい
    """
    tree = cKDTree(points)
    n = len(points)
    for s in range(0, n, chunk):
        idx = np.arange(s, min(n, s + chunk))
        rec = cKDTree(points[idx]).sparse_distance_matrix(tree, radius_nm, output_type="ndarray")
        i, j = idx[rec["i"]], rec["j"].astype(np.int64)
        if roles is None:
            keep = labels[i] < labels[j]
        else:
            keep = (labels[i] != labels[j]) & (roles[i] == AXON)
            if not post_any:
                keep &= roles[j] == DENDRITE
        if keep.any():
            yield i[keep], j[keep]


def contacts(cloud: dict, radius_nm: float, roles: np.ndarray | None = None,
             col_a: str = "a", col_b: str = "b", post_any: bool = False) -> pd.DataFrame:
    """close_point_pairs + aggregate_contacts と同じ結果を、点ペアを溜めずに少しずつ集計して返す。"""
    P, lab, w = cloud["points"], cloud["label"], cloud["weight"]
    codes, uniq = pd.factorize(lab)
    N, n = len(uniq), len(P)
    ki, kj, parts = [], [], []
    for i, j in _iter_pair_chunks(P, lab, radius_nm, roles, post_any=post_any):
        key = codes[i].astype(np.int64) * N + codes[j]
        ki.append(np.unique(key * n + i))
        kj.append(np.unique(key * n + j))
        mid = (P[i] + P[j]) / 2
        parts.append(pd.DataFrame({"key": key, "cx": mid[:, 0], "cy": mid[:, 1], "cz": mid[:, 2]})
                     .groupby("key").agg(n=("cx", "size"), sx=("cx", "sum"), sy=("cy", "sum"), sz=("cz", "sum")))
    cols = [col_a, col_b, "contact_a", "contact_b", "n_point_pairs", "cx", "cy", "cz"]
    if not parts:
        return pd.DataFrame(columns=cols)
    agg = pd.concat(parts).groupby(level=0).sum()
    ui = np.unique(np.concatenate(ki))
    uj = np.unique(np.concatenate(kj))
    ca = pd.Series(w[ui % n]).groupby(ui // n).sum()
    cb = pd.Series(w[uj % n]).groupby(uj // n).sum()
    key = agg.index.to_numpy()
    return pd.DataFrame({
        col_a: uniq[key // N], col_b: uniq[key % N],
        "contact_a": ca.reindex(key).to_numpy(), "contact_b": cb.reindex(key).to_numpy(),
        "n_point_pairs": agg.n.to_numpy(),
        "cx": (agg.sx / agg.n).to_numpy(), "cy": (agg.sy / agg.n).to_numpy(), "cz": (agg.sz / agg.n).to_numpy(),
    })[cols]


def aggregate_contacts(i: np.ndarray, j: np.ndarray, cloud: dict,
                       col_a: str = "a", col_b: str = "b") -> pd.DataFrame:
    """点ペアを細胞ペアごとに集計する。

    contact_a = 相手に近い a 側の点（重複なし）の重みの合計（ケーブル長 or 面積）
    contact_b = 同じく b 側
    n_point_pairs, 接触点の重心 (cx, cy, cz) も返す（可視化・目視確認用）。
    """
    if len(i) == 0:
        return pd.DataFrame(columns=[col_a, col_b, "contact_a", "contact_b", "n_point_pairs",
                                     "cx", "cy", "cz"])
    lab, w, P = cloud["label"], cloud["weight"], cloud["points"]
    df = pd.DataFrame({col_a: lab[i], col_b: lab[j], "pi": i, "pj": j})
    g = [col_a, col_b]
    ca = df.drop_duplicates(g + ["pi"]).assign(w=lambda d: w[d.pi]).groupby(g).w.sum()
    cb = df.drop_duplicates(g + ["pj"]).assign(w=lambda d: w[d.pj]).groupby(g).w.sum()
    mid = (P[i] + P[j]) / 2
    df[["cx", "cy", "cz"]] = mid
    agg = df.groupby(g).agg(n_point_pairs=("pi", "size"), cx=("cx", "mean"),
                            cy=("cy", "mean"), cz=("cz", "mean"))
    agg["contact_a"] = ca
    agg["contact_b"] = cb
    return agg.reset_index()[[col_a, col_b, "contact_a", "contact_b", "n_point_pairs",
                              "cx", "cy", "cz"]]


def contact_patches(i: np.ndarray, j: np.ndarray, cloud: dict, link_nm: float) -> np.ndarray:
    """接触している点ペアを「接触パッチ」に分ける（同じ細胞ペア内で、中点どうしが
    link_nm 以内につながるものを単連結クラスタリング）。点ペアごとのパッチ番号を返す。"""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    if len(i) == 0:
        return np.empty(0, dtype=int)
    lab = cloud["label"]
    key = pd.factorize(pd.MultiIndex.from_arrays([lab[i], lab[j]]))[0]
    mid = (cloud["points"][i] + cloud["points"][j]) / 2
    # 細胞ペアごとに x 方向へ大きくずらして、別ペアの点どうしがつながらないようにする
    mid = mid + np.c_[key * 1e9, np.zeros(len(key)), np.zeros(len(key))]
    pr = cKDTree(mid).query_pairs(link_nm, output_type="ndarray")
    n = len(mid)
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n)) if len(pr) else coo_matrix((n, n))
    _, comp = connected_components(g, directed=False)
    return comp
