"""全アプローチ共通のデータ構造（ライブラリ非依存）。

FlyWire / neuPrint / 合成データのどれから読み込んでも、この形に変換してから各アプローチに渡す。
"""
from __future__ import annotations

import pickle
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

UNKNOWN, AXON, DENDRITE = 0, 1, 2


@dataclass
class Neuron:
    id: int
    nodes: np.ndarray                      # (N,3) スケルトンのノード座標 [nm]
    parents: np.ndarray                    # (N,) 親ノードの index、根は -1
    radius: np.ndarray | None = None       # (N,) [nm]
    compartment: np.ndarray | None = None  # (N,) int8: 0 不明, 1 軸索, 2 樹状突起（既知の場合のみ）
    pre: np.ndarray | None = None          # (P,3) この細胞がシナプス前になる位置 [nm]
    post: np.ndarray | None = None         # (Q,3) この細胞がシナプス後になる位置 [nm]
    mesh_vertices: np.ndarray | None = None  # (M,3) [nm]
    mesh_faces: np.ndarray | None = None     # (F,3)
    meta: dict = field(default_factory=dict)  # cell type など

    @property
    def n_pre(self) -> int:
        return 0 if self.pre is None else len(self.pre)

    @property
    def n_post(self) -> int:
        return 0 if self.post is None else len(self.post)


@dataclass
class Dataset:
    neurons: dict[int, Neuron]
    edges: pd.DataFrame                    # 列: pre, post, weight（有向、シナプス数）
    name: str = "dataset"
    truth: pd.DataFrame | None = None      # 合成データのみ: 正解ペア (pre, post)
    synapses: pd.DataFrame | None = None   # 個々のシナプス: pre, post, x, y, z [nm]（アプローチ F で使用）
    meta: dict = field(default_factory=dict)

    def directed_weights(self) -> dict[tuple[int, int], int]:
        e = self.edges
        return {(int(a), int(b)): int(w) for a, b, w in zip(e.pre, e.post, e.weight)}

    def save(self, path: str | Path) -> None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(path: str | Path) -> "Dataset":
        with open(path, "rb") as f:
            return pickle.load(f)


def attach_weights(df: pd.DataFrame, ds: Dataset, a: str, b: str, directed: bool) -> pd.DataFrame:
    """候補ペア表にシナプス数の列を付ける。

    directed=True : weight = a→b, weight_reverse = b→a
    directed=False: weight = a→b + b→a（両方向の合計）
    """
    w = ds.directed_weights()
    ab = np.array([w.get((int(x), int(y)), 0) for x, y in zip(df[a], df[b])], dtype=int)
    ba = np.array([w.get((int(y), int(x)), 0) for x, y in zip(df[a], df[b])], dtype=int)
    df = df.copy()
    if directed:
        df["weight"] = ab
        df["weight_reverse"] = ba
    else:
        df["weight"] = ab + ba
        df["weight_ab"] = ab
        df["weight_ba"] = ba
    return df
