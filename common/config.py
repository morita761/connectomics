"""全アプローチ共通の閾値。

座標は内部的にすべて nm（FlyWire と同じ）。パラメータ名に _um / _um2 が付くものは µm 単位。
値は文献・慣例に基づく「目安」なので、データセットや目的に合わせて変更すること。
CLI からは各 run.py の --set key=value で上書きできる。
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields


@dataclass
class Thresholds:
    # ------------------------------------------------------------------
    # シナプス結合の強さ
    # ------------------------------------------------------------------
    # 「弱い結合」とみなす最大シナプス数（これ以下 = 結合なし or 弱い）。
    # FlyWire 全脳コネクトーム (Dorkenwald et al. 2024, Nature) では 5 シナプス以上を
    # 既定の「結合」としており、hemibrain (Scheffer et al. 2020) でも数シナプス未満の結合は
    # 自動検出の誤りを含みやすいとされる。→ 4 以下を weak とする。
    weak_synapse_max: int = 4

    # ------------------------------------------------------------------
    # 近接性
    # ------------------------------------------------------------------
    # スケルトン中心線同士の距離。中心線はふつう神経突起の太さ (半径 0.1–0.5 µm) のぶん
    # 膜表面から離れているので、膜同士の接触を拾うには約 1 µm が目安になる。
    # （哺乳類の Peters' rule 研究では spine の届く範囲として約 2 µm を使う。
    #  Stepanyants & Chklovskii 2005）
    skeleton_contact_um: float = 1.0
    # 膜表面（メッシュ頂点／半径から作るチューブ表面）同士の距離。
    # シナプス間隙は約 20 nm だが、メッシュの簡略化誤差を見込んで 0.2 µm。
    surface_contact_um: float = 0.2
    # ペアを「近接している」とみなす最小の接触ケーブル長（ノイズ的な 1 点だけの接触を除く）
    min_contact_cable_um: float = 5.0
    # 同じく最小の接触面積 (µm²)
    min_contact_area_um2: float = 2.0

    # ------------------------------------------------------------------
    # 幾何の前処理
    # ------------------------------------------------------------------
    cable_spacing_um: float = 0.5        # スケルトンを点群にするときの間隔
    default_radius_nm: float = 150.0     # 半径情報がないスケルトン用（FlyWire の L2 スケルトンなど）
    surface_points_per_ring: int = 8     # チューブ表面の 1 周あたりの点数

    # ------------------------------------------------------------------
    # 軸索／樹状突起の推定（アプローチ C, D）
    # ------------------------------------------------------------------
    compartment_k: int = 15              # 近くの何個のシナプスで多数決するか
    compartment_radius_um: float = 10.0  # 多数決に使うシナプスの探索半径

    # ------------------------------------------------------------------
    # 統計モデル（アプローチ D）
    # ------------------------------------------------------------------
    fdr_alpha: float = 0.05              # Benjamini–Hochberg の FDR
    min_expected_synapses: float = 3.0   # 「本来ならこれくらい結合しているはず」の下限
    robust_refit_p: float = 0.01         # 外れ値を除いて再フィットするときの p 値

    # ------------------------------------------------------------------
    # 細胞タイプの一貫性（アプローチ E）
    # ------------------------------------------------------------------
    # 同じタイプの組の「他のペア」が何組以上あれば、タイプ単位の期待値を信用するか
    min_type_pairs: int = 3

    # ------------------------------------------------------------------
    # 接触パッチ（アプローチ F）
    # ------------------------------------------------------------------
    # 接触点どうしがこの距離以内でつながっていれば同じパッチ（単連結クラスタリング）
    patch_link_um: float = 2.0
    # パッチ上の点からこの距離以内にあるシナプスを、そのパッチのシナプスとみなす
    patch_synapse_radius_um: float = 1.0
    # パッチ単位で評価する最小の接触ケーブル長
    min_patch_contact_um: float = 3.0
    # パッチに「本来ならこれくらいシナプスがあるはず」の下限
    min_patch_expected: float = 2.0

    # ------------------------------------------------------------------
    # 出力部位（アプローチ G）
    # ------------------------------------------------------------------
    # 自分のシナプス前部位 (T-bar) からこの距離以内のケーブルを「出力部位」とみなす。
    # hemibrain の T-bar は軸索上で数 µm おきに並ぶので 3 µm 程度でほぼ連続になる
    output_site_radius_um: float = 3.0
    # ペアの他のパッチからシナプス密度を推定するときの事前分布の強さ（擬似シナプス数）。
    # 他のパッチが小さいペアほど、全体平均の密度に引き寄せられる（経験ベイズ的な縮小推定）
    patch_prior_synapses: float = 1.0

    def update(self, overrides: dict[str, str]) -> "Thresholds":
        types = {f.name: f.type for f in fields(self)}
        for k, v in overrides.items():
            if k not in types:
                raise KeyError(f"unknown threshold: {k}  (available: {sorted(types)})")
            cur = getattr(self, k)
            setattr(self, k, type(cur)(v))
        return self

    def as_dict(self) -> dict:
        return asdict(self)


def parse_overrides(pairs: list[str] | None) -> dict[str, str]:
    out = {}
    for p in pairs or []:
        k, _, v = p.partition("=")
        out[k.strip()] = v.strip()
    return out
