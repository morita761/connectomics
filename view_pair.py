"""候補ペアを 3D で見る（plotly の HTML を出力。ブラウザで回転・拡大できる）。

  # ペアを直接指定
  python view_pair.py --data data/hemibrain_mb.pkl --pair 784204967 673509195
  # 結果 CSV の上位 N 件をまとめて
  python view_pair.py --data data/hemibrain_mb.pkl --csv results/hemibrain_mb/approach_a_skeleton_kdtree.csv --top 3

表示内容:
  - 2 つのニューロンのスケルトン（青系 = 1 つ目、橙系 = 2 つ目。濃い色 = 軸索、薄い色 = 樹状突起、灰 = 不明）
  - 接触箇所（赤い点）: 中心線どうしが skeleton_contact_um 以内の場所
  - シナプス（緑 = 1 つ目→2 つ目、紫 = 2 つ目→1 つ目）
凡例をクリックすると、その要素の表示/非表示を切り替えられる。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.spatial import cKDTree

from common.compartments import navis_split, synapse_vote
from common.config import Thresholds, parse_overrides
from common.geometry import cable_points
from common.model import AXON, DENDRITE, Dataset
from common.runner import ROOT

UM = 1000.0
COLORS = {0: {AXON: "#1f4e9c", DENDRITE: "#7fb2ff", 0: "#9aa5b1"},
          1: {AXON: "#b5480b", DENDRITE: "#ffb26b", 0: "#b8a99a"}}
NAMES = {AXON: "軸索", DENDRITE: "樹状突起", 0: "不明"}


def skeleton_traces(n, labels, k: int) -> list[go.Scatter3d]:
    traces = []
    child = np.flatnonzero(n.parents >= 0)
    for comp in (AXON, DENDRITE, 0):
        c = child[labels[child] == comp]
        if len(c) == 0:
            continue
        a, b = n.nodes[c] / UM, n.nodes[n.parents[c]] / UM
        xyz = np.full((len(c) * 3, 3), np.nan)
        xyz[0::3], xyz[1::3] = a, b
        name = f"{n.meta.get('type') or ''} {n.id}: {NAMES[comp]}"
        traces.append(go.Scatter3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], mode="lines", name=name,
                                   line=dict(color=COLORS[k][comp], width=2), hoverinfo="name"))
    return traces


def contact_points(na, nb, th: Thresholds) -> tuple[np.ndarray, float]:
    """接触箇所の点と、接触ケーブル長 [µm]（A と同じ定義: 両側の平均）を返す。"""
    ca = cable_points(na, th.cable_spacing_um * UM, th.default_radius_nm)
    cb = cable_points(nb, th.cable_spacing_um * UM, th.default_radius_nm)
    r = th.skeleton_contact_um * UM
    da, ka = cKDTree(cb["points"]).query(ca["points"], distance_upper_bound=r)
    db, _ = cKDTree(ca["points"]).query(cb["points"], distance_upper_bound=r)
    ha, hb = np.isfinite(da), np.isfinite(db)
    length = (ca["weight"][ha].sum() + cb["weight"][hb].sum()) / 2 / UM
    return (ca["points"][ha] + cb["points"][ka[ha]]) / 2, length


def pair_synapses(ds: Dataset, a: int, b: int) -> np.ndarray:
    if ds.synapses is None:
        return np.empty((0, 3))
    s = ds.synapses
    return s.loc[(s.pre == a) & (s.post == b), ["x", "y", "z"]].to_numpy(float)


def render(ds: Dataset, a: int, b: int, th: Thresholds, out: Path, note: str = "") -> Path:
    na, nb = ds.neurons[a], ds.neurons[b]
    traces = []
    for k, n in enumerate((na, nb)):
        lab = navis_split(n)
        if lab is None:  # 巨大なニューロンなどは、近くのシナプスの多数決で色分け
            lab = synapse_vote(n, n.nodes.astype(float), th)
        traces += skeleton_traces(n, lab, k)

    cp, contact_len = contact_points(na, nb, th)
    cp = cp / UM
    if len(cp):
        traces.append(go.Scatter3d(x=cp[:, 0], y=cp[:, 1], z=cp[:, 2], mode="markers",
                                   name=f"接触箇所（{th.skeleton_contact_um} µm 以内）",
                                   marker=dict(size=2, color="#e3242b", opacity=0.6), hoverinfo="name"))
    w = ds.directed_weights()
    for (p, q, color) in ((a, b, "#1a9e3f"), (b, a, "#7b2cbf")):
        s = pair_synapses(ds, p, q) / UM
        if len(s):
            traces.append(go.Scatter3d(x=s[:, 0], y=s[:, 1], z=s[:, 2], mode="markers",
                                       name=f"シナプス {p}→{q}（{len(s)}）",
                                       marker=dict(size=5, color=color, symbol="diamond",
                                                   line=dict(color="white", width=1))))
    title = (f"{na.meta.get('type')} {a}  ×  {nb.meta.get('type')} {b}<br>"
             f"<sup>接触 約 {contact_len:.0f} µm ／ シナプス {a}→{b}: {w.get((a, b), 0)}、"
             f"{b}→{a}: {w.get((b, a), 0)} {note}</sup>")
    fig = go.Figure(traces)
    fig.update_layout(title=title, template="plotly_white", height=850,
                      scene=dict(aspectmode="data", xaxis_title="x (µm)", yaxis_title="y (µm)",
                                 zaxis_title="z (µm)"),
                      legend=dict(itemsizing="constant"), margin=dict(l=0, r=0, t=70, b=0))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out, include_plotlyjs="cdn")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--pair", nargs=2, type=int, action="append", default=[], metavar=("A", "B"))
    ap.add_argument("--csv", default=None, help="アプローチの結果 CSV（上位 --top 件を表示）")
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()

    th = Thresholds().update(parse_overrides(args.set))
    ds = Dataset.load(args.data)
    out_dir = Path(args.out_dir or ROOT / "results" / Path(args.data).stem / "3d")

    jobs = [(a, b, "") for a, b in args.pair]
    if args.csv:
        df = pd.read_csv(args.csv)
        x, y = ("pre", "post") if "pre" in df else ("a", "b")
        src = Path(args.csv).stem
        for r, row in enumerate(df.head(args.top).itertuples()):
            jobs.append((int(getattr(row, x)), int(getattr(row, y)), f"／ {src} の {r + 1} 位"))
    if not jobs:
        ap.error("--pair か --csv を指定してください")
    for a, b, note in jobs:
        p = render(ds, a, b, th, out_dir / f"pair_{a}_{b}.html", note)
        print(p)


if __name__ == "__main__":
    main()
