"""1 つのドライバーが標識する細胞をまとめて 3D で見る。

  python view_group.py --data data/hemibrain_mb_all.pkl \
      --source PPL105 --others MBON01 MBON03 MBON04 \
      --out results/zero_3d/hemibrain_PPL105_vs_SL86063.html

なぜ必要か:
  ドライバーが複数の型を標識するとき（SL86063 は MBON01・MBON03・MBON04）、
  目標の相手以外にも人工分子が乗る。その細胞が上流の近くを通っていれば、
  そこにも意図しないシナプスができる。
  **「狙っていない細胞が、どれだけ上流の近くにいるか」を目で見るための図。**

色: 上流（--source）は赤系、ほかは型ごとに別の色。
    接触箇所（既定 1 µm 以内）は上流との組ごとに点で出す。
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

from common.config import Thresholds, parse_overrides
from common.model import Dataset
from common.runner import ROOT
from view_pair import contact_points

UM = 1000.0
# 型ごとの色（上流は最初の色）
PALETTE = ["#e3242b", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]


def cells_of(ds: Dataset, t: str) -> list[int]:
    return [i for i, n in ds.neurons.items() if n.meta.get("type") == t]


def skeleton(n, color: str, name: str, width: int = 2) -> go.Scatter3d:
    child = np.flatnonzero(n.parents >= 0)
    a, b = n.nodes[child] / UM, n.nodes[n.parents[child]] / UM
    xyz = np.full((len(child) * 3, 3), np.nan)
    xyz[0::3], xyz[1::3] = a, b
    return go.Scatter3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], mode="lines", name=name,
                        line=dict(color=color, width=width), hoverinfo="name")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--source", required=True, help="上流の型（これとの近さを測る）")
    ap.add_argument("--others", nargs="+", required=True, help="ドライバーが標識する型")
    ap.add_argument("--out", required=True)
    ap.add_argument("--contact-um", type=float, nargs="+", default=[1.0],
                    help="接触とみなす中心線どうしの距離。複数指定するとまとめて 1 枚に出す")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()

    ds = Dataset.load(args.data)
    base = Thresholds(**parse_overrides(args.set))
    src = cells_of(ds, args.source)
    if not src:
        print(f"{args.source} がデータに無い")
        return

    # 骨格は 1 度だけ描く（しきい値ごとに描き直すとファイルが倍々に大きくなる）
    traces = [skeleton(ds.neurons[i], PALETTE[0], f"{args.source} {i}（上流）", width=3)
              for i in src]
    for k, t in enumerate(args.others, start=1):
        for i in cells_of(ds, t):
            traces.append(skeleton(ds.neurons[i], PALETTE[k % len(PALETTE)], f"{t} {i}"))

    lines = []
    for ui, um_th in enumerate(sorted(args.contact_um)):
        th = replace(base, skeleton_contact_um=um_th)
        got = []
        for k, t in enumerate(args.others, start=1):
            allpts, total = [], 0.0
            for p in src:
                for i in cells_of(ds, t):
                    cp, um = contact_points(ds.neurons[p], ds.neurons[i], th)
                    total += um
                    if len(cp):
                        allpts.append(cp)
            if allpts:
                cp = np.vstack(allpts) / UM
                traces.append(go.Scatter3d(
                    x=cp[:, 0], y=cp[:, 1], z=cp[:, 2], mode="markers",
                    name=f"{um_th:g} µm 以内: {args.source} × {t}（{total:.0f} µm）",
                    marker=dict(size=3, color=PALETTE[k % len(PALETTE)], symbol="x", opacity=0.9),
                    visible=True if ui == 0 else "legendonly"))
            got.append(f"{t} {total:.0f} µm")
            print(f"  {um_th:g} µm 以内  {args.source} × {t}: 接触 {total:.1f} µm")
        lines.append(f"{um_th:g} µm 以内 — " + " ／ ".join(got))

    title = (f"{args.source} と、同じドライバーが標識する細胞の位置関係<br>"
             + "<br>".join(f"<sup>{x}</sup>" for x in lines))
    fig = go.Figure(traces)
    fig.update_layout(title=title, template="plotly_white", height=850,
                      scene=dict(aspectmode="data", xaxis_title="x (µm)",
                                 yaxis_title="y (µm)", zaxis_title="z (µm)"),
                      legend=dict(itemsizing="constant"), margin=dict(l=0, r=0, t=80, b=0))
    out = Path(args.out) if Path(args.out).is_absolute() else ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
