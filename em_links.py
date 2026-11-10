"""results/em_gap/summary.csv の測定点から neuroglancer のリンク（3D と EM 断面）を作る。

  python em_links.py --pairs MBON04:PPL105 PPL105:MBON04 MBON04:PPL107 \
      --out results/em_links_MBON04_PPL105.md

各組について summary.csv の中で直接隣接ボクセルが最大の行（= 最もよく接している測定点）を選び、
その座標に飛ぶリンクを 2 本出す。色は 赤 = 1 つ目のタイプ / 水色 = 2 つ目。
"""
from __future__ import annotations

import argparse
import json
import urllib.parse

import pandas as pd

from common.runner import ROOT

BASE = "https://neuroglancer-demo.appspot.com/#!"
EM = "precomputed://gs://neuroglancer-janelia-flyem-hemibrain/emdata/clahe_yz/jpeg"
SEG = "precomputed://gs://neuroglancer-janelia-flyem-hemibrain/v1.2/segmentation"
RED, CYAN = "#ff3b30", "#00e5ff"


def link(a: int, b: int, pos: list[int], with_em: bool) -> str:
    seg = {"type": "segmentation", "source": SEG, "name": "seg",
           "segments": [str(a), str(b)],
           "segmentColors": {str(a): RED, str(b): CYAN},
           "objectAlpha": 0.45 if with_em else 1.0}
    layers = ([{"type": "image", "source": EM, "name": "EM"}] if with_em else []) + [seg]
    state = {"dimensions": {d: [8e-09, "m"] for d in "xyz"},
             "position": [int(v) for v in pos],
             "crossSectionScale": 0.5 if with_em else 0.55,
             "projectionScale": 22000,
             "layers": layers,
             "layout": "xy-3d" if with_em else "3d",
             "showSlices": False,
             "selectedLayer": {"layer": "seg"}}
    return BASE + urllib.parse.quote(json.dumps(state, separators=(",", ":")), safe="")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", nargs="+", required=True, help="t1:t2 の形で指定")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", default=str(ROOT / "results" / "em_gap" / "summary.csv"))
    args = ap.parse_args()

    d = pd.read_csv(args.summary)
    lines = ["# 接触部位の 3D と EM（hemibrain v1.2, 8 nm 等方）", "",
             "色は 赤 = 1 つ目のタイプ / 水色 = 2 つ目。3D はメッシュ（実際の膜の形）。", ""]
    for p in args.pairs:
        t1, t2 = p.split(":")
        rows = d[(d.t1 == t1) & (d.t2 == t2)]
        if not len(rows):
            print(f"{p}: summary.csv に測定がありません")
            continue
        r = rows.loc[rows.touch_voxels.idxmax()]
        pos = [r.x, r.y, r.z]
        lines += [f"## {t1}（赤） → {t2}（水色）  直接隣接 {int(r.touch_voxels)} voxel / "
                  f"間に何も無い表面点 {int(r.direct_pct)}%",
                  f"- 3D: {link(int(r.pre), int(r.post), pos, False)}",
                  f"- 接触部位の EM: {link(int(r.pre), int(r.post), pos, True)}", ""]
    out = ROOT / args.out if not str(args.out).startswith(("/", "C:")) else args.out
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"-> {out}")


if __name__ == "__main__":
    main()
