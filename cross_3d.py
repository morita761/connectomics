"""別個体（male-CNS など）での候補ペアの 3D を作る。

  python cross_3d.py --data data/male_cns_top.pkl --dataset male-cns:v1.0 \
      --pair MBON04 PPL105 --pair MBON07 PPL106 --pair MBON31 PAM11

hemibrain の 3D は results/zero_3d/ にあるが、「別個体でも同じ場所で接している」ことを
目で確かめるには、確認先の個体の 3D も要る。左右の半球それぞれで 1 枚ずつ作る。

出力: results/cross_3d/<データセット名>_<上流>_<下流>_<L|R>.html
（最終候補のページ results/shortlist.html がこの名前を探してリンクを張る）
"""
from __future__ import annotations

import argparse
from pathlib import Path

from common.config import Thresholds, parse_overrides
from common.model import Dataset
from common.runner import ROOT
from view_pair import render

OUT = ROOT / "results" / "cross_3d"


def side_of(meta: dict) -> str:
    """左右どちらの細胞か。somaSide が無ければ instance 名の末尾から拾う。"""
    s = meta.get("somaSide") or meta.get("side")
    if isinstance(s, str) and s[:1] in ("L", "R"):
        return s[0]
    inst = str(meta.get("instance") or "")
    if inst.endswith("_L"):
        return "L"
    if inst.endswith("_R"):
        return "R"
    return "?"


def cells_by_side(ds: Dataset, t: str) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    for i, n in ds.neurons.items():
        if n.meta.get("type") == t:
            out.setdefault(side_of(n.meta), []).append(i)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="確認先の個体の pickle（fetch_data.py で取得）")
    ap.add_argument("--dataset", default="male-cns:v1.0", help="ファイル名に使うデータセット名")
    ap.add_argument("--pair", nargs=2, action="append", default=[], metavar=("T1", "T2"))
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()

    ds = Dataset.load(args.data)
    th = Thresholds(**parse_overrides(args.set))
    tag = args.dataset.replace(":", "_").replace(".", "_")
    OUT.mkdir(parents=True, exist_ok=True)

    for t1, t2 in args.pair:
        a_side, b_side = cells_by_side(ds, t1), cells_by_side(ds, t2)
        if not a_side or not b_side:
            print(f"{t1} → {t2}: データに細胞が無い（{list(a_side)} / {list(b_side)}）")
            continue
        for side in ("L", "R"):
            a, b = a_side.get(side), b_side.get(side)
            if not a or not b:
                print(f"  {t1} → {t2} {side}: 片方が無い")
                continue
            # その側で最も接している組み合わせを選ぶ
            best, best_um = None, -1.0
            from view_pair import contact_points
            for x in a:
                for y in b:
                    _, um = contact_points(ds.neurons[x], ds.neurons[y], th)
                    if um > best_um:
                        best, best_um = (x, y), um
            if best is None:
                continue
            out = OUT / f"{tag}_{t1}_{t2}_{side}.html"
            render(ds, best[0], best[1], th, out,
                   note=f"{args.dataset} {side} 半球・接触 {best_um:.1f} µm")
            print(f"  {t1} → {t2} {side}: body {best[0]} → {best[1]}  接触 {best_um:.1f} µm -> {out.name}")


if __name__ == "__main__":
    main()
