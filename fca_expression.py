"""Fly Cell Atlas のクラス単位の発現を引く（FlyBase の集計ファイル）。

  python fca_expression.py kirre sns hbs rst RhoGAP100F Rdl
  python fca_expression.py --all-classes kirre sns          # 22 クラス全部

なぜこれを使うか:
  受容体や接着分子が「その細胞にあるか」を知りたいが、scRNA-seq は
  PPL105 や MBON04 のようなサブタイプを分離できない。
  一方「成虫ニューロン全般にあるか」はクラス単位の問いなので、これで解ける。
  FlyBase が FlyCellAtlas を 22 の高次クラスに集計したファイルを配っている。

値の読み方:
  「平均発現量 : 発現している細胞の割合(%)」。割合を見る。
  床（本当に無い遺伝子）は 1% 前後 — Act88F 1%、twi 0%、Hml 1%。
  一方 Mhc は 16%、Yp1 は 20% のニューロンで「検出」される。
  転写産物が大量にある遺伝子は漏れこむので、**15〜20% 以下は無いと見る。**
  逆にグリアの repo でも 14% しか出ない。この集計は感度が低く、保守的。
"""
from __future__ import annotations

import argparse
import gzip
import urllib.request
from pathlib import Path

from common.runner import ROOT

URL = ("https://s3ftp.flybase.org/releases/current/precomputed_files/genes/"
       "FlyCellAtlas_slimmed_gene_expression_fb_2026_03.tsv.gz")
CACHE = ROOT / "data" / "flycellatlas_slimmed.tsv.gz"
# 既定で見るクラス。--all-classes で 22 クラス全部
MAIN = ["neuron", "interneuron", "motor_neuron", "sensory_neuron",
        "glial_cell", "muscle_cell"]
# 較正用。床と、漏れこみの天井を示す
CALIB = ["Act88F", "twi", "Hml", "Mhc", "Yp1", "repo", "elav"]


def fetch() -> Path:
    if not CACHE.exists():
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        print(f"取得中: {URL}")
        urllib.request.urlretrieve(URL, CACHE)
    return CACHE


def table(genes: list[str]) -> tuple[list[str], dict[str, list[str]]]:
    hdr: list[str] = []
    rows: dict[str, list[str]] = {}
    want = set(genes)
    with gzip.open(fetch(), "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#gene_id"):
                hdr = line.rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 2 and p[1] in want:
                rows[p[1]] = p
    return hdr, rows


def pct(v: str) -> int | None:
    try:
        return int(v.split(":")[1].rstrip("%"))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="+")
    ap.add_argument("--all-classes", action="store_true", help="22 クラス全部を出す")
    ap.add_argument("--no-calib", action="store_true", help="較正用の遺伝子を出さない")
    args = ap.parse_args()

    genes = list(dict.fromkeys(args.genes + ([] if args.no_calib else CALIB)))
    hdr, rows = table(genes)
    cols = hdr[2:] if args.all_classes else MAIN
    idx = [hdr.index(c) for c in cols]
    short = [c.replace("_cell", "").replace("specialized_Malpighian_tubule", "Malp")[:11]
             for c in cols]

    print("発現している細胞の割合(%)。15〜20% 以下は「無い」と読む（docstring 参照）\n")
    print(f"{'遺伝子':16s}" + "".join(f"{s:>13s}" for s in short))
    print("-" * (16 + 13 * len(cols)))
    for g in genes:
        r = rows.get(g)
        if r is None:
            print(f"{g:16s}（このファイルに無い）")
            continue
        vals = [pct(r[i]) if i < len(r) else None for i in idx]
        print(f"{g:16s}" + "".join(
            f"{(str(v) + '%' if v is not None else '-'):>13s}" for v in vals))
    print("\n較正: " + " ／ ".join(CALIB) + "（--no-calib で消せる）")


if __name__ == "__main__":
    main()
