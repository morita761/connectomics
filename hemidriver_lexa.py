"""split-GAL4 のヘミドライバー断片に LexA 版があるかを調べる。

  python hemidriver_lexa.py

なぜこれを見るか:
  split-GAL4 は 2 つのエンハンサー断片の積でできている（例 SS49882 = VT045661-GAL4.DBD ∩ VT039475-p65.AD）。
  その片方の断片に LexA 版があれば、**その断片だけで標的の細胞を含む集合を標識できる**。
  特異性は split より落ちる（積ではなく片側だけになる）が、
  「LexA が無いので発現系が 2 つにならない」という壁はそれで外れる。

やること:
  1. VirtualFlyBrain から、split-GAL4 の試料（JRC_SSxxxxx_…）が属する
     「P{A-GAL4.DBD} ∩ P{B-p65.AD} expression pattern」クラスを引き、断片名を取る。
  2. FlyBase の lexA allele 名（Ecol\\lexA[...]）を全部取り、断片名を抜き出す。
     GMR・VT のどちらも、GAD/p65/DBD などの表記ゆれを吸収する。
  3. 候補ペアの split-GAL4 の断片と突き合わせる。

出力: results/hemidriver_lexa.csv と標準出力。
"""
from __future__ import annotations

import argparse
import re

import pandas as pd
import requests

from common.runner import ROOT

NEO4J = "https://pdb.virtualflybrain.org/db/data/transaction/commit"
SOLR = "https://solr.virtualflybrain.org/solr/ontology/select"
# Ecol\lexA[GAD.VT019018] / Ecol\lexA[GMR58E02] / Ecol\lexA[p65.VT033006] / Ecol\lexA[VT045661]
LEXA_ALLELE = re.compile(r"lexA\[(?:[A-Za-z0-9]+\.)?((?:GMR|VT)[A-Za-z0-9]+?)(?:\.[A-Z]{2})?\]")
FRAG = re.compile(r"\{((?:GMR|VT)?[A-Za-z0-9]+)-(GAL4\.DBD|p65\.AD|GAL4\.AD|GAL4\.DBD\.?)\}")


def lexa_fragments() -> set[str]:
    """LexA 版が存在するエンハンサー断片名の集合（GMRxxxxx / VTxxxxxx）。"""
    out = set()
    for pat in ("*lexA*GMR*", "*lexA*VT*"):
        r = requests.get(SOLR, params={"q": f"label:{pat}", "rows": 5000,
                                       "wt": "json", "fl": "label"}, timeout=180)
        r.raise_for_status()
        for d in r.json()["response"]["docs"]:
            m = LEXA_ALLELE.search(d.get("label", ""))
            if m:
                out.add(m.group(1))
    return out


def hemidrivers(lines: list[str]) -> dict[str, list[tuple[str, str]]]:
    """split-GAL4 の系統名 → [(断片名, AD/DBD), ...]"""
    cy = ("MATCH (n:Individual:Split)-[:INSTANCEOF]->(c) "
          "WHERE any(l IN $lines WHERE n.label CONTAINS l) "
          "RETURN DISTINCT n.label AS sample, c.label AS combo")
    r = requests.post(NEO4J, auth=("neo4j", "vfb"), timeout=300,
                      json={"statements": [{"statement": cy, "parameters": {"lines": lines}}]})
    r.raise_for_status()
    res = r.json()
    if res.get("errors"):
        raise RuntimeError(res["errors"][:1])
    out: dict[str, list[tuple[str, str]]] = {}
    for sample, combo in (d["row"] for d in res["results"][0]["data"]):
        frags = [(a if a.startswith(("GMR", "VT")) else "GMR" + a.lstrip("R") if re.fullmatch(r"R\d{2}[A-H]\d{2}", a) else a,
                  "DBD" if "DBD" in b else "AD")
                 for a, b in FRAG.findall(combo)]
        if not frags:
            continue
        for l in lines:
            if l in sample:
                out.setdefault(l, frags)
    return out


def split_lines_of(cell: str) -> list[str]:
    """shortlist_final.csv の「SS49882〈GAL4〉(+VES019)[特異性2] / …」から系統名だけ取る。"""
    return re.findall(r"\b([A-Z]{2,3}\d{4,5})\b", cell or "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shortlist", default="results/shortlist_final.csv")
    args = ap.parse_args()

    d = pd.read_csv(ROOT / args.shortlist)
    blocked = d[(d.both_drivers == True) & (d.two_systems != True)].copy()  # noqa: E712
    good = blocked[(blocked.contact_reproduced == True) & (blocked.ok_active == True)]  # noqa: E712
    print(f"LexA が無くて落ちた組 {len(blocked)}（うち接触が再現し上流も出力している組 {len(good)}）")

    want: dict[str, set[str]] = {}
    for r in blocked.itertuples():
        for cell, t in ((r.up_gal4, r.t1), (r.down_gal4, r.t2)):
            for l in split_lines_of(str(cell)):
                want.setdefault(l, set()).add(t)
    print(f"調べる split-GAL4 系統: {len(want)}")

    lexa = lexa_fragments()
    print(f"LexA 版がある断片: {len(lexa)} 件")
    hemi = hemidrivers(sorted(want))
    print(f"ヘミドライバーが分かった系統: {len(hemi)} / {len(want)}\n")

    rows = []
    for line, frags in sorted(hemi.items()):
        hit = [(f, side) for f, side in frags if f in lexa]
        rows.append({"split_line": line,
                     "types": ";".join(sorted(want[line])),
                     "fragments": " ∩ ".join(f"{f}({s})" for f, s in frags),
                     "lexa_fragment": ";".join(f for f, _ in hit),
                     "lexa_side": ";".join(s for _, s in hit)})
    out = pd.DataFrame(rows).sort_values(["lexa_fragment", "split_line"], ascending=[False, True])
    f = ROOT / "results" / "hemidriver_lexa.csv"
    out.to_csv(f, index=False, encoding="utf-8")

    hits = out[out.lexa_fragment != ""]
    print(f"★ 片方の断片に LexA 版がある系統: {len(hits)} / {len(out)}")
    for r in hits.itertuples():
        print(f"  {r.split_line}（{r.types}）  {r.fragments}")
        print(f"      → {r.lexa_fragment}-LexA が使える（{r.lexa_side} 側の断片）")
    print(f"\n-> {f}")

    # どの組が救えるか
    print("\n【救える組】")
    for r in blocked.itertuples():
        up = [l for l in split_lines_of(str(r.up_gal4)) if l in set(hits.split_line)]
        dn = [l for l in split_lines_of(str(r.down_gal4)) if l in set(hits.split_line)]
        if not (up or dn):
            continue
        mark = "★" if (r.contact_reproduced == True and r.ok_active == True) else " "  # noqa: E712
        side = f"上流 {','.join(up)}" if up else ""
        side += ("／" if up and dn else "") + (f"下流 {','.join(dn)}" if dn else "")
        print(f" {mark} {r.t1} → {r.t2}  接触 {r.contact_src_um:.0f}µm  別個体 "
              f"{r.cross_contact_min_um:.0f}〜{r.cross_contact_max_um:.0f}µm  [{side}]")


if __name__ == "__main__":
    main()
