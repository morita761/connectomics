"""GMR-LexA（Rubin のエンハンサー断片 LexA）を候補ペアの細胞型と突き合わせる。

  python gmr_lexa_match.py
  python gmr_lexa_match.py --refresh-lines   # VFB から LexA 系統の一覧を取り直す

なぜ必要か:
  細胞型特異な split-LexA は 20 系統・12 型しかなく、すべてキノコ体の細胞だった。
  そのため「上流・下流の両方に split-GAL4 はあるのに、LexA が無いので発現系が 2 つにならない」
  という理由だけで落ちた組が 32 組ある。GMR-LexA（約 1,400 系統）が使えれば、この壁は外れる。

やること:
  1. GMR-LexA として実在が確認できる系統名の一覧を持つ（data/gmr_lexa_lines.txt）
     FlyBase の allele 名 Ecol\\lexA[GMRxxxxx] を VirtualFlyBrain の solr から取得したもの。
  2. 候補ペアに出てくる細胞型の「論文上の名前」を作る（neuPrint の instance 名と
     data/driver_lines.csv の記述から）。
  3. data/literature_cache/ の本文を文単位で走査し、
     「GMR-LexA が実在する R 系統名」と「細胞型の名前」が同じ文に出る組を拾う。
     → その断片がその細胞を標識するという、論文で確認された根拠になる。

出力: results/gmr_lexa_match.csv（組ごと）と標準出力の一覧。

注意: GMR-LexA は 1 断片が数個〜数十の細胞型を標識する。本実験の条件
「上流が標識する型の集合 U と下流の集合 D について U × D のすべてが 0」を満たすかは、
FlyLight の画像で標識細胞を確定してから判定する必要がある。ここで出るのは候補の絞り込みまで。
"""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from common.runner import ROOT

CACHE = ROOT / "data" / "literature_cache"
LINES_FILE = ROOT / "data" / "gmr_lexa_lines.txt"
VFB_SOLR = ("https://solr.virtualflybrain.org/solr/ontology/select"
            "?q=label:*lexA*GMR*&rows=3000&wt=json&fl=label,short_form")
RLINE = re.compile(r"\bR[0-9]{2}[A-H][0-9]{2}\b")


def refresh_lines() -> list[str]:
    """VirtualFlyBrain（FlyBase 由来）から Ecol\\lexA[GMRxxxxx] の一覧を取る。"""
    import requests
    r = requests.get(VFB_SOLR, timeout=120)
    r.raise_for_status()
    out = set()
    for d in r.json()["response"]["docs"]:
        m = re.search(r"GMR(\w+)\]", d.get("label", ""))
        if m:
            out.add("R" + m.group(1))
    lines = sorted(out)
    LINES_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"GMR-LexA 系統 {len(lines)} 件 -> {LINES_FILE}")
    return lines


VFB_NEO4J = "https://pdb.virtualflybrain.org/db/data/transaction/commit"
# FlyBase の解剖学名 → hemibrain の型名（VFB が Gen1 系統に付けている注釈を読み替える）
FBBT = [
    (re.compile(r"^mushroom body output neuron (\d+)$"), lambda m: "MBON%02d" % int(m.group(1))),
    (re.compile(r"^dopaminergic PAM neuron (\d+)$"), lambda m: "PAM%02d" % int(m.group(1))),
    (re.compile(r"^dopaminergic PPL1 neuron (\d+)$"), lambda m: "PPL1%02d" % int(m.group(1))),
]


def vfb_celltypes(lexa_lines: list[str]) -> pd.DataFrame:
    """VirtualFlyBrain に登録された「この GMR 系統はどの細胞を標識するか」を取る。

    Gen1 MCFO の単一細胞画像に付いた細胞型の注釈を使う。LexA 版は GAL4 版と同じ
    エンハンサー断片なので、GAL4 版の注釈も証拠として扱う。
    """
    import requests
    labels = [f"{pre}GMR{l[1:]}-{sysname}}} expression pattern"
              for l in lexa_lines for sysname in ("GAL4", "lexA") for pre in ("P{", "PBac{")]
    cy = ("MATCH (n:Individual)-[:part_of]->(ep:Class) WHERE ep.label IN $labs "
          "MATCH (n)-[:INSTANCEOF]->(c:Class) "
          "RETURN ep.label AS ep, c.label AS cell, count(*) AS n")
    try:
        r = requests.post(VFB_NEO4J, auth=("neo4j", "vfb"), timeout=300,
                          json={"statements": [{"statement": cy, "parameters": {"labs": labels}}]})
        r.raise_for_status()
        res = r.json()
        if res.get("errors"):
            raise RuntimeError(res["errors"][:1])
        rows = [d["row"] for d in res["results"][0]["data"]]
    except Exception as e:
        print(f"  VFB に問い合わせできなかった: {e}")
        return pd.DataFrame(columns=["type", "lexa_line", "system", "fbbt"])
    out = []
    for ep, cell, _ in rows:
        t = next((f(m) for p, f in FBBT if (m := p.fullmatch(cell))), None)
        if not t:
            continue
        m = re.search(r"GMR(\w+)-(GAL4|lexA)", ep)
        out.append({"type": t, "lexa_line": "R" + m.group(1), "system": m.group(2), "fbbt": cell})
    d = pd.DataFrame(out).drop_duplicates()
    print(f"  VFB の注釈 {len(rows)} 件のうち、型に読み替えられたもの {len(d)} 件")
    return d


def norm(s: str) -> str:
    """ギリシャ文字・引用符の表記ゆれをそろえ、小文字にする。"""
    s = s.replace("′", "'").replace("’", "'").replace("ß", "β").replace("&gt;", ">")
    s = s.replace("β", "b").replace("γ", "y").replace("α", "a").replace("δ", "d").replace("′", "'")
    return s.lower()


def aliases(types: list[str]) -> dict[str, list[str]]:
    """細胞型 → 論文で使われる名前の候補。

    neuPrint の instance 名（例 PPL105(a'2a2)）から区画名を取り、
    DAN/MBON は「PPL1-a'2a2」「MBON-a1」の形に組み立てる。
    それ以外の型は型名そのものを使う（LAL098 などは論文でもその名前で呼ばれる）。
    """
    from common.loaders import neuprint_client, _retry
    from neuprint import fetch_neurons, NeuronCriteria as NC
    out: dict[str, list[str]] = {t: [norm(t)] for t in types}
    try:
        c = neuprint_client("hemibrain:v1.2.1")
        n, _ = _retry(fetch_neurons, NC(type=types, regex=False), client=c)
    except Exception as e:
        print(f"  neuPrint から instance 名を取れなかったので型名だけで照合する: {e}")
        return out
    for t, g in n.groupby("type"):
        for inst in g.instance.dropna().unique():
            m = re.match(r"^[^(]+\(([^)]+)\)", str(inst))
            if not m:
                continue
            comp = norm(m.group(1))
            if t.startswith("MBON"):
                out[t] += [f"mbon-{comp}", f"mbon {comp}"]
            elif t.startswith(("PPL", "PAM", "PAL")):
                grp = {"PPL": "ppl1", "PAM": "pam", "PAL": "pal"}[t[:3]]
                out[t] += [f"{grp}-{comp}"]
    return {t: sorted(set(v)) for t, v in out.items()}


def mine(lexa_lines: set[str], alias: dict[str, list[str]]):
    """本文を文単位で見て、R 系統名と細胞名が同じ文にある組を集める。"""
    hits: dict[tuple[str, str], set[str]] = defaultdict(set)
    sents: dict[tuple[str, str], str] = {}
    files = sorted(CACHE.glob("*.txt"))
    for f in files:
        try:
            txt = norm(f.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not txt:
            continue
        for sent in re.split(r"(?<=[.;])\s+|\n", txt):
            rs = {r for r in RLINE.findall(sent.upper()) if r in lexa_lines}
            if not rs:
                continue
            for t, als in alias.items():
                if any(a in sent for a in als if len(a) > 3):
                    for r in rs:
                        hits[(t, r)].add(f.stem)
                        sents.setdefault((t, r), sent.strip()[:300])
    return hits, sents, len(files)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh-lines", action="store_true")
    ap.add_argument("--shortlist", default="results/shortlist_final.csv")
    args = ap.parse_args()

    lexa = set(refresh_lines() if args.refresh_lines or not LINES_FILE.exists()
               else LINES_FILE.read_text(encoding="utf-8").split())
    print(f"GMR-LexA として実在が確認できる系統: {len(lexa)} 件")

    d = pd.read_csv(ROOT / args.shortlist)
    blocked = d[(d.both_drivers == True) & (d.two_systems != True)].copy()  # noqa: E712
    print(f"「両方に split-GAL4 はあるが LexA が無い」ために落ちた組: {len(blocked)}")
    good = blocked[(blocked.contact_reproduced == True) & (blocked.ok_active == True)]  # noqa: E712
    print(f"  うち別個体でも接触し、上流がその領域で出力している組: {len(good)}")

    print("\n【証拠 1】VirtualFlyBrain に登録された細胞型の注釈")
    vfb = vfb_celltypes(sorted(lexa))
    if len(vfb):
        vfb.to_csv(ROOT / "results" / "gmr_lexa_vfb.csv", index=False, encoding="utf-8")
        for t, g in vfb.groupby("type"):
            have = sorted({r for r in g.lexa_line if r in lexa})
            print(f"  {t}: {' '.join(have) if have else '（LexA 版なし）'}")
        blocked_types = set(blocked.t1) | set(blocked.t2)
        cover = sorted(blocked_types & set(vfb.type))
        print(f"  → 落ちた組の細胞型 {len(blocked_types)} のうち、VFB に注釈があるのは "
              f"{len(cover)} 件 {cover}")

    print("\n【証拠 2】論文本文での共起")
    types = sorted(set(blocked.t1) | set(blocked.t2))
    print(f"照合する細胞型: {len(types)}")
    alias = aliases(types)
    hits, sents, nfiles = mine(lexa, alias)
    print(f"文献キャッシュ {nfiles} 本を走査、該当 {len(hits)} 組\n")

    rows = []
    for (t, r), pmcs in sorted(hits.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        rows.append({"type": t, "lexa_line": r, "n_papers": len(pmcs),
                     "papers": ";".join(sorted(pmcs)), "sentence": sents[(t, r)]})
    out = pd.DataFrame(rows)
    f = ROOT / "results" / "gmr_lexa_match.csv"
    out.to_csv(f, index=False, encoding="utf-8")

    for t in sorted({r["type"] for r in rows}):
        sub = [r for r in rows if r["type"] == t]
        sides = blocked[(blocked.t1 == t) | (blocked.t2 == t)]
        role = "上流" if (blocked.t1 == t).any() else "下流"
        print(f"### {t}（{role}・関係する組 {len(sides)}）")
        for r in sub[:5]:
            print(f"    {r['lexa_line']}-LexA  論文 {r['n_papers']}  {r['papers'][:40]}")
            print(f"      > {r['sentence'][:200]}")
    print(f"\n-> {f}")


if __name__ == "__main__":
    main()
