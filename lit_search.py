"""細胞タイプごとに論文を検索し、「行動実験に使われているか」を本文から拾う（Europe PMC の公開 API）。

  python lit_search.py                     # TARGETS の全タイプ
  python lit_search.py MBON01 MBON11       # 一部だけ
  python lit_search.py --max-papers 80

やること:
  1. Europe PMC で「細胞の論文上の名前（別名）または系統名」AND Drosophila を検索
     （オープンアクセス論文は本文全体、それ以外はタイトル・要旨が検索対象）
  2. オープンアクセス論文は本文 XML を取得し、それ以外は要旨を使う
  3. 文ごとに「細胞名・系統名」「操作（光遺伝学・活動阻害など）」「行動（学習・選好など）」の語を探し、
     3 つそろう文を「行動実験の根拠」として数える
出力: results/literature/<タイプ>.csv（論文ごと）と results/literature/summary.md（根拠の文の抜粋）
取得した本文は data/literature_cache/ に保存し、2 回目以降は再利用する。

注意: 語の一致で拾っているだけなので、最終的な判断は論文本体で確認すること。
"""
from __future__ import annotations

import argparse
import json
import re
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

from common.runner import ROOT

API = "https://www.ebi.ac.uk/europepmc/webservices/rest"
CACHE = ROOT / "data" / "literature_cache"
OUT = ROOT / "results" / "literature"

# hemibrain タイプ → 論文上の名前（data/driver_lines.csv の Descriptive Cell Type から）と主な系統
TARGETS = {
    # キノコ体案で判定が上位の 4 組に出てくる細胞（2026-11-10 追加）
    "PPL105": {"names": ["PPL1-α'2α2", "PPL1-a'2a2", "PPL1-α′2α2", "MB-V1"],
               "lines": ["MB059B", "MB320B", "MB569B", "MB099C", "SS48817", "MB328B"]},
    "PPL106": {"names": ["PPL1-α3", "PPL1-a3"],
               "lines": ["SL86050", "SL86062", "MB287C", "MB328B", "MB099C", "SS01246"]},
    "MBON07": {"names": ["MBON-α1", "MBON-a1"],
               "lines": ["MB319B", "MB319C", "MB323B", "SL86055", "SS01145", "SS01180"]},
    "PAM11": {"names": ["PAM-α1", "PAM-a1"],
              "lines": ["MB325B", "SL86065", "SL86059", "SL86053", "MB043C", "MB299B"]},
    "MBON33": {"names": ["MBON-γ2γ3", "MBON-y2y3", "MBON33"],
               "lines": ["SS75199", "SS67221"]},
    # 新 6 領域の候補の下流（2026-11 追加）
    "oviIN": {"names": ["oviIN", "oviposition inhibitory neuron", "oviposition inhibitory"],
              "lines": ["SS65422", "SS65423"]},
    "FB4N": {"names": ["FB4N", "fan-shaped body FB4N"], "lines": ["SS66004", "SS95185"]},
    "FB4X": {"names": ["FB4X"], "lines": ["SS70711"]},
    "DNa01": {"names": ["DNa01"], "lines": ["SS00731"]},
    "DNa02": {"names": ["DNa02"], "lines": ["SS00730"]},
    "DNb01": {"names": ["DNb01"], "lines": ["SS02383", "SS02542"]},
    "MBON30": {"names": ["MBON-γ2α'1", "MBON30"], "lines": ["SS76094"]},
    "MBON31": {"names": ["MBON-α'1", "MBON31"], "lines": ["SS52050", "SS80896"]},
    "MBON21": {"names": ["MBON-γ'4γ'5", "MBON21"], "lines": ["SS46348", "SS81353"]},
    # 下流（MBON）
    "MBON01": {"names": ["MBON-γ5β'2a", "MBON-γ5β′2a"],
               "lines": ["SL86047", "SL86046", "SS01308", "MB011B", "MB210B", "SS01141", "SS01148"]},
    "MBON11": {"names": ["MBON-γ1pedc", "MBON-γ1pedc>α/β", "MVP2"],
               "lines": ["SL86058", "MB112C", "MB085C", "SS01359", "MB262B"]},
    "MBON04": {"names": ["MBON-β'2mp_bilateral", "β'2mp_bilateral", "MBON-β'2mp"],
               "lines": ["SL86063", "SS01125", "SS01150", "MB011B"]},
    # 上流（DAN）
    "PAM07": {"names": ["PAM-γ4<γ1γ2", "PAM-γ1γ2>γ4"], "lines": ["MB197B"]},
    "PPL103": {"names": ["PPL1-γ2α'1", "PPL1-γ2α′1"], "lines": ["MB302B", "MB296B", "MB099C"]},
    "PAM08": {"names": ["PAM-γ4"], "lines": ["MB518C", "MB312B"]},
    # 3 個体で確かめた候補（2026-11-09）
    "VL1_vPN": {"names": ["VL1_vPN", "VL1 ventral PN", "vPN VL1"], "lines": ["SS86025", "SS34947"]},
    "LHAV1a3": {"names": ["LHAV1a3", "AV1a3"], "lines": ["LH1543"]},
    "LHAV1a1": {"names": ["LHAV1a1", "AV1a1"], "lines": ["LH728", "LH1983", "LH1989"]},
    "LHAV2m1": {"names": ["LHAV2m1", "AV2m1"], "lines": []},
    "LHPD2a1": {"names": ["LHPD2a1", "PD2a1", "PD2a1/b1"], "lines": ["LH989", "LH991", "LH990"]},
    "ExR3": {"names": ["ExR3", "ER3 extrinsic", "extrinsic ring neuron 3"], "lines": ["SS82798", "SS82773"]},
    "ExR6": {"names": ["ExR6"], "lines": ["SS53617", "SS58754"]},
    "FB8F": {"names": ["FB8F", "FB8F_a"], "lines": ["SS75533", "SS69168"]},
    "PFL2": {"names": ["PFL2", "PB1-4FB1,2,4,5LAL"], "lines": ["SS82329", "SS82331"]},
    "PFL3": {"names": ["PFL3"], "lines": ["SS82335", "SS60291"]},
    "M_vPNml87": {"names": ["M_vPNml87", "multiglomerular vPN"], "lines": ["SS00587", "SS00586"]},
}

MANIP = r"cschrimson|chrimson|reachr|chr2|channelrhodopsin|gtacr|trpa1|shibire|shi\[?ts|kir2\.1|\btnt\b|tetanus|" \
        r"optogenetic|thermogenetic|silenc|block(ed|ing)|inhibit(ed|ing|ion of)|activat(ed|ing|ion of)|ablat|knock-?down|rnai"
BEHAV = r"t-maze|preference|avoid|attract|approach|learning|memory|memories|conditioning|extinction|reversal|" \
        r"sleep|feeding|courtship|odor choice|performance index|walking|arena|behavio"

S = requests.Session()


def norm(s: str) -> str:
    """名前の表記ゆれ（′/’ → '、ß → β、ガンマの代わりの y など）をそろえる。"""
    s = s.replace("′", "'").replace("’", "'").replace("ß", "β").replace("&gt;", ">")
    return s.lower()


def alias_pattern(names: list[str], lines: list[str]) -> re.Pattern:
    alts = set()
    for n in names:
        n = norm(n)
        alts.add(re.escape(n))
        alts.add(re.escape(n.replace("γ", "y")))  # 補足表などで γ の代わりに y を使う場合
    alts |= {r"\b" + re.escape(l.lower()) + r"\b" for l in lines}
    return re.compile("|".join(sorted(alts, key=len, reverse=True)))


def search(names, lines, max_papers) -> list[dict]:
    terms = [f'"{n}"' for n in names] + [f'"{l}"' for l in lines]
    query = f"({' OR '.join(terms)}) AND (Drosophila OR fly OR flies)"
    out, cursor = [], "*"
    while len(out) < max_papers:
        for wait in (0, 5, 20, 60):  # 503 などの一時的なエラーは待って再試行
            time.sleep(wait)
            r = S.get(f"{API}/search", params={"query": query, "format": "json", "resultType": "core",
                                               "pageSize": 100, "cursorMark": cursor}, timeout=60)
            if r.status_code < 500:
                break
        r.raise_for_status()
        j = r.json()
        res = j["resultList"]["result"]
        out += res
        nxt = j.get("nextCursorMark")
        if not res or nxt == cursor:
            break
        cursor = nxt
    return out[:max_papers]


def full_text(pmcid: str) -> str | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{pmcid}.txt"
    if f.exists():
        return f.read_text(encoding="utf-8") or None  # 空ファイル = 本文なし（前回取得できなかった）
    try:
        r = requests.get(f"{API}/{pmcid}/fullTextXML", timeout=45)
        if r.status_code != 200:
            f.write_text("", encoding="utf-8")
            return None
        root = ET.fromstring(r.content)
    except Exception:
        return None
    body = root.find(".//body")
    parts = [" ".join(p.itertext()) for p in (body if body is not None else root).iter("p")]
    # 図の説明文も含める（系統名は図の凡例にあることが多い）
    parts += [" ".join(c.itertext()) for c in root.iter("caption")]
    text = "\n".join(parts)
    f.write_text(text, encoding="utf-8")
    return text


def sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(])|\n", text) if len(s.strip()) > 20]


def scan(text: str, pat: re.Pattern) -> tuple[list[str], int]:
    strong, mention = [], 0
    for s in sentences(text):
        ns = norm(s)
        if not pat.search(ns):
            continue
        mention += 1
        if re.search(MANIP, ns) and re.search(BEHAV, ns):
            strong.append(s)
    return strong, mention


def run_target(t: str, spec: dict, max_papers: int) -> pd.DataFrame:
    pat = alias_pattern(spec["names"], spec["lines"])
    papers = search(spec["names"], spec["lines"], max_papers)
    # 本文の取得は 1 本 10 秒以上かかるので並列に
    oa_ids = [p["pmcid"] for p in papers if p.get("isOpenAccess") == "Y" and p.get("pmcid")]
    with ThreadPoolExecutor(8) as ex:
        texts = dict(zip(oa_ids, ex.map(full_text, oa_ids)))
    rows = []
    for p in papers:
        pmcid = p.get("pmcid")
        text = texts.get(pmcid)
        src = "full text" if text else "abstract"
        if not text:
            text = f"{p.get('title', '')}. {p.get('abstractText', '')}"
        text = re.sub(r"<[^>]+>", "", text)
        strong, mention = scan(text, pat)
        lines_found = sorted({l for l in spec["lines"] if re.search(r"\b" + re.escape(l) + r"\b", text)})
        rows.append({
            "type": t, "year": p.get("pubYear"), "title": re.sub(r"<[^>]+>", "", p.get("title", "")),
            "authors": (p.get("authorString") or "")[:80], "journal": (p.get("journalInfo") or {}).get("journal", {}).get("title"),
            "doi": p.get("doi"), "pmid": p.get("pmid"), "pmcid": pmcid, "source": src,
            "mentions": mention, "behavior_sentences": len(strong), "lines_found": ";".join(lines_found),
            "evidence": " ||| ".join(strong[:3]),
        })
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["behavior_sentences", "mentions", "year"], ascending=False)
        # 同じ論文が PubMed と PMC（やプレプリント）で重複して返ることがあるので、タイトルで 1 本にまとめる
        key = df.title.str.lower().str.replace(r"[^a-z0-9]", "", regex=True)
        df = df[~key.duplicated()]
    return df


def auto_target(t: str) -> dict:
    """TARGETS に無い細胞型を、型名と data/driver_lines.csv の系統名から組み立てる。

    LAL098 や DNb05 のように、論文でも hemibrain の型名で呼ばれる細胞のため。
    別名が分かっているものは TARGETS に書いたほうが拾いが良くなる。
    """
    names, lines = [t], []
    f = ROOT / "data" / "driver_lines.csv"
    if f.exists():
        d = pd.read_csv(f)
        hit = d[d.em_types.fillna("").str.split(";").apply(lambda xs: t in xs)]
        lines = sorted(hit.line.dropna().unique())[:12]
        for desc in hit.description.dropna():
            # 記述の先頭（「PPL1-α'2α2, ...」の 1 つ目）を論文上の名前として足す
            first = str(desc).split(",")[0].strip()
            if 2 < len(first) < 40 and first.lower() != t.lower():
                names.append(first)
    names = list(dict.fromkeys(names))[:4]
    print(f"  {t}: TARGETS に無いので自動生成 — 名前 {names} / 系統 {len(lines)} 件")
    return {"names": names, "lines": lines}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("types", nargs="*", help=f"既定: {' '.join(TARGETS)}")
    ap.add_argument("--max-papers", type=int, default=150)
    ap.add_argument("--snippets", type=int, default=6, help="summary.md に載せる論文数（タイプごと）")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# 論文検索の結果（lit_search.py）", "",
          "「細胞名・系統名」「操作」「行動」の語が同じ文にそろう論文を上から並べた。語の一致なので、本文で確認すること。", ""]
    for t in args.types or TARGETS:
        spec = TARGETS.get(t) or auto_target(t)
        df = run_target(t, spec, args.max_papers)
        df.to_csv(OUT / f"{t}.csv", index=False, encoding="utf-8")
        n_beh = int((df.behavior_sentences > 0).sum()) if len(df) else 0
        if not len(df):  # 1 本も見つからないタイプ（新しい名前など）
            md += [f"## {t}（{' / '.join(spec['names'][:2])}）", "", "論文が見つからなかった。", ""]
            continue
        print(f"{t}: 論文 {len(df)} 本、行動実験の根拠となる文がある論文 {n_beh} 本 -> {OUT / (t + '.csv')}")
        md += [f"## {t}（{' / '.join(spec['names'][:2])}）", "",
               f"検索で見つかった論文 {len(df)} 本、うち行動実験の根拠となる文がある論文 {n_beh} 本。", ""]
        for r in df[df.behavior_sentences > 0].head(args.snippets).itertuples():
            link = f"https://doi.org/{r.doi}" if r.doi else f"https://europepmc.org/article/MED/{r.pmid}"
            md.append(f"- **{r.year} {r.title}**（{r.authors}）[{r.source}] 該当文 {r.behavior_sentences}"
                      + (f"、系統 {r.lines_found.replace(';', ', ')}" if r.lines_found else "") + f" — {link}")
            ev = r.evidence.split(" ||| ")[0]
            md.append(f"  - 例: {ev[:300]}{'…' if len(ev) > 300 else ''}")
        md.append("")
    (OUT / "summary.md").write_text("\n".join(md), encoding="utf-8")
    print(f"-> {OUT / 'summary.md'}")


if __name__ == "__main__":
    main()
