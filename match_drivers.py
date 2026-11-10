"""候補のタイプの組に、ドライバー（split-GAL4 / split-LexA）が両側にあるかを照合する。

  # 1) ドライバー表を読み込む（初回だけ。eLife 94168 Supplementary file 1）
  python match_drivers.py --import-elife path/to/elife-94168-supp1-v1.xlsx
  # 2) 照合
  python match_drivers.py --data data/hemibrain_mb_all.pkl

ドライバー表 data/driver_lines.csv の列:
  line, system (split-GAL4 / split-LexA), em_types（hemibrain の細胞タイプ、";" 区切り）, n_types, n_cells_lm,
  specificity (2 理想 / 1 可 / 0 行動実験には不向き), primary, description, source
drivers.csv（手で調べた分）も合わせて使う。

照合の方針:
  - 上流側は GAL4 系（split-GAL4 / GAL4）、下流側は LexA 系（split-LexA / LexA）が理想
  - 各タイプについて、標識するタイプ数が少なく、特異性が高い系統を優先
出力: results/<データ名>/driver_matches.csv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common.runner import ROOT

LINES = ROOT / "data" / "driver_lines.csv"
ELIFE = "https://elifesciences.org/articles/94168"


def import_elife(xlsx: str) -> pd.DataFrame:
    x = pd.read_excel(xlsx)
    x.columns = [re.sub(r"\s+", " ", c).strip() for c in x.columns]
    em = "EM Cell Type (best match in hemibrain v1.2.1)"
    rows = []
    for r in x.itertuples(index=False):
        d = r._asdict() if hasattr(r, "_asdict") else dict(zip(x.columns, r))
        d = dict(zip(x.columns, r))
        line = str(d["Split-GAL4 Line"]).strip()
        types = [t.strip() for t in re.split(r"[\n,;]+", str(d[em])) if t.strip() and t.strip() != "nan"]
        if not types:
            continue
        spec = d["Specificity (Ideal=2; OK=1; Not for Behavior=0)"]
        rows.append({
            "line": line,
            "system": "split-LexA" if line.upper().startswith("SL") else "split-GAL4",
            "em_types": ";".join(types), "n_types": len(types),
            "n_cells_lm": str(d["Number of Cells (LM)"]),
            "specificity": None if pd.isna(spec) else int(spec),
            "primary": d["Primary driver=1"] == 1,
            "description": str(d["Descriptive Cell Type"]).replace("\n", " "),
            "note": "" if pd.isna(d["Note"]) else str(d["Note"]),
            "source": ELIFE,
        })
    df = pd.DataFrame(rows)
    LINES.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(LINES, index=False, encoding="utf-8")
    print(f"{len(df)} lines ({(df.system == 'split-LexA').sum()} split-LexA) -> {LINES}")
    return df


def _save(rows: list[dict], name: str) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    f = LINES.parent / f"driver_lines_{name}.csv"
    df.to_csv(f, index=False, encoding="utf-8")
    print(f"{len(df)} lines, {df.em_types.str.split(';').explode().nunique()} types -> {f}")
    return df


def import_cx(xlsx: str) -> pd.DataFrame:
    """Wolff et al. eLife 104764（中心複合体）Figure 2—source data 1。細胞タイプ名は hemibrain と同じ。"""
    x = pd.read_excel(xlsx)
    x.columns = [str(c).strip() for c in x.columns]
    x["Cell type"] = x["Cell type"].ffill()
    x = x.dropna(subset=["GAL4 line number"])
    rows = []
    for r in x.itertuples(index=False):
        d = dict(zip(x.columns, r))
        t = re.sub(r"\s*\(was .*\)", "", str(d["Cell type"])).strip()
        types = [s.strip() for s in re.split(r"/", t)] if t.startswith(("PFL1", "LCNO")) else [t]
        note = "" if pd.isna(d.get("Comments")) else str(d["Comments"])
        rows.append({"line": str(d["GAL4 line number"]).strip(), "system": "split-GAL4",
                     "em_types": ";".join(types), "n_types": len(types), "n_cells_lm": "",
                     "specificity": np.nan, "primary": False, "description": str(d.get("Cell type long name", "")),
                     "note": note, "source": "https://elifesciences.org/articles/104764"})
    return _save(rows, "cx")


def import_dn(xlsx: str) -> pd.DataFrame:
    """Namiki et al. eLife 34272（下行ニューロン）Supplementary file 3。

    「DN type」は "DNa05, DNa07" のようにカンマ区切り。sample/stochastic/bg は標識の質。
    hemibrain の型名と同じ綴りなので、そのまま使える。"""
    x = pd.read_excel(xlsx, sheet_name="Sheet1")
    x.columns = [str(c).strip() for c in x.columns]
    rows = []
    for r in x.itertuples(index=False):
        d = dict(zip(x.columns, r))
        line = str(d["LineName"]).strip()
        types = [t.strip() for t in re.split(r"[,;/]+", str(d["DN type"])) if t.strip() not in ("", "nan")]
        if not line.startswith("SS") or not types:
            continue
        # bg（背景の標識）が少ないほど行動実験に向く → 特異性 2/1/0 に読み替える
        bg = d.get("bg*3")
        spec = max(0, 2 - int(bg)) if str(bg).strip().isdigit() else np.nan
        rows.append({"line": line, "system": "split-GAL4", "em_types": ";".join(types),
                     "n_types": len(types), "n_cells_lm": "", "specificity": spec, "primary": False,
                     "description": str(d.get("Genotype stock name", "")),
                     "note": f"stochastic={d.get('stochastic pattern*2')}",
                     "source": "https://elifesciences.org/articles/34272"})
    return _save(rows, "dn")


def import_cns(xlsx: str) -> pd.DataFrame:
    """Meissner et al. eLife 98405 Figure 1—source data 1。中枢神経系全体の split-GAL4 3,000 系統超。

    「Cell types」は改行・カンマ区切り。Quality は 1 が最良（1〜5）→ 特異性に読み替える。
    VNC only の系統は脳の候補には使えないので落とす。"""
    x = pd.read_excel(xlsx, sheet_name="Cell type lines")
    x.columns = [re.sub(r"\s+", " ", str(c)).strip() for c in x.columns]
    rows = []
    for r in x.itertuples(index=False):
        d = dict(zip(x.columns, r))
        if str(d.get("Adult/ Larval", "Adult")).strip() != "Adult":
            continue
        if str(d.get("VNC only", "No")).strip().lower() == "yes":
            continue
        line = str(d["Line"]).strip()
        types = [t.strip() for t in re.split(r"[\n,;]+", str(d["Cell types"])) if t.strip() not in ("", "nan")]
        if not types:
            continue
        q = d.get("Quality")
        # Quality 1 が最良。特異性は大きいほど良い向きなので反転する
        spec = np.nan if pd.isna(q) or not str(q).strip().isdigit() else max(0, 3 - int(q))
        rows.append({"line": line, "system": "split-LexA" if line.upper().startswith("SL") else "split-GAL4",
                     "em_types": ";".join(types), "n_types": len(types), "n_cells_lm": "",
                     "specificity": spec, "primary": False,
                     "description": str(d.get("Alias", "")),
                     "note": f"{d.get('First author(s) & year', '')} quality={q}",
                     "source": "https://elifesciences.org/articles/98405"})
    return _save(rows, "cns")


def import_lh(xlsx: str) -> pd.DataFrame:
    """Dolan et al. eLife 43079（外側角）Supplementary table 2。
    名前は「PV4a1」のような LH の解剖学的な名前 → hemibrain では先頭に LH が付く（LHPV4a1）。
    「PV4a1:5」は PV4a1〜PV4a5、「PV2a1/b1」は PV2a1 と PV2b1。系統名 L188 は LH188 として公開されている。"""
    x = pd.read_excel(xlsx, header=1)
    q = {"Ideal": 2, "Good": 1, "Combinatorial": 0}
    rows = []
    for r in x.itertuples(index=False):
        line, ct, qual = str(r[0]).strip(), str(r[3]).strip(), str(r[4]).strip()
        if not line.startswith("L") or ct in ("nan", ""):
            continue
        types = []
        for part in ct.split("/"):
            part = part.strip()
            m = re.fullmatch(r"([A-Z]{2}\d+[a-z])(\d+):(\d+)", part)
            if m:
                types += [f"{m[1]}{k}" for k in range(int(m[2]), int(m[3]) + 1)]
            elif re.fullmatch(r"[a-z]\d+", part) and types:  # 「PV2a1/b1」の後半
                types.append(re.match(r"[A-Z]{2}\d+", types[-1])[0] + part)
            else:
                types.append(part)
        types = [("LH" + t if re.fullmatch(r"(AV|PV|AD|PD)\d+[a-z]\d+", t) else t) for t in types]
        rows.append({"line": "LH" + line[1:], "system": "split-GAL4", "em_types": ";".join(types),
                     "n_types": len(types), "n_cells_lm": "", "specificity": q.get(qual, np.nan),
                     "primary": False, "description": ct, "note": f"{qual}; AD {r[1]} / DBD {r[2]}",
                     "source": "https://elifesciences.org/articles/43079"})
    return _save(rows, "lh")


def import_lc(docx: str) -> pd.DataFrame:
    """Wu et al. eLife 21022（視覚投射ニューロン LC・LPLC）Supplementary file 1 の 3 番目の表（全系統）。
    行動実験で使われた系統（2 番目の表）を primary にする。"""
    import xml.etree.ElementTree as ET
    import zipfile
    W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    root = ET.fromstring(zipfile.ZipFile(docx).read("word/document.xml"))
    tables = [[["".join(t.text or "" for t in tc.iter(W + "t")).strip() for tc in tr.iter(W + "tc")]
               for tr in tbl.iter(W + "tr")] for tbl in root.iter(W + "tbl")]
    # 2 番目の表の最初のブロック（対照 + 各タイプ 1 系統）が、論文の主な行動実験に使われた系統
    behav = {row[0].split()[0] for row in tables[1][2:24] if row and row[0]}
    rows, cur = [], None
    for row in tables[2][1:]:
        if len(row) < 2 or not row[1]:
            continue
        if row[0]:
            cur = row[0].split("(")[0].strip()
        types = ["LC10a", "LC10b", "LC10c", "LC10d"] if cur.startswith("LC10") else [cur]
        rows.append({"line": row[1], "system": "split-GAL4", "em_types": ";".join(types), "n_types": len(types),
                     "n_cells_lm": "", "specificity": np.nan, "primary": row[1] in behav,
                     "description": cur, "note": f"AD {row[2]} / DBD {row[3] if len(row) > 3 else ''}"
                                                 + ("; 行動実験に使用" if row[1] in behav else ""),
                     "source": "https://elifesciences.org/articles/21022"})
    return _save(rows, "lc")


def load_lines() -> pd.DataFrame:
    parts = []
    for f in sorted(LINES.parent.glob("driver_lines*.csv")):  # driver_lines.csv と driver_lines_<領域>.csv
        parts.append(pd.read_csv(f))
    manual = ROOT / "drivers.csv"
    if manual.exists():  # 手で調べた分（確認済みのみ。系統名未確認の行は除く）
        m = pd.read_csv(manual)
        m = m[(m.status == "確認済み") & ~m.driver.str.contains("未確認")]
        parts.append(pd.DataFrame({"line": m.driver, "system": m.system, "em_types": m.type, "n_types": 1,
                                   "n_cells_lm": "", "specificity": np.nan, "primary": False,
                                   "description": m.expression, "note": "drivers.csv", "source": m.source}))
    d = pd.concat(parts, ignore_index=True)
    return d.drop_duplicates(["line", "em_types"])


def is_lexa(system: str) -> bool:
    return "lexa" in str(system).lower()


def base_type(t: str) -> str:
    """hemibrain のサブタイプ接尾辞（_a, _b など）を外した名前（LHPV2a1_c → LHPV2a1）。"""
    return re.sub(r"_[a-z]$", "", t)


def best_lines(lines: pd.DataFrame, cell_type: str, want_lexa: bool | None, k: int = 3) -> list[dict]:
    # 完全一致を優先。なければサブタイプを区別しない表記（例: 表に LHPV2a1、hemibrain は LHPV2a1_c）でも拾う
    ts_all = lines.em_types.astype(str).str.split(";")
    hit = lines[ts_all.apply(lambda ts: cell_type in ts)]
    if hit.empty and base_type(cell_type) != cell_type:
        hit = lines[ts_all.apply(lambda ts: base_type(cell_type) in ts)]
        hit = hit.assign(note=hit.note.fillna("").astype(str) + "（サブタイプ未区別）")
    if want_lexa is not None:
        hit = hit[hit.system.apply(is_lexa) == want_lexa]
    if hit.empty:
        return []
    hit = hit.assign(_spec=hit.specificity.fillna(-1), _prim=hit.primary.fillna(False).astype(int))
    hit = hit.sort_values(["n_types", "_spec", "_prim"], ascending=[True, False, False])
    out = []
    for r in hit.head(k).itertuples(index=False):
        others = [t for t in str(r.em_types).split(";") if t != cell_type]
        out.append({"line": r.line, "system": r.system, "others": others,
                    "spec": None if pd.isna(r.specificity) else int(r.specificity)})
    return out


def fmt(ls: list[dict]) -> str:
    return " / ".join(f"{l['line']}" + (f"(+{','.join(l['others'])})" if l["others"] else "")
                      + (f"[特異性{l['spec']}]" if l["spec"] is not None else "") for l in ls)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--import-elife", default=None, help="eLife 94168 Supplementary file 1 (.xlsx)")
    ap.add_argument("--import-cx", default=None, help="Wolff et al. eLife 104764 Figure 2—source data 1 (.xlsx)")
    ap.add_argument("--import-lh", default=None, help="Dolan et al. eLife 43079 Supplementary file 2 (.xlsx)")
    ap.add_argument("--import-lc", default=None, help="Wu et al. eLife 21022 Supplementary file 1 (.docx)")
    ap.add_argument("--import-dn", default=None, help="Namiki 2018 の下行ニューロン表 (xlsx)")
    ap.add_argument("--import-cns", default=None, help="Meissner 2024 の中枢神経系全体の表 (xlsx)")
    ap.add_argument("--data", default=None)
    ap.add_argument("--kind", default="output", choices=["output", "directed", "undirected"])
    ap.add_argument("--all", action="store_true", help="flagged でない組も含める")
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()

    if args.import_elife:
        import_elife(args.import_elife)
    if args.import_cx:
        import_cx(args.import_cx)
    if args.import_lh:
        import_lh(args.import_lh)
    if args.import_lc:
        import_lc(args.import_lc)
    if args.import_dn:
        import_dn(args.import_dn)
    if args.import_cns:
        import_cns(args.import_cns)
    if not args.data:
        return
    lines = load_lines()
    res = ROOT / "results" / Path(args.data).stem
    tp = pd.read_csv(res / f"type_pairs_{args.kind}.csv")
    if not args.all:
        tp = tp[tp.flagged]
    rows = []
    for r in tp.itertuples(index=False):
        up_gal4 = best_lines(lines, r.t1, want_lexa=False)
        up_lexa = best_lines(lines, r.t1, want_lexa=True)
        dn_lexa = best_lines(lines, r.t2, want_lexa=True)
        dn_gal4 = best_lines(lines, r.t2, want_lexa=False)
        # 使いやすさ: 3 = 上流 GAL4 + 下流 LexA、2 = 逆向き（上流 LexA + 下流 GAL4）、
        #             1 = 両方 GAL4 系のみ、0.5 = 片側だけ、0 = なし
        if up_gal4 and dn_lexa:
            avail = 3
        elif up_lexa and dn_gal4:
            avail = 2
        elif up_gal4 and dn_gal4:
            avail = 1
        elif up_gal4 or dn_gal4 or up_lexa or dn_lexa:
            avail = 0.5
        else:
            avail = 0
        rows.append({**r._asdict(), "availability": avail,
                     "up_gal4": fmt(up_gal4), "down_lexa": fmt(dn_lexa),
                     "up_lexa": fmt(up_lexa), "down_gal4": fmt(dn_gal4),
                     "up_best_n_types": up_gal4[0]["others"].__len__() + 1 if up_gal4 else None,
                     "down_best_n_types": dn_lexa[0]["others"].__len__() + 1 if dn_lexa else None})
    out = pd.DataFrame(rows)
    out = out.sort_values(["availability", "priority" if "priority" in out else "score"], ascending=False)
    f = res / "driver_matches.csv"
    out.to_csv(f, index=False, encoding="utf-8")
    print(f"{len(out)} 組を照合 -> {f}")
    print(out.availability.value_counts().sort_index(ascending=False).rename("組数").to_string())
    cols = ["t1", "t2", "contact_total_um", "synapses", "expected", "feedforward", "memory", "availability",
            "up_gal4", "down_lexa", "down_gal4"]
    with pd.option_context("display.width", 300, "display.max_colwidth", 70):
        print(out[out.availability >= 1][cols].head(args.top).round(1).to_string(index=False))


if __name__ == "__main__":
    main()
