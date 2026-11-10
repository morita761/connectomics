"""結果の一覧サイトを作る（results/index.html → データセットごとのページ → 3D 表示）。

  python build_site.py              # 全データセットのページを作る（3D は上位ペアの分だけ作成）
  python build_site.py --n3d 50     # 3D を作るペア数を増やす
  ./tasks.ps1 server               # サイトを作ってローカルサーバーで開く（serve.py と同じ）

対象: results/<名前>/ に approach_*.csv があり、data/<名前>.pkl があるもの。
ページの中身:
  - データセットの概要、アプローチ間の重なり
  - ペアの表（一致したアプローチ数・上位何 % の基準・細胞タイプで絞り込み、列クリックで並べ替え）
  - 細胞名（type / instance）と、drivers.csv に登録したドライバー
  - 3D 表示へのリンク
"""
from __future__ import annotations

import argparse
import html
import json
import re
from pathlib import Path
from urllib.parse import quote_plus

import numpy as np
import pandas as pd

from common.config import Thresholds
from common.model import Dataset
from common.runner import ROOT
from compare_methods import SHORT, load_pairs

RESULTS = ROOT / "results"
METHOD_INFO = {
    "A": ("スケルトン", "中心線どうしの距離で近接を判定（形だけを見る）"),
    "B": ("表面接触", "膜表面どうしの接触面積（形だけを見る）"),
    "C": ("軸索→樹状突起", "軸索と樹状突起の接触だけを数える（有向）"),
    "D": ("統計モデル", "接触量から期待シナプス数を回帰し、期待より少ないか"),
    "E": ("細胞タイプ", "同じタイプの組の他のペアと比べて少ないか"),
    "F": ("接触パッチ", "同じペアの他の接触部位と比べてシナプスが少ない部位"),
    "G": ("出力部位→相手", "上流の T-bar 付近（出力部位）が下流のどこかに接するのに結合しないか（有向）"),
}
MAX_PCT = 5.0  # ページに載せるのは、いずれかのアプローチで上位 5% 以内のペア


def pair_table(res: Path, ds: Dataset) -> tuple[pd.DataFrame, dict]:
    tables = {s: load_pairs(res / f"{name}.csv") for name, s in SHORT.items() if (res / f"{name}.csv").exists()}
    idx = sorted(set().union(*[set(t.index) for t in tables.values()]))
    out = pd.DataFrame(index=pd.MultiIndex.from_tuples(idx, names=["n1", "n2"]))
    for s, t in tables.items():
        out[f"pct_{s}"] = t["rank"] / len(t) * 100
        out[f"rank_{s}"] = t["rank"]
        out[f"flag_{s}"] = t["flagged"].astype("boolean")
    pcts = out[[f"pct_{s}" for s in tables]]
    out = out[pcts.min(axis=1) <= MAX_PCT].copy()
    for col, src in (("contact_um", "A"), ("contact_area_um2", "B"), ("expected", "D")):
        if src in tables and col in tables[src]:
            out[col] = tables[src][col]
    if "C" in tables:
        out["direction"] = tables["C"]["direction"]
    w = ds.directed_weights()
    n1, n2 = out.index.get_level_values(0), out.index.get_level_values(1)
    out["w12"] = [w.get((i, j), 0) for i, j in zip(n1, n2)]
    out["w21"] = [w.get((j, i), 0) for i, j in zip(n1, n2)]
    out["n_top1"] = (out[[f"pct_{s}" for s in tables]] <= 1).sum(axis=1)
    out["mean_pct"] = out[[f"pct_{s}" for s in tables]].mean(axis=1)
    out = out.sort_values(["n_top1", "mean_pct"], ascending=[False, True]).reset_index()
    sizes = {s: len(t) for s, t in tables.items()}
    return out, sizes


def neuron_info(ds: Dataset) -> dict:
    info = {}
    for i, n in ds.neurons.items():
        m = n.meta
        info[str(i)] = {k: (None if pd.isna(v) else v) for k, v in
                        {"type": m.get("type"), "instance": m.get("instance"),
                         "status": m.get("status"), "cropped": m.get("cropped"), "group": m.get("group")}.items()
                        if v is not None}
    return info


def make_3d(ds: Dataset, pairs: list[tuple[int, int]], out_dir: Path, th: Thresholds) -> None:
    from view_pair import render
    for a, b in pairs:
        f = out_dir / f"pair_{a}_{b}.html"
        g = out_dir / f"pair_{b}_{a}.html"
        if not f.exists() and not g.exists():
            print(f"  3D: {a} - {b}")
            render(ds, a, b, th, f)


def driver_rows(types: set[str]) -> list[dict]:
    parts = []
    f = ROOT / "drivers.csv"  # 手で調べた分
    if f.exists():
        parts.append(pd.read_csv(f))
    # match_drivers.py --import-elife / --import-cx / --import-lh / --import-lc で取り込んだ表
    for g in sorted((ROOT / "data").glob("driver_lines*.csv")):
        d = pd.read_csv(g)
        d["type"] = d.em_types.astype(str).str.split(";")
        d = d.explode("type")
        others = d.apply(lambda r: [t for t in str(r.em_types).split(";") if t != r.type], axis=1)
        d["expression"] = [
            (f"{desc}" + (f"（ほかに {', '.join(o)} も標識）" if o else "（このタイプのみ）")
             + (f"、LM {n} 個" if str(n) not in ("nan", "", "-") else "")
             + (f"、特異性 {int(sp)}" if sp == sp else ""))
            for desc, o, n, sp in zip(d.description, others, d.n_cells_lm, d.specificity)]
        d["driver"], d["status"] = d.line, "確認済み"
        parts.append(d[["type", "driver", "system", "expression", "status", "source"]])
    if not parts:
        return []
    d = pd.concat(parts, ignore_index=True).drop_duplicates(["type", "driver"])
    return d[d.type.isin(types)].fillna("").to_dict("records")


def lookup_links(t: str) -> list[tuple[str, str]]:
    q = quote_plus
    return [("Google（split-GAL4）", f"https://www.google.com/search?q={q(t + ' split-GAL4 Drosophila')}"),
            ("Google（LexA）", f"https://www.google.com/search?q={q(t + ' LexA Drosophila')}"),
            ("Virtual Fly Brain", f"https://www.google.com/search?q={q(t + ' site:virtualflybrain.org')}"),
            ("FlyLight split-GAL4", "https://splitgal4.janelia.org/cgi-bin/splitgal4.cgi")]


def build_dataset_page(stem: str, n3d: int, th: Thresholds, n_type3d: int = 15) -> dict | None:
    res = RESULTS / stem
    data = ROOT / "data" / f"{stem}.pkl"
    if not data.exists() or not any(res.glob("approach_*.csv")):
        return None
    print(f"[{stem}]")
    ds = Dataset.load(data)
    tab, sizes = pair_table(res, ds)
    methods = list(sizes)

    out3d = res / "3d"
    want = list(zip(tab.n1.head(n3d), tab.n2.head(n3d)))
    for s in methods:  # 各アプローチの上位 5 件も
        top = tab.sort_values(f"rank_{s}").head(5)
        want += list(zip(top.n1, top.n2))
    type_pairs = {}
    for key in ("output", "directed", "undirected"):
        f = res / f"type_pairs_{key}.csv"
        if f.exists():
            tp = pd.read_csv(f)
            tp = tp.sort_values(["flagged", "priority" if "priority" in tp else "score"], ascending=False)
            type_pairs[key] = tp
            ex = tp.dropna(subset=["ex_a"]).head(n_type3d)
            want += [(int(x), int(y)) for x, y in zip(ex.ex_a, ex.ex_b)]
    make_3d(ds, list(dict.fromkeys(want)), out3d, th)
    have = {p.stem for p in out3d.glob("pair_*.html")}

    def link_for(a, b):
        return (f"3d/pair_{a}_{b}.html" if f"pair_{a}_{b}" in have else
                f"3d/pair_{b}_{a}.html" if f"pair_{b}_{a}" in have else None)

    tp_rows = {}
    for key, tp in type_pairs.items():
        recs = []
        for r in tp.head(400).itertuples(index=False):
            d = {k: (None if pd.isna(v) else (round(float(v), 3) if isinstance(v, float) else v))
                 for k, v in r._asdict().items()}
            if d.get("ex_a") is not None:
                d["ex_a"], d["ex_b"] = str(int(d["ex_a"])), str(int(d["ex_b"]))
                d["link"] = link_for(int(d["ex_a"]), int(d["ex_b"]))
            d["flagged"] = bool(r.flagged)
            recs.append(d)
        tp_rows[key] = recs

    rows = []
    for r in tab.itertuples(index=False):
        a, b = int(r.n1), int(r.n2)
        link = (f"3d/pair_{a}_{b}.html" if f"pair_{a}_{b}" in have else
                f"3d/pair_{b}_{a}.html" if f"pair_{b}_{a}" in have else None)
        row = {"a": str(a), "b": str(b), "w12": int(r.w12), "w21": int(r.w21), "link": link,
               "pct": {s: (None if pd.isna(getattr(r, f"pct_{s}")) else round(float(getattr(r, f"pct_{s}")), 3))
                       for s in methods},
               "rank": {s: (None if pd.isna(getattr(r, f"rank_{s}")) else int(getattr(r, f"rank_{s}")))
                        for s in methods},
               "flag": {s: bool(getattr(r, f"flag_{s}")) if not pd.isna(getattr(r, f"flag_{s}")) else False
                        for s in methods}}
        for col in ("contact_um", "contact_area_um2", "expected"):
            if col in tab and not pd.isna(getattr(r, col)):
                row[col] = round(float(getattr(r, col)), 2)
        if "direction" in tab and isinstance(r.direction, str):
            row["direction"] = r.direction
        rows.append(row)

    neurons = neuron_info(ds)
    types = {v.get("type") for v in neurons.values() if v.get("type")}
    drivers = driver_rows(types)
    payload = {
        "name": region_label(stem)[0], "label_desc": region_label(stem)[1], "stem": stem,
        "dataset": ds.name, "methods": methods, "sizes": sizes,
        "method_info": {s: METHOD_INFO[s] for s in methods},
        "neurons": neurons, "rows": rows, "drivers": drivers,
        "lookup": {t: lookup_links(t) for t in sorted(types)},
        "n_neurons": len(ds.neurons), "n_edges": int(len(ds.edges)),
        "has_instance": any("instance" in v for v in neurons.values()),
        "type_pairs": tp_rows,
    }
    page = DATASET_TEMPLATE.replace("__TITLE__", html.escape(stem)).replace(
        "__DATA__", json.dumps(payload, ensure_ascii=False).replace("</", "<\\/"))
    (res / "index.html").write_text(page, encoding="utf-8")
    best = tab.head(1)
    return {"name": stem, "dataset": ds.name, "n_neurons": len(ds.neurons), "methods": "".join(methods),
            "n_rows": len(rows), "n_all": int((tab.n_top1 == len(methods)).sum()),
            "types": ", ".join(sorted(types))[:120]}


REGION_LABELS = {
    "hemibrain_mb_all": ("キノコ体まわり 全タイプ", "記憶の中枢。KC・MBON・DAN・投射ニューロンなど"),
    "hemibrain_mb": ("キノコ体 KCg-m × MBON01", "最初の試験用。2 タイプだけ"),
    "hemibrain_lh": ("外側角（LH）", "匂いの生得的な価値づけ。投射ニューロン → 外側角ニューロン"),
    "hemibrain_cx": ("中心複合体（EB・FB・PB・NO）", "方向と航行。リング → 柱状 → 扇状体 → 出力"),
    "hemibrain_vis": ("視覚の投射先（AOTU・PVLP）", "視覚投射ニューロン（LC・LT）とその下流"),
    "hemibrain_al": ("触角葉（AL）", "嗅覚の入口。ORN → 投射ニューロン → 局所介在ニューロン"),
    "hemibrain_smp": ("上前大脳（SMP・SIP・SLP）", "MBON の出力先。下行ニューロンや神経ペプチド細胞"),
    "hemibrain_lal": ("側方付属葉（LAL・VES・GOR）", "前運動。中心複合体 → 下行ニューロン"),
    "hemibrain_vis_olf": ("視覚系 × 嗅覚系", "系をまたぐ接触の調査"),
}


def region_label(stem: str) -> tuple[str, str]:
    if stem in REGION_LABELS:
        return REGION_LABELS[stem]
    if stem.startswith("synthetic"):
        return ("合成データ（正解つき）", "手法の比較・検証用。実際の脳ではない")
    return (stem, "")


def find_3d(t1: str, t2: str) -> list[tuple[str, str]]:
    """最終候補ページから張る 3D のリンク（results/zero_3d と cross_3d を探す）。

    3D のファイル名は作ったときの向きで決まるので、t1_t2 と t2_t1 の両方を探す。
    同じ 2 細胞なら向きが逆でも中身は同じ（接触点と両方のスケルトン）。
    """
    seen, out = set(), []
    for d in (RESULTS / "zero_3d", RESULTS / "cross_3d"):
        if not d.exists():
            continue
        for pat in (f"*{t1}_{t2}*.html", f"*{t2}_{t1}*.html"):
            for f in sorted(d.glob(pat)):
                who = "hemibrain" if f.stem.startswith("hemibrain") else (
                    "male-cns" if "male" in f.stem else "別個体")
                if who in seen:  # 同じ個体の 3D は 1 つだけ
                    continue
                seen.add(who)
                out.append((f"{d.name}/{f.name}", f"{who} の 3D"))
    return out


def find_docs(t1: str, t2: str) -> list[tuple[str, str]]:
    """その組について書いた報告書・EM リンク集（results/docs/ の HTML）。

    build_docs.py が results/*.md を results/docs/*.html にしてある前提。
    publish_site.py で公開するときは docs/ -> reports/ に読み替える。
    """
    out = []
    for a, b in ((t1, t2), (t2, t1)):
        for stem, label in ((f"candidate_{a}_{b}", "検討"), (f"em_links_{a}_{b}", "EM")):
            f = RESULTS / "docs" / f"{stem}.html"
            if f.exists() and (f"docs/{stem}.html", label) not in out:
                out.append((f"docs/{stem}.html", label))
    return out


def em_gap_of() -> dict[tuple[str, str], tuple[int, int]]:
    """EM で測った膜の密着（results/em_gap/summary.csv）を組ごとに 1 つにまとめる。

    1 組につき測定点が複数あるので、最もよく接している点を代表にする。
    """
    f = RESULTS / "em_gap" / "summary.csv"
    if not f.exists():
        return {}
    d = pd.read_csv(f)
    out: dict[tuple[str, str], tuple[int, int]] = {}
    for r in d.itertuples(index=False):
        k = (str(r.t1), str(r.t2))
        v = (int(r.touch_voxels), int(r.direct_pct))
        if k not in out or v[0] > out[k][0]:
            out[k] = v
    return out


GRADES = [
    ("◎", "上流と下流に<b>別々の発現系</b>（GAL4 と LexA）がそろい、両者がその領域で働いている"
          " … 上流と下流に別の人工分子を入れられる"),
    ("○", "両側にドライバーはあるが<b>どちらも GAL4</b>。2 つ目の発現系（LexA・QF など）を作る必要がある"),
    ("△", "片側しかドライバーが無い、または無い"),
    ("要注意", "上流がその領域でほとんど出力していない（放出できない場所での接触）"),
    ("接触未確認", "別個体での接触をまだ計算していない"),
]


def _i(v) -> int:
    return 0 if v != v else int(v)  # NaN は 0 扱い


def grade_of(r) -> tuple[str, str]:
    if r.cross_contact_min_um != r.cross_contact_min_um:  # NaN
        return "接触未確認", "別個体での接触をまだ計算していない"
    if r.pre_out_in_main != r.pre_out_in_main:
        return "接触未確認", "接触領域をまだ調べていない（inspect_pairs.py 未実行）"
    if not r.ok_active:
        return "要注意", f"{r.main_roi} で上流の出力が {_i(r.pre_out_in_main)} 個しかない"
    note = ""
    if _i(r.pre_out_in_main) < 300:
        note = f"（上流の出力は {_i(r.pre_out_in_main)} 個とやや少ない）"
    if getattr(r, "two_systems", False):
        return "◎", "上流 GAL4 + 下流 LexA がそろう" + note
    if r.both_drivers:
        return "○", "両側にドライバーあり。ただし<b>どちらも GAL4</b>" + note
    if r.up_driver != "-" or r.down_driver != "-":
        return "△", "片側だけドライバーあり" + note
    return "△", "ドライバーが見つかっていない" + note


def build_shortlist_page() -> dict | None:
    """3 個体で確かめた最終候補のページ（results/shortlist.html）。"""
    f = RESULTS / "shortlist_final.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d["ok_active"] = d.ok_active.fillna(False).astype(bool)
    d["both_drivers"] = d.both_drivers.fillna(False).astype(bool)
    rank = {g: i for i, (g, _) in enumerate(GRADES)}
    EM = em_gap_of()
    rows = []
    for r in d.itertuples(index=False):
        g, why = grade_of(r)
        links = " ".join(f'<a href="{html.escape(u)}">{html.escape(lab)}</a>' for u, lab in find_3d(r.t1, r.t2))
        docs = " ".join(f'<a class="rep" href="{html.escape(u)}">{html.escape(lab)}</a>'
                        for u, lab in find_docs(r.t1, r.t2))
        g_em = EM.get((r.t1, r.t2)) or EM.get((r.t2, r.t1))
        em_txt = f"{g_em[0]:,} / {g_em[1]}%" if g_em else ""
        em_n = g_em[0] if g_em else -1
        def cell(gal4, lexa):
            parts = []
            for val, label, cls in ((gal4, "GAL4", ""), (lexa, "LexA", ' class="lexa"')):
                if isinstance(val, str) and val != "-":
                    txt = val.replace("〈GAL4〉", "").replace("〈LexA〉", "")  # 行頭のラベルと重複するので消す
                    parts.append(f'<div{cls}>{label}: {html.escape(txt)}</div>')
            return "".join(parts) or "—"
        drv_up = cell(getattr(r, "up_gal4", "-"), getattr(r, "up_lexa", "-"))
        drv_dn = cell(getattr(r, "down_gal4", "-"), getattr(r, "down_lexa", "-"))
        # 「flywire:783(合計 0/0), male-cns:v1.0(合計 0/0)」→「FlyWire 0/0・male-cns 0/0」
        ev = re.sub(r"(\w[\w-]*):[\w.]+\(合計\s*", r"\1 ", str(r.datasets or "")).replace(")", "")
        ev = ev.replace("flywire", "FlyWire").replace(", ", "・")
        ev = f"hemibrain {int(r.syn_src_fwd)}/{int(r.syn_src_rev)}・{ev}"
        cross = re.sub(r"male-cns\s*", "", str(r.cross_contact_detail or "")).replace("µm", "").strip()
        rows.append({
            "_g": rank.get(g, 9), "grade": g, "why": why,
            "pair": f'<span class="cell">{html.escape(r.t1)}</span> → <span class="cell">{html.escape(r.t2)}</span>',
            "t1": r.t1, "t2": r.t2, "region": r.main_roi or "", "rois": r.contact_rois or "",
            "contact": round(float(r.contact_src_um), 1),
            "contact_txt": (f"{float(r.contact_src_um):.0f} → "
                            f"{float(r.cross_contact_min_um):.0f}–{float(r.cross_contact_max_um):.0f} µm"),
            "cross_contact": f"別個体 {cross}",
            "cross_min": round(float(r.cross_contact_min_um), 1),
            "syn": int(r.total_syn), "datasets": ev,
            "pre_out": _i(r.pre_out_in_main), "post_in": _i(r.post_in_in_main),
            "active": (f"{_i(r.pre_out_in_main)} / {_i(r.post_in_in_main)}"
                       if r.pre_out_in_main == r.pre_out_in_main else "未確認"),
            "comp": f"{r.comp_type} {float(r.comp_rate):.2f}" if isinstance(r.comp_type, str) else "",
            "up": drv_up, "down": drv_dn, "two": bool(getattr(r, "two_systems", False)),
            "sysnote": ("上流 GAL4 + 下流 LexA" if getattr(r, "two_systems", False)
                        else ("両方 GAL4 のみ" if r.both_drivers else "ドライバー不足")),
            "readout": r.readout if isinstance(r.readout, str) and r.readout != "-" else "",
            "links": links + ("  " if links and docs else "") + docs,
            "em": em_txt, "em_n": em_n,
        })
    rows.sort(key=lambda x: (x["_g"], -x["contact"]))
    legend = "".join(f'<li><b>{html.escape(g)}</b> … {t}</li>' for g, t in GRADES)
    page = (SHORTLIST_TEMPLATE
            .replace("__ROWS__", json.dumps(rows, ensure_ascii=False))
            .replace("__LEGEND__", legend)
            .replace("__N__", str(len(rows)))
            .replace("__NBEST__", str(sum(1 for x in rows if x["grade"] == "◎"))))
    (RESULTS / "shortlist.html").write_text(page, encoding="utf-8")
    print(f"-> {RESULTS / 'shortlist.html'}")
    return {"n": len(rows), "n_best": sum(1 for x in rows if x["grade"] == "◎")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n3d", type=int, default=30, help="一致数の多い順に、3D を作るペア数")
    args = ap.parse_args()
    th = Thresholds()
    short = build_shortlist_page()
    entries = []
    for d in sorted(p for p in RESULTS.iterdir() if p.is_dir()):
        e = build_dataset_page(d.name, args.n3d, th)
        if e:
            entries.append(e)
    top = ""
    if short:
        top = (f'<a class="card hero" href="shortlist.html"><h2>★ 最終候補 — 3 個体で確かめた組（{short["n"]} 組）</h2>'
               f'<p>hemibrain・male-cns・FlyWire の <b>3 個体すべてで両方向のシナプスが 0</b> で、'
               f'別個体でも接触が再現した組。<b>{short["n_best"]} 組</b>は両側にドライバーがある。</p>'
               f'<p class="muted">実験に進めるならここから選ぶ。下のデータセット別ページは探索の元データ。</p></a>')
    docs_index = RESULTS / "docs" / "index.html"
    if docs_index.exists():
        n_docs = len(list((RESULTS / "docs").glob("*.html"))) - 1
        top += (f'<a class="card hero" href="docs/index.html"><h2>◇ 報告書（{n_docs} 件）</h2>'
                f'<p>候補の検討・ドライバー・文献のまとめ。'
                f'results/*.md を HTML にしたもの（build_docs.py）。</p></a>')

    items = []
    for e in entries:
        label, desc = region_label(e["name"])
        items.append(
            f'<a class="card" href="{html.escape(e["name"])}/index.html"><h2>{html.escape(label)}</h2>'
            f'<p class="muted">{html.escape(desc)}</p>'
            f'<p>{e["n_neurons"]} ニューロン・アプローチ {e["methods"]}'
            f'<span class="sub"> （データ: {html.escape(e["name"])} ／ {html.escape(e["dataset"])}）</span></p>'
            f'<p><b>{e["n_all"]}</b> 組が全アプローチの上位 1% に一致 ／ 一覧 {e["n_rows"]} 組</p></a>')
    body = top + "\n".join(items) if items else "<p>結果がありません。run_all.py を実行してください。</p>"
    (RESULTS / "index.html").write_text(INDEX_TEMPLATE.replace("__ITEMS__", body), encoding="utf-8")
    print(f"-> {RESULTS / 'index.html'}")


BASE_CSS = """
:root{--bg:#f7f8fa;--fg:#1d2330;--muted:#5b6475;--card:#fff;--line:#e2e5ea;--accent:#2457c5;--good:#1a7f45;--warn:#a15c00;--chip:#eef2fb}
@media (prefers-color-scheme:dark){:root{--bg:#14171c;--fg:#e6e9ef;--muted:#9aa3b2;--card:#1c2027;--line:#2d333d;--accent:#7aa2ff;--good:#5cc98a;--warn:#e0a650;--chip:#232a37}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.55 system-ui,-apple-system,"Segoe UI","Hiragino Sans","Yu Gothic UI",sans-serif}
main{max-width:1280px;margin:0 auto;padding:24px 16px 64px}h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:28px 0 10px}
a{color:var(--accent)}.muted{color:var(--muted)}.card{display:block;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;margin:12px 0;text-decoration:none;color:inherit}
.card:hover{border-color:var(--accent)}.card h2{margin:0 0 6px}.card p{margin:4px 0}
"""

INDEX_TEMPLATE = """<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>コネクトミクス結果一覧</title><style>""" + BASE_CSS + """
.card.hero{border-color:var(--accent);border-width:2px}.card.hero h2{color:var(--accent)}
</style></head><body><main>
<h1>「近いのに結合しない」ペア探索 — 結果一覧</h1>
<p class="muted">上の「最終候補」が結論です。下は領域ごとの探索結果（元データ）で、絞り込み前の候補が入っています。</p>
__ITEMS__
</main></body></html>"""

SHORTLIST_TEMPLATE = """<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>最終候補 — 3 個体で確かめた組</title><style>""" + BASE_CSS + """
.tablewrap{overflow-x:auto;background:var(--card);border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:7px 9px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none;white-space:nowrap}
th.sorted::after{content:" ▼"}th.sorted.asc::after{content:" ▲"}
td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.cell{font-weight:600}.sub{color:var(--muted);font-size:12px}
.g{display:inline-block;min-width:3.6em;text-align:center;padding:1px 8px;border-radius:999px;font-weight:700;white-space:nowrap}
td:first-child{min-width:11em;max-width:15em}td:nth-child(2){min-width:11em}
td:nth-child(3){max-width:12em}td:nth-child(8),td:nth-child(9){max-width:19em;white-space:normal}
td .sub{white-space:normal}td:nth-child(4),td:nth-child(5){min-width:9em}
.g.best{background:var(--good);color:#fff}.g.ok{background:var(--accent);color:#fff}
.g.maybe{background:var(--chip);color:var(--fg)}.g.warn{background:var(--warn);color:#fff}
.controls{display:flex;flex-wrap:wrap;gap:12px 20px;align-items:center;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px;margin:10px 0}
.controls label{display:flex;gap:6px;align-items:center}
select,input{font:inherit;padding:4px 6px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg)}
.note{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 16px}
.note ul{margin:6px 0 0 18px;padding:0}.note li{margin:3px 0}
.drv{color:var(--good)}.drv.no{color:var(--muted)}.sub.okdrv{color:var(--good);font-weight:700}
.lexa{color:var(--accent);font-weight:600}
.g.none{background:var(--chip);color:var(--muted)}
</style></head><body><main>
<p><a href="index.html">← 結果一覧</a></p>
<h1>最終候補 — 3 個体で確かめた組（__N__ 組）</h1>
<p class="muted">hemibrain（雌・右半球）・male-cns（雄・全神経系）・FlyWire（雌・全脳）の
<b>3 個体すべてで両方向のシナプスが 0</b> だった組。人工分子でシナプスを作らせる実験を想定した並び。<br>
別個体での接触と接触領域まで確認した組に判定を付けている（「接触未確認」はこれから調べる分）。
<b>上流と下流に別の分子を入れるには 2 つの発現系（GAL4 と LexA）が要る</b>ので、そこを判定の軸にしている。</p>
<div class="note">
<b>判定の意味</b><ul>__LEGEND__</ul>
<p class="sub" style="margin-top:8px">「上流出力 / 下流入力」は主な接触領域にある上流のシナプス前部位と下流のシナプス後部位の個数。
どちらも多いほど「通り道が重なっているだけ」ではない。<br>
「比較」は同じ下流タイプに結合している別の上流タイプとその密度（シナプス/µm）。ここが大きいほど「結合できるのにしていない」と言える。</p>
</div>
<div class="controls">
  <label>判定 <select id="fg"><option value="">すべて</option><option value="◎">◎ のみ</option><option value="◎○">◎ と ○</option></select></label>
  <label><input type="checkbox" id="fs"> <b>2 つの発現系がそろう組だけ</b>（上流 GAL4 + 下流 LexA）</label>
  <label><input type="checkbox" id="fd"> ドライバーが両側にある組だけ</label>
  <label><input type="checkbox" id="fr"> 行動の読み出しがある組だけ</label>
  <label>細胞名で絞り込み <input id="q" placeholder="例: LHPV6c1"></label>
  <span class="sub" id="cnt"></span>
</div>
<div class="tablewrap"><table id="t"><thead><tr>
<th data-k="_g">判定</th><th data-k="t1">上流 → 下流</th><th data-k="region">主な接触領域</th>
<th data-k="contact" class="num">接触 元 → 別個体</th>
<th data-k="syn" class="num">シナプス<br><span class="sub">3 個体合計</span></th>
<th data-k="pre_out" class="num">上流出力 / 下流入力</th>
<th data-k="em_n" class="num">膜の密着<br><span class="sub">EM・直接隣接 / 割合</span></th>
<th data-k="comp">比較（同じ下流へ）</th><th data-k="up">上流ドライバー</th><th data-k="down">下流ドライバー</th><th>3D・報告書</th>
</tr></thead><tbody></tbody></table></div>
<script>
const ROWS = __ROWS__;
const GC = {"◎":"best","○":"ok","△":"maybe","要注意":"warn","接触未確認":"none"};
let sortK = "_g", asc = true;
function view(){
  const g = document.getElementById('fg').value, dOnly = document.getElementById('fd').checked;
  const rOnly = document.getElementById('fr').checked, q = document.getElementById('q').value.trim().toLowerCase();
  const sOnly = document.getElementById('fs').checked;
  let rs = ROWS.filter(r => (!g || g.includes(r.grade))
    && (!sOnly || r.two)
    && (!dOnly || (r.up !== '—' && r.down !== '—'))
    && (!rOnly || r.readout)
    && (!q || (r.t1 + ' ' + r.t2).toLowerCase().includes(q)));
  rs.sort((a, b) => {
    const x = a[sortK], y = b[sortK];
    const c = (typeof x === 'number' && typeof y === 'number') ? x - y : String(x).localeCompare(String(y), 'ja');
    return asc ? c : -c;
  });
  document.querySelector('#t tbody').innerHTML = rs.map(r => `<tr>
    <td><span class="g ${GC[r.grade] || ''}">${r.grade}</span><div class="sub">${r.why}</div></td>
    <td>${r.pair}${r.readout ? `<div class="sub">読み出し: ${r.readout}</div>` : ''}</td>
    <td><div class="sub">${r.rois}</div></td>
    <td class="num">${r.contact_txt}<div class="sub">${r.cross_contact}</div></td>
    <td class="num">${r.syn}<div class="sub">${r.datasets}</div></td>
    <td class="num">${r.active}</td>
    <td class="num">${r.em || '<span class="sub">未測定</span>'}</td>
    <td>${r.comp}</td>
    <td class="${r.up === '—' ? 'drv no' : 'drv'}">${r.up}</td>
    <td class="${r.down === '—' ? 'drv no' : 'drv'}">${r.down}
      <div class="sub ${r.two ? 'okdrv' : ''}">${r.two ? '✓ ' : ''}${r.sysnote}</div></td>
    <td>${r.links || '<span class="sub">—</span>'}</td></tr>`).join('');
  document.getElementById('cnt').textContent = `${rs.length} / ${ROWS.length} 組`;
  document.querySelectorAll('#t th').forEach(th => {
    th.classList.toggle('sorted', th.dataset.k === sortK);
    th.classList.toggle('asc', th.dataset.k === sortK && asc);
  });
}
document.querySelectorAll('#t th[data-k]').forEach(th => th.onclick = () => {
  const k = th.dataset.k;
  if (k === sortK) { asc = !asc; } else { sortK = k; asc = !['contact','cross_min','pre_out','post_in'].includes(k); }
  view();
});
['fg','fs','fd','fr','q'].forEach(id => document.getElementById(id).addEventListener('input', view));
view();
</script>
</main></body></html>"""

DATASET_TEMPLATE = """<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__ — 結果</title><style>""" + BASE_CSS + """
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px}.stat{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px}
.stat b{display:block;font-size:20px}.controls{display:flex;flex-wrap:wrap;gap:12px 20px;align-items:center;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px;margin:10px 0}
.controls label{display:flex;gap:6px;align-items:center}.controls label.chip{display:inline-flex;gap:3px;padding:2px 8px;cursor:pointer}select,input{font:inherit;padding:4px 6px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg)}
.tablewrap{overflow-x:auto;background:var(--card);border:1px solid var(--line);border-radius:10px}table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap;vertical-align:top}th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none}
th.sorted::after{content:" ▼"}th.sorted.asc::after{content:" ▲"}td.num{text-align:right;font-variant-numeric:tabular-nums}
.chip{display:inline-block;padding:0 6px;border-radius:999px;background:var(--chip);margin-right:3px;font-size:12px}.chip.on{background:var(--accent);color:#fff}
.cell{font-weight:600}.sub{color:var(--muted);font-size:12px}.drv{color:var(--good);font-size:12px}.drv.no{color:var(--warn)}
.note{background:var(--card);border-left:3px solid var(--warn);padding:8px 12px;border-radius:6px;margin:10px 0}
.methods td:first-child{font-weight:700}.mat td,.mat th{text-align:right;cursor:default}
</style></head><body><main>
<p><a href="../index.html">← 結果一覧</a> ／ <a href="../shortlist.html">★ 最終候補（3 個体で確かめた組）</a></p>
<p class="muted">このページは<b>絞り込み前</b>の探索結果です。実験に進める候補は「最終候補」を見てください。</p>
<h1 id="title"></h1><p class="muted" id="subtitle"></p>
<div class="grid" id="stats"></div>

<div id="tpsec" style="display:none">
<h2>細胞タイプの組（ドライバーで実験できる単位）</h2>
<div class="controls">
  <label>接触の種類 <select id="tpkind"><option value="output">出力部位→相手（G、有向）</option><option value="directed">軸索→樹状突起（C、有向）</option><option value="undirected">中心線の接触（A、無向）</option></select></label>
  <label><input type="checkbox" id="tpff"> 上流→下流のみ</label>
  <label><input type="checkbox" id="tpmem"> 記憶系（キノコ体）を含む</label>
  <label><input type="checkbox" id="tpflag" checked> flagged のみ</label>
  <label><input type="checkbox" id="tpcross" checked> 異なるグループ間のみ</label>
  <label>タイプ <input id="tpq" placeholder="例: LC10" size="10"></label>
  <span class="muted" id="tpcount"></span>
</div>
<div class="tablewrap"><table id="tptbl"><thead><tr><th>タイプ 1</th><th>タイプ 2</th><th>細胞数</th><th>接触ペア / 結合ペア</th><th>接触合計 µm</th><th>シナプス</th><th>期待値</th><th>観測/期待</th><th>ドライバー</th><th>代表ペア（最も長く接し、結合が弱い）</th></tr></thead><tbody></tbody></table></div>
<p class="muted">期待値 = データ全体の「シナプス数 / 接触長」× その組の接触長。flagged = 期待値 10 以上で、観測が期待値の 0.2 倍以下。
有向（C）では「タイプ 1 の軸索 → タイプ 2 の樹状突起」、G では「タイプ 1 の出力部位（T-bar 付近）→ タイプ 2 のどこか」。シナプス数は接触していないペアも含めた、その組の全ペアの合計。</p>
</div>

<h2>ペアの一覧</h2>
<div class="controls">
  <label>「上位」の基準 <select id="pct"><option value="0.5">上位 0.5%</option><option value="1" selected>上位 1%</option><option value="2">上位 2%</option><option value="5">上位 5%</option></select></label>
  <label>一致したアプローチ数 ≥ <select id="minm"></select></label>
  <span style="display:flex;gap:6px;align-items:center">必ず含むアプローチ <span id="must" style="display:inline-flex;gap:4px;flex-wrap:wrap"></span></span>
  <label>細胞タイプ/ID <input id="q" placeholder="例: MBON01" size="14"></label>
  <label><input type="checkbox" id="only3d"> 3D ありのみ</label>
  <span class="muted" id="count"></span>
</div>
<div class="tablewrap"><table id="tbl"><thead></thead><tbody></tbody></table></div>
<p class="muted">順位の数字は「各アプローチの中で上位何 %」（小さいほど怪しい）。色付きのチップ = 選んだ基準で上位に入ったアプローチ。
w(1→2) / w(2→1) はシナプス数。3D がないペアは <code>python view_pair.py --data data/__TITLE__.pkl --pair ID1 ID2</code> で作れます。</p>

<h2>細胞タイプとドライバー</h2>
<div class="note">ドライバーは <b>細胞タイプ単位</b>。同じタイプのニューロンどうしのペア（例: KCg-m × KCg-m）は、遺伝学的に個別に操作できない点に注意。
ドライバーの情報は <code>drivers.csv</code> に追記すると反映されます。「要確認」は文献の原表で未確認のもの。</div>
<div class="tablewrap"><table id="drv"><thead><tr><th>細胞タイプ</th><th>ドライバー</th><th>種類</th><th>発現</th><th>確認</th><th>出典</th></tr></thead><tbody></tbody></table></div>

<h2>アプローチ</h2>
<div class="tablewrap"><table class="methods" id="meth"><thead><tr><th></th><th>名前</th><th>内容</th><th>ペア数</th></tr></thead><tbody></tbody></table></div>
<div class="note">C〜F は同じ「軸索–樹状突起の接触」計算を土台にしているので互いに独立ではない。実質的に独立なのは「形だけの A・B」と「区画を考える C〜F」の 2 系統。</div>
<h2>アプローチ間の重なり（選んだ基準で上位に入ったペア数）</h2>
<div class="tablewrap"><table class="mat" id="mat"></table></div>

<script>
const D = __DATA__;
const $ = s => document.querySelector(s);
const M = D.methods;
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const name = id => { const n = D.neurons[id] || {}; return n.instance || n.type || "?"; };
const drvTypes = new Set(D.drivers.map(d => d.type));
$("#title").textContent = D.name;
$("#subtitle").textContent = `${D.label_desc ? D.label_desc + " ／ " : ""}${D.dataset} ・ ${D.n_neurons} ニューロン ・ 結合 ${D.n_edges} 組 ・ データ ${D.stem}` + (D.has_instance ? "" : " ・ instance 名なし（fetch_data.py neuprint-meta で追加可）");
document.title = D.name + " — 結果";

let sortKey = "n", sortAsc = false;
function nTop(r, pct){ return M.filter(s => r.pct[s] !== null && r.pct[s] <= pct).length; }
function meanPct(r){ const v = M.map(s => r.pct[s]).filter(x => x !== null); return v.reduce((a,b)=>a+b,0)/v.length; }

$("#minm").innerHTML = M.map((_, i) => `<option value="${i+1}" ${i+1===Math.min(3,M.length)?"selected":""}>${i+1}</option>`).join("");
$("#must").innerHTML = M.map(s => `<label class="chip"><input type="checkbox" data-m="${s}"> ${s}</label>`).join("");

function stats(){
  const pct = +$("#pct").value;
  const all = D.rows.filter(r => nTop(r,pct) === M.length).length;
  const ge3 = D.rows.filter(r => nTop(r,pct) >= 3).length;
  const w3d = D.rows.filter(r => r.link).length;
  $("#stats").innerHTML = [
    ["全アプローチで上位 " + pct + "%", all + " 組"], ["3 つ以上で上位 " + pct + "%", ge3 + " 組"],
    ["一覧のペア数", D.rows.length + " 組"], ["3D 表示あり", w3d + " 組"]
  ].map(([k,v]) => `<div class="stat"><span class="muted">${k}</span><b>${v}</b></div>`).join("");
  const sets = Object.fromEntries(M.map(s => [s, new Set(D.rows.filter(r => r.pct[s] !== null && r.pct[s] <= pct).map(r => r.a+"-"+r.b))]));
  $("#mat").innerHTML = "<tr><th></th>" + M.map(s=>`<th>${s}</th>`).join("") + "</tr>" +
    M.map(x => `<tr><th>${x}</th>` + M.map(y => `<td>${[...sets[x]].filter(k => sets[y].has(k)).length}</td>`).join("") + "</tr>").join("");
}

const COLS = [
  ["n", "一致数"], ["pair", "ペア（細胞名 / bodyId）"], ["drv", "ドライバー"],
  ...M.map(s => ["m_"+s, s]), ["contact", "接触 µm (A)"], ["area", "面積 µm² (B)"], ["w12", "w(1→2)"], ["w21", "w(2→1)"], ["dir", "C の向き"], ["link", "3D"]
];
function head(){
  $("#tbl thead").innerHTML = "<tr>" + COLS.map(([k,l]) => `<th data-k="${k}" class="${k===sortKey?"sorted"+(sortAsc?" asc":""):""}" title="${k.startsWith("m_") ? esc(D.method_info[k.slice(2)][1]) : ""}">${l}</th>`).join("") + "</tr>";
  document.querySelectorAll("#tbl th").forEach(th => th.onclick = () => {
    const k = th.dataset.k; if (sortKey === k) sortAsc = !sortAsc; else { sortKey = k; sortAsc = !(k === "n"); }
    render();
  });
}
function val(r, k, pct){
  if (k === "n") return nTop(r,pct) * 1000 - meanPct(r);
  if (k.startsWith("m_")) { const v = r.pct[k.slice(2)]; return v === null ? 1e9 : v; }
  return {contact: r.contact_um ?? -1, area: r.contact_area_um2 ?? -1, w12: r.w12, w21: r.w21,
          pair: name(r.a)+name(r.b), dir: r.direction || "", link: r.link ? 0 : 1,
          drv: (drvTypes.has((D.neurons[r.a]||{}).type)?0:1) + (drvTypes.has((D.neurons[r.b]||{}).type)?0:1)}[k];
}
function cellHtml(id){
  const n = D.neurons[id] || {};
  const extra = [n.group, n.instance && n.type && n.instance !== n.type ? n.type : null, n.cropped ? "cropped" : null, n.status && n.status !== "Traced" ? n.status : null].filter(Boolean).join("・");
  return `<span class="cell">${esc(name(id))}</span> <span class="sub">${id}${extra ? "・"+esc(extra) : ""}</span>`;
}
function drvHtml(id){
  const t = (D.neurons[id]||{}).type; const ds = D.drivers.filter(d => d.type === t);
  return ds.length ? `<span class="drv">${esc(t)}: ${ds.map(d => esc(d.driver) + (d.status === "要確認" ? "?" : "")).join(", ")}</span>` : `<span class="drv no">${esc(t||"?")}: 未登録</span>`;
}
function render(){
  const pct = +$("#pct").value, minm = +$("#minm").value, q = $("#q").value.trim().toLowerCase(), only3d = $("#only3d").checked;
  const must = [...document.querySelectorAll("#must input:checked")].map(i => i.dataset.m);
  let rows = D.rows.filter(r => nTop(r,pct) >= minm && must.every(s => r.pct[s] !== null && r.pct[s] <= pct) && (!only3d || r.link));
  if (q) rows = rows.filter(r => [r.a, r.b, name(r.a), name(r.b), (D.neurons[r.a]||{}).type, (D.neurons[r.b]||{}).type].some(x => String(x||"").toLowerCase().includes(q)));
  rows.sort((x,y) => { const a = val(x,sortKey,pct), b = val(y,sortKey,pct); return (a<b?-1:a>b?1:0) * (sortAsc?1:-1); });
  $("#count").textContent = `${rows.length} 組` + (rows.length > 300 ? "（先頭 300 組を表示）" : "");
  head();
  $("#tbl tbody").innerHTML = rows.slice(0,300).map(r => {
    const n = nTop(r,pct);
    return `<tr><td class="num"><b>${n}</b>/${M.length}</td>
      <td>${cellHtml(r.a)}<br>${cellHtml(r.b)}</td><td>${drvHtml(r.a)}<br>${drvHtml(r.b)}</td>` +
      M.map(s => { const v = r.pct[s]; const on = v !== null && v <= pct;
        return `<td class="num" title="${s}: ${r.rank[s] ?? "-"} 位${r.flag[s] ? "・flagged" : ""}">${v === null ? "–" : `<span class="chip ${on?"on":""}">${v < 0.1 ? v.toFixed(3) : v.toFixed(1)}%</span>`}</td>`; }).join("") +
      `<td class="num">${r.contact_um ?? "–"}</td><td class="num">${r.contact_area_um2 ?? "–"}</td>
      <td class="num">${r.w12}</td><td class="num">${r.w21}</td><td class="sub">${esc(r.direction || "")}</td>
      <td>${r.link ? `<a href="${r.link}" target="_blank">3D を開く</a>` : '<span class="sub">なし</span>'}</td></tr>`;
  }).join("");
}
$("#drv tbody").innerHTML = Object.keys(D.lookup).map(t => {
  const ds = D.drivers.filter(d => d.type === t);
  const look = D.lookup[t].map(([l,u]) => `<a href="${u}" target="_blank">${esc(l)}</a>`).join(" ・ ");
  if (!ds.length) return `<tr><td class="cell">${esc(t)}</td><td colspan="4" class="drv no">未登録 — 調べる: ${look}</td><td></td></tr>`;
  return ds.map((d,i) => `<tr><td class="cell">${i===0?esc(t):""}</td><td><b>${esc(d.driver)}</b></td><td>${esc(d.system)}</td><td style="white-space:normal;min-width:260px">${esc(d.expression)}</td><td class="${d.status==="確認済み"?"drv":"drv no"}">${esc(d.status)}</td><td><a href="${esc(d.source)}" target="_blank">出典</a></td></tr>`).join("") +
    `<tr><td></td><td colspan="5" class="sub">ほかを調べる: ${look}</td></tr>`;
}).join("");
function tpRender(){
  const kind = $("#tpkind").value, rows0 = (D.type_pairs||{})[kind] || [];
  const q = $("#tpq").value.trim().toLowerCase();
  let rows = rows0.filter(r => (!$("#tpflag").checked || r.flagged) && (!$("#tpcross").checked || (r.g1 && r.g2 && r.g1 !== r.g2))
                         && (!$("#tpff").checked || r.feedforward) && (!$("#tpmem").checked || r.memory));
  if (q) rows = rows.filter(r => (r.t1+" "+r.t2).toLowerCase().includes(q));
  $("#tpcount").textContent = rows.length + " 組";
  const dr = t => { const ds = D.drivers.filter(d => d.type === t); return ds.length ? `<span class="drv">${ds.map(d => esc(d.driver)+(d.status==="要確認"?"?":"")).join(", ")}</span>` : `<a class="drv no" href="https://www.google.com/search?q=${encodeURIComponent(t+" split-GAL4 Drosophila")}" target="_blank">調べる</a>`; };
  const g = x => x ? `<span class="sub">${esc(x)}</span>` : "";
  $("#tptbl tbody").innerHTML = rows.slice(0,300).map(r => `<tr>
    <td class="cell">${esc(r.t1)} ${g(r.g1)}</td><td class="cell">${esc(r.t2)} ${g(r.g2)}${r.feedforward ? ' <span class="chip">上流→下流</span>' : ''}${r.memory ? ' <span class="chip">記憶系</span>' : ''}</td>
    <td class="num">${r.n_cells_1} × ${r.n_cells_2}</td><td class="num">${r.n_pairs_contact} / ${r.n_pairs_connected}</td>
    <td class="num">${(r.contact_total_um??0).toFixed(1)}</td><td class="num">${r.synapses}</td><td class="num">${(r.expected??0).toFixed(1)}</td>
    <td class="num"><span class="chip ${r.flagged?"on":""}">${(r.obs_over_exp??0).toFixed(2)}</span></td>
    <td>${esc(r.t1)}: ${dr(r.t1)}<br>${esc(r.t2)}: ${dr(r.t2)}</td>
    <td>${r.ex_a ? `<span class="sub">${r.ex_a} – ${r.ex_b}（${(r.ex_contact??0).toFixed(1)} µm, ${r.ex_weight} シナプス）</span><br>` + (r.link ? `<a href="${r.link}" target="_blank">3D を開く</a>` : '<span class="sub">3D なし</span>') : '<span class="sub">–</span>'}</td></tr>`).join("");
}
if (D.type_pairs && Object.keys(D.type_pairs).length) {
  $("#tpsec").style.display = "";
  // flagged が多い方の接触の種類を最初に表示する
  const nf = k => ((D.type_pairs[k]||[]).filter(r => r.flagged && r.g1 && r.g2 && r.g1 !== r.g2)).length;
  const kinds = ["output", "directed", "undirected"].filter(k => D.type_pairs[k]);
  $("#tpkind").value = kinds.sort((a, b) => nf(b) - nf(a))[0];
  [...$("#tpkind").options].forEach(o => { if (!D.type_pairs[o.value]) o.disabled = true; });
  ["#tpkind","#tpflag","#tpcross","#tpq","#tpff","#tpmem"].forEach(s => $(s).addEventListener("input", tpRender));
  tpRender();
}
$("#meth tbody").innerHTML = M.map(s => `<tr><td>${s}</td><td>${esc(D.method_info[s][0])}</td><td style="white-space:normal">${esc(D.method_info[s][1])}</td><td class="num">${D.sizes[s]}</td></tr>`).join("");
["#pct","#minm","#q","#only3d"].forEach(s => $(s).addEventListener("input", () => { stats(); render(); }));
document.querySelectorAll("#must input").forEach(i => i.addEventListener("change", render));
stats(); render();
</script>
</main></body></html>"""


if __name__ == "__main__":
    main()
