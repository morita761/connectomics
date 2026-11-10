"""人工分子でシナプスを作らせる実験に向けた最終候補を、複数の個体（データセット）を合わせて作る。

  python shortlist.py --max-total 2

条件（強制シナプス形成の実験を想定）:
  1. 膜が接している（接着分子が働くには向かい合っている必要がある）
  2. 全データセット・両方向でシナプスがほぼ 0 → 新しくできたシナプスを 0 を基準に数えられる
  3. 別個体でも接触が再現する（接触が偶然でないこと）
  4. 両側にドライバーがある（前側・後側の分子を別々に入れるため）

元データ（hemibrain）に加えて、validate_cross.py で作った各データセットの
cross_<データセット>.csv・cross_<データセット>_contact.csv をすべて読み、1 行にまとめる。
照合できたデータセットが多い組（n_datasets）を優先する。

出力: results/shortlist.csv と results/shortlist.md
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common.runner import ROOT

REGIONS = {"mb_all": "キノコ体", "lh": "外側角", "cx": "中心複合体", "vis": "視覚", "al": "触角葉",
           "smp": "上前大脳", "lal": "側方付属葉", "wed": "楔形領域", "ps": "後斜面",
           "cre": "クランプ", "vlp": "腹外側前大脳", "lexa": "LexA 12 型を軸に再探索", "mb_two": "キノコ体の中だけ"}

# 下流タイプに行動・生理の読み出しがあるか（lit_search.py と Dolan 2019 補足表 3 で確認した分）
READOUTS = {
    "LHAV1a1": "活性化で忌避・抑制でジオスミン応答消失（Dolan 2019）",
    "LHAV1a3": "光遺伝学で検査済みだが有意差なし（Dolan 2019）",
    "LHPD2a1": "学習した匂いの想起に必要（Dolan 2018）",
    "LHPD2b1": "同上",
    "LHPV5a1": "活性化で誘引（PI +0.38）",
    "LHPV6a3": "活性化で誘引（PI +0.42）",
    "LHAD1b2": "活性化で誘引（PI +0.53）",
    "LHPV4a1": "有意な表現型あり", "LHPV4a2": "有意な表現型あり",
    "PFL2": "目標指向の操舵（Westeinde 2024）", "PFL3": "目標指向の操舵（Westeinde 2024）",
    "MBON01": "活性化で忌避（Owald 2015・Aso 2014。MB011B で 3 タイプまとめて操作）",
    "MBON04": "活性化で忌避（MB011B で MBON01・03・04 をまとめて操作した実験のみ）",
    "MBON03": "活性化で忌避（MB011B で MBON01・03・04 をまとめて操作した実験のみ）",
    "MBON11": "活性化で誘引・抑制で 2 時間の嫌悪記憶が低下（Aso 2014・Perisse 2016）",
    # 下行ニューロン: 行動の読み出しが直接とれる
    "DNa02": "活性化で同側への旋回。歩行中の角速度と相関（Rayshubskiy 2024）",
    "DNa01": "角速度と相関し、活性化で旋回（Rayshubskiy 2024）",
    "DNp09": "活性化で歩行停止・フリージング（Zacarias 2018）",
    "DNb01": "光遺伝学で行動の変化を記述済み（Cande 2018、系統 SS02542）",
    # MBON（lit_search.py で確認）
    "MBON30": "活動低下で嫌悪記憶の表出・抑制で食物報酬記憶が低下。睡眠の促進にも関与"
              "（Aso 2014・Sitaraman 2015・Jacob 2021）",
    # 読み出しが無いことが分かったもの（記録として残す）
    "oviIN": "論文 27 本を検索したが行動実験の記述は見つからなかった（形態の記載のみ）",
    "FB4N": "論文 3 本のみ。行動実験の記述なし",
    "FB4X": "論文 3 本のみ。行動実験の記述なし",
}


def drivers() -> pd.DataFrame:
    d = pd.concat([pd.read_csv(f) for f in sorted((ROOT / "data").glob("driver_lines*.csv"))], ignore_index=True)
    d["type"] = d.em_types.astype(str).str.split(";")
    d = d.explode("type")
    d["n_types"] = pd.to_numeric(d.n_types, errors="coerce").fillna(1)
    return d


def is_lexa(system) -> bool:
    return "lexa" in str(system).lower()


def best(d: pd.DataFrame, t: str, k: int = 2, want_lexa: bool | None = None) -> str:
    """そのタイプの系統を、標識タイプ数が少なく特異性が高い順に。want_lexa で LexA 系 / GAL4 系を絞る。"""
    h = d[d.type == t]
    if want_lexa is not None:
        h = h[h.system.apply(is_lexa) == want_lexa]
    if h.empty:
        return "-"
    h = h.assign(_s=h.specificity.fillna(-1)).sort_values(["n_types", "_s"], ascending=[True, False])
    out = []
    for r in h.head(k).itertuples(index=False):
        extra = [x for x in str(r.em_types).split(";") if x != t]
        sysname = "LexA" if is_lexa(r.system) else "GAL4"
        out.append(f"{r.line}〈{sysname}〉" + (f"(+{','.join(extra)})" if extra else "")
                   + ("" if pd.isna(r.specificity) else f"[特異性{int(r.specificity)}]"))
    return " / ".join(out)


def collect(region: str) -> pd.DataFrame | None:
    """1 つの領域について、全データセットの結合・接触をタイプの組ごとに 1 行にまとめる。"""
    # 領域ごとのデータは results/hemibrain_<領域>/、領域をまたぐ探索は results/<名前>/
    res = ROOT / "results" / f"hemibrain_{region}"
    if not res.exists():
        res = ROOT / "results" / region
    conn_files = [f for f in sorted(res.glob("cross_*.csv")) if not f.name.endswith("_contact.csv")]
    if not conn_files:
        return None
    conn = pd.concat([pd.read_csv(f) for f in conn_files], ignore_index=True)
    cont_parts = []
    for f in sorted(res.glob("cross_*_contact.csv")):
        c = pd.read_csv(f)
        if "dataset" not in c:  # 古い形式のファイルはファイル名からデータセット名を補う
            c["dataset"] = f.name[len("cross_"):-len("_contact.csv")].replace("_", ":", 1)
        cont_parts.append(c)
    cont = pd.concat(cont_parts, ignore_index=True) if cont_parts else None

    rows = []
    for (t1, t2), g in conn.groupby(["t1", "t2"]):
        ok = g[g.n_sides > 0]
        src = g.iloc[0]
        row = {
            "region": REGIONS.get(region, region), "t1": t1, "t2": t2,
            "contact_src_um": src.contact_um_src, "expected_src": src.expected_src,
            "syn_src_fwd": src.syn_src_fwd, "syn_src_rev": src.syn_src_rev,
            "n_datasets": len(ok),  # hemibrain 以外で照合できた個体の数
            "datasets": ", ".join(f"{r.dataset}({r.detail})" for r in ok.itertuples(index=False)),
            "missing": ", ".join(g[g.n_sides == 0].dataset),
            "syn_cross": int(ok.cross_total.sum()) if len(ok) else -1,
            "comp_type": src.comp_type, "comp_rate": src.comp_rate,
        }
        row["total_syn"] = row["syn_src_fwd"] + row["syn_src_rev"] + max(row["syn_cross"], 0)
        if cont is not None:
            c = cont[(cont.t1 == t1) & (cont.t2 == t2)]
            row["cross_contact_detail"] = ", ".join(
                f"{r.dataset.split(':')[0]} {r.side} {r.contact_um:.0f}µm" for r in c.itertuples(index=False))
            row["cross_contact_min_um"] = float(c.contact_um.min()) if len(c) else np.nan
            row["cross_contact_max_um"] = float(c.contact_um.max()) if len(c) else np.nan
            row["n_sides_with_contact"] = int(len(c))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-total", type=int, default=2, help="全データセット・両方向の合計シナプス数の上限")
    ap.add_argument("--min-cross-contact", type=float, default=10.0,
                    help="別個体での接触長の下限 [µm]（これ以上を「再現」とみなす）")
    args = ap.parse_args()

    dr = drivers()
    parts = [p for r in REGIONS if (p := collect(r)) is not None]
    if not parts:
        print("cross_*.csv が見つかりません（validate_cross.py を先に実行）")
        return
    out = pd.concat(parts, ignore_index=True)
    # 同じ組が複数の領域の探索に出てくることがある（領域が隣接している・LexA 軸の再探索と重なる）。
    # 1 組 1 行にまとめ、領域名は「/」でつなぐ
    label = out.groupby(["t1", "t2"]).region.apply(lambda v: " / ".join(dict.fromkeys(v)))
    out = (out.sort_values(["n_datasets", "contact_src_um"], ascending=False)
           .drop_duplicates(["t1", "t2"]))
    out["region"] = [label[(a, b)] for a, b in zip(out.t1, out.t2)]
    out["up_driver"] = [best(dr, t) for t in out.t1]
    out["down_driver"] = [best(dr, t) for t in out.t2]
    out["both_drivers"] = (out.up_driver != "-") & (out.down_driver != "-")
    # 上流と下流に別々の分子を入れるには 2 つの発現系が要る。上流 GAL4 + 下流 LexA が理想。
    out["up_gal4"] = [best(dr, t, want_lexa=False) for t in out.t1]
    out["down_lexa"] = [best(dr, t, want_lexa=True) for t in out.t2]
    out["up_lexa"] = [best(dr, t, want_lexa=True) for t in out.t1]
    out["down_gal4"] = [best(dr, t, want_lexa=False) for t in out.t2]
    out["two_systems"] = ((out.up_gal4 != "-") & (out.down_lexa != "-")) | \
                         ((out.up_lexa != "-") & (out.down_gal4 != "-"))
    out["readout"] = out.t2.map(READOUTS).fillna("-")
    if "cross_contact_min_um" in out:
        out["contact_reproduced"] = out.cross_contact_min_um >= args.min_cross_contact
    else:
        out["contact_reproduced"] = False
    out = out[(out.total_syn <= args.max_total) & (out.n_datasets > 0)]
    # 並べ方: 接触が再現 → 照合できた個体数 → 両側ドライバー → シナプス数 → 元の接触長
    out = out.sort_values(["contact_reproduced", "n_datasets", "both_drivers", "total_syn", "contact_src_um"],
                          ascending=[False, False, False, True, False])
    f = ROOT / "results" / "shortlist.csv"
    out.to_csv(f, index=False, encoding="utf-8")

    # 接触領域の情報（inspect_pairs.py）を足して、最終候補 shortlist_final.csv にする
    ins = [pd.read_csv(p) for p in ROOT.glob("results/hemibrain_*/inspect_output.csv")]
    if ins:
        ins = pd.concat(ins, ignore_index=True).drop_duplicates(["t1", "t2"])
        fin = out.merge(ins[["t1", "t2", "contact_rois", "main_roi",
                             "pre_out_in_main", "post_in_in_main"]], on=["t1", "t2"], how="left")
        # 上流がその領域で出力し、下流が入力を受けているか（通り道が重なっただけでないか）
        fin["ok_active"] = (fin.pre_out_in_main.fillna(0) >= 100) & (fin.post_in_in_main.fillna(0) >= 100)
        fin["score"] = (fin.both_drivers.astype(int) * 2 + fin.two_systems.astype(int) * 2
                        + fin.ok_active.astype(int) * 2 + (fin.cross_contact_min_um >= 20).astype(int))
        fin = fin.sort_values(["score", "contact_src_um"], ascending=False)
        g = ROOT / "results" / "shortlist_final.csv"
        fin.to_csv(g, index=False, encoding="utf-8")
        print(f"最終候補 {len(fin)} 組（2 つの発現系がそろう組 {int(fin.two_systems.sum())}）-> {g}")

    md = ["# 強制シナプス形成に向けた最終候補（shortlist.py）", "",
          f"条件: 合計シナプス {args.max_total} 個以下。「接触の再現」は別個体の各半球で "
          f"{args.min_cross_contact:.0f} µm 以上の接触があること。",
          "シナプス数は hemibrain の順/逆と、照合できた各個体の両方向の合計。", "",
          "| 領域 | 上流 → 下流 | hemibrain 接触 | シナプス合計 | 照合した個体 | 別個体の接触 | 上流ドライバー | 下流ドライバー |",
          "|---|---|---|---|---|---|---|---|"]
    for r in out.itertuples(index=False):
        md.append(f"| {r.region} | {r.t1} → {r.t2} | {r.contact_src_um:.0f} µm | {int(r.total_syn)} | "
                  f"{r.datasets} | {getattr(r, 'cross_contact_detail', '') or '未確認'} | "
                  f"{r.up_driver} | {r.down_driver} |")
    (ROOT / "results" / "shortlist.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    cols = ["region", "t1", "t2", "contact_src_um", "total_syn", "n_datasets",
            "cross_contact_min_um", "cross_contact_max_um", "up_driver", "down_driver"]
    with pd.option_context("display.width", 260, "display.max_colwidth", 30):
        print(out[[c for c in cols if c in out]].round(1).to_string(index=False))
    print(f"\n{len(out)} 組 -> {f} / {ROOT / 'results' / 'shortlist.md'}")


if __name__ == "__main__":
    main()
