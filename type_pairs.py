"""ニューロンのペアの結果を「細胞タイプの組」ごとにまとめる（ドライバーで実験できる単位）。

  python type_pairs.py --data data/hemibrain_vis_olf.pkl --cross visual olfactory
  python type_pairs.py --data data/hemibrain_mb.pkl            # グループ指定なし = 全タイプの組

入力: run_all.py が results/<データ名>/ に出力した A（スケルトン接触）・C（軸索→樹状突起接触）・
      G（出力部位→相手）の結果。
      B（表面接触）があれば接触面積も足す。

タイプの組ごとに:
  contact_um      その組のニューロンペアの接触長の合計（A: 無向 / C: 軸索→樹状突起の有向）
  synapses        その組のニューロンペアのシナプス数の合計（接触していないペアも含む全ペア）
  expected        基準の密度 × その組の接触長。基準の密度 = 表に残った組全体の「シナプス数 / 接触長」
                  （--cross を付けるとグループ間の組だけで計算する）
  obs_over_exp    synapses / expected（小さいほど「接触のわりに結合しない」）
  n_pairs_contact 接触しているニューロンペア数、n_pairs_connected そのうち weak_synapse_max を超えて結合している数
flagged: 期待値が min_type_expected 以上なのに、観測が期待値の max_obs_over_exp 倍以下
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from common.config import Thresholds, parse_overrides
from common.model import Dataset
from common.runner import ROOT


def summarize(df: pd.DataFrame, ds: Dataset, a: str, b: str, directed: bool, th: Thresholds,
              contact_col: str, min_expected: float, max_ratio: float, cross: list[str] | None) -> pd.DataFrame:
    meta = {i: n.meta for i, n in ds.neurons.items()}
    t = lambda i: meta[i].get("type") or f"untyped_{i}"
    g = lambda i: meta[i].get("group")
    df = df[df[contact_col] >= th.min_contact_cable_um].copy() if contact_col == "contact_um" else df.copy()
    df["t1"], df["t2"] = df[a].map(t), df[b].map(t)
    df["g1"], df["g2"] = df[a].map(g), df[b].map(g)
    if not directed:  # 無向: タイプ名の順に並べる
        sw = df.t1 > df.t2
        for x, y in (("t1", "t2"), ("g1", "g2")):
            df.loc[sw, [x, y]] = df.loc[sw, [y, x]].values
    df["connected"] = df.weight > th.weak_synapse_max
    # 代表ペア: その組の中で、結合が弱い（weak 以下）ペアのうち最も長く接しているもの（3D 確認用）
    weak = df[~df.connected]
    ex = weak.loc[weak.groupby(["t1", "t2"])[contact_col].idxmax()].set_index(["t1", "t2"])
    agg = df.groupby(["t1", "t2"]).agg(g1=("g1", "first"), g2=("g2", "first"),
                                        n_pairs_contact=(contact_col, "size"),
                                        n_pairs_connected=("connected", "sum"),
                                        contact=(contact_col, "sum"))

    # シナプス数は全ペアから（接触していないペアのシナプスも含める）
    e = ds.edges.assign(t1=ds.edges.pre.map(t), t2=ds.edges.post.map(t))
    if not directed:
        sw = e.t1 > e.t2
        e.loc[sw, ["t1", "t2"]] = e.loc[sw, ["t2", "t1"]].values
    agg["ex_a"] = ex[a].reindex(agg.index)
    agg["ex_b"] = ex[b].reindex(agg.index)
    agg["ex_contact"] = ex[contact_col].reindex(agg.index)
    agg["ex_weight"] = ex["weight"].reindex(agg.index)
    agg["synapses"] = e.groupby(["t1", "t2"]).weight.sum().reindex(agg.index).fillna(0).astype(int)
    n_types = {k: v for k, v in pd.Series({i: t(i) for i in meta}).value_counts().items()}
    agg = agg.reset_index()
    agg["n_cells_1"] = agg.t1.map(n_types)
    agg["n_cells_2"] = agg.t2.map(n_types)
    if cross:
        x, y = cross
        keep = (agg.g1 == x) & (agg.g2 == y)
        if not directed:
            keep |= (agg.g1 == y) & (agg.g2 == x)
        agg = agg[keep].copy()
    # 期待値の基準: 表に残った組（--cross なら グループ間の組だけ）の「シナプス数 / 接触長」。
    # 同じ系どうしの接触は束になって並走する軸索が多く、密度が低く出るので、それを基準にしない。
    rate = agg.synapses.sum() / agg.contact.sum()
    agg["baseline_rate_per_um"] = rate
    agg["expected"] = rate * agg.contact
    agg["obs_over_exp"] = agg.synapses / agg.expected
    agg["score"] = np.log2((agg.expected + 1) / (agg.synapses + 1))
    agg["p_poisson"] = stats.poisson.cdf(agg.synapses, agg.expected)
    agg["flagged"] = (agg.expected >= min_expected) & (agg.obs_over_exp <= max_ratio)
    agg = agg.rename(columns={"contact": contact_col.replace("contact", "contact_total")})
    agg.attrs["rate"] = rate
    return agg.sort_values("score", ascending=False).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--cross", nargs=2, default=None, metavar=("GROUP1", "GROUP2"),
                    help="この 2 グループの間の組だけ出す（fetch_data.py --group で付けた名前）")
    ap.add_argument("--min-type-expected", type=float, default=5.0)
    ap.add_argument("--max-obs-over-exp", type=float, default=0.2)
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--levels", nargs="*", default=["PN=0", "KC=1", "DAN=1", "MBIN=1", "MBON=2", "LHN=2"],
                    metavar="GROUP=LEVEL", help="グループの階層（小さいほど上流）。上流→下流の組を優先する")
    ap.add_argument("--memory-groups", nargs="*", default=["KC", "MBON", "DAN", "MBIN"],
                    help="記憶に関わるグループ（キノコ体）。これを含む組を優先する")
    args = ap.parse_args()

    th = Thresholds().update(parse_overrides(args.set))
    ds = Dataset.load(args.data)
    res = ROOT / "results" / Path(args.data).stem
    # グループ（fetch_data.py --group）があれば、グループの組ごとに分けて集計する。
    # 同じ系どうし（束になって並走する軸索が多い）と異なる系の間では接触あたりの結合密度が大きく違うため。
    groups = sorted({n.meta.get("group") for n in ds.neurons.values() if n.meta.get("group")})
    if args.cross:
        combos = [args.cross]
    elif groups:
        combos = [[x, y] for k, x in enumerate(groups) for y in groups[k:]]
    else:
        combos = [None]
    # 有向の集計は、グループの順序付きの組ごと（例: KC→MBON と MBON→KC は別の基準）
    combos_dir = [[x, y] for x in groups for y in groups] if groups and not args.cross else combos
    levels = {k: float(v) for k, v in (x.split("=") for x in args.levels)}
    memory = set(args.memory_groups)
    outs = {}
    for key, name, a, b, directed in (("undirected", "approach_a_skeleton_kdtree", "a", "b", False),
                                      ("directed", "approach_c_axon_dendrite", "pre", "post", True),
                                      ("output", "approach_g_output_sites", "pre", "post", True)):
        f = res / f"{name}.csv"
        if not f.exists():
            print(f"skip {key}: {f} がありません（run_all.py を先に実行）")
            continue
        df = pd.read_csv(f)
        bdf = (pd.read_csv(res / "approach_b_surface_contact.csv")
               if key == "undirected" and (res / "approach_b_surface_contact.csv").exists() else None)
        parts = []
        for combo in (combos_dir if directed else combos):  # グループの組ごとに、それぞれの基準の密度で期待値を計算する
            t = summarize(df, ds, a, b, directed, th, "contact_um",
                          args.min_type_expected, args.max_obs_over_exp, combo)
            if bdf is not None:
                bt = summarize(bdf, ds, "a", "b", False, th, "contact_area_um2", 0, 1, combo)
                t = t.merge(bt[["t1", "t2", "contact_total_area_um2"]], on=["t1", "t2"], how="left")
            parts.append(t)
        tab = pd.concat(parts, ignore_index=True)
        # 優先度: 上流→下流（階層が下がる向き）と、記憶系（キノコ体）を含む組を優先
        tab["level1"] = tab.g1.map(levels)
        tab["level2"] = tab.g2.map(levels)
        tab["feedforward"] = directed & (tab.level1 < tab.level2)
        tab["memory"] = tab.g1.isin(memory) | tab.g2.isin(memory)
        tab["priority"] = tab.score + tab.feedforward.astype(int) + tab.memory.astype(int)
        tab = tab.sort_values(["flagged", "priority"], ascending=False).reset_index(drop=True)
        tab.to_csv(res / f"type_pairs_{key}.csv", index=False)
        outs[key] = tab
        label = {"undirected": "接触（無向）", "directed": "軸索→樹状突起（有向）",
                 "output": "出力部位→相手（有向）"}[key]
        print(f"\n=== タイプの組: {label}  {len(tab)} 組、flagged {int(tab.flagged.sum())} 組 "
              f"-> {res / f'type_pairs_{key}.csv'}")
        if groups:
            print(tab.groupby(["g1", "g2"]).agg(n=("flagged", "size"), flagged=("flagged", "sum"),
                                                baseline_per_um=("baseline_rate_per_um", "first"))
                  .round(3).to_string())
        cols = ["t1", "t2", "g1", "g2", "n_cells_1", "n_cells_2", "n_pairs_contact", "n_pairs_connected",
                "contact_total_um", "synapses", "expected", "obs_over_exp", "p_poisson", "feedforward", "memory", "flagged"]
        with pd.option_context("display.width", 250, "display.max_columns", 30):
            print(tab[cols].head(args.top).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
