"""候補のタイプの組が、別の個体のコネクトームでも「接触するのに結合しない」かを確かめる。

人工分子でシナプスを作らせる実験に向けて、「両方向ともシナプスが 1 個もない」組を探す用途を想定している。

  # 1) 結合だけを確かめる（スケルトン不要。速い）
  python validate_cross.py --data data/hemibrain_cx.pkl --datasets male-cns:v1.0 optic-lobe:v1.1 --top 150
  # 2) 残った組について、別個体でも接触しているかを確かめる（スケルトンを取得する。遅い）
  python validate_cross.py --data data/hemibrain_cx.pkl --stage contact --top 10 --max-total 2

段階 1（conn）:
  候補（candidates_ranked.csv のうち、元データで t1→t2 のシナプスが 0 の組）について、
  各データセットで両方向・半球ごとのシナプス数を数える（neuPrint の fetch_adjacencies）。
  元データの逆向き（t2→t1）も数え、合計を total_syn に入れる。
段階 2（contact）:
  段階 1 を通った組について、各データセットの代表ニューロンのスケルトンを取得し、
  アプローチ G（出力部位 → 相手）の接触長を計算する。接触があり、かつシナプスが 0 なら再現とみなす。

半球は somaSide、無ければ instance の末尾 _R / _L を使う。どちらも無いデータセット（optic-lobe など）は
半球を分けず "all" として 1 つにまとめる。

出力: results/<データ名>/cross_<データセット名>.csv と cross_<データセット名>_contact.csv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common.config import Thresholds, parse_overrides
from common.model import Dataset
from common.runner import ROOT, save_merged

UM = 1000.0
DEFAULT_DATASETS = ["male-cns:v1.0", "optic-lobe:v1.1"]


def tag_of(dataset: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", dataset)


from common.loaders import _retry  # 接続エラー・タイムアウトに対する再試行


def side_ids(client, types: list[str]) -> tuple[dict[tuple[str, str], list[int]], dict[str, int]]:
    """タイプ名 → {(タイプ, 半球): bodyId の一覧} と、タイプごとの注釈済みシナプス数（pre+post）。

    半球が分からないデータセットは "all"。
    シナプス数（pre+post）は参考値。実際に照合できるかは visible_weight() で判定する。
    """
    from neuprint import NeuronCriteria as NC, fetch_neurons
    df, _ = _retry(fetch_neurons, NC(type=types), client=client)
    out: dict[tuple[str, str], list[int]] = {}
    if not len(df):
        return out, {}
    nsyn = (df.assign(_n=df.get("pre", 0).fillna(0) + df.get("post", 0).fillna(0))
            .groupby(df.type.astype(str))._n.sum().astype(int).to_dict())  # 参考値（注釈数）
    if "somaSide" in df.columns:
        side = df.somaSide.astype(str).replace({"nan": "?", "None": "?", "": "?"})
    elif "instance" in df.columns:
        side = df.instance.astype(str).str.extract(r"_([RL])$")[0].fillna("?")
    else:
        side = pd.Series("?", index=df.index)
    if not side.isin(["R", "L"]).any():
        side = pd.Series("all", index=df.index)
    for (t, s), g in df.groupby([df.type.astype(str), side]):
        if s in ("R", "L", "all"):
            out[(t, s)] = [int(b) for b in g.bodyId]
    return out, nsyn


def visible_weight(client, ids: list[int]) -> tuple[int, int]:
    """そのニューロン群について、このデータセットで「見える」結合の重みの合計（出力, 入力）。

    optic-lobe のように一部の領域だけを公開したデータセットでは、中枢脳での結合が 0 件に見える。
    「もし結合があれば見えるはずか」を確かめないと、見かけ上の 0 を本物の 0 と取り違える。
    """
    from neuprint import NeuronCriteria as NC, fetch_adjacencies
    if not ids:
        return 0, 0
    out = in_ = 0
    _, c1 = _retry(fetch_adjacencies, NC(bodyId=ids), None, client=client)
    out = int(c1.weight.sum()) if len(c1) else 0
    _, c2 = _retry(fetch_adjacencies, None, NC(bodyId=ids), client=client)
    in_ = int(c2.weight.sum()) if len(c2) else 0
    return out, in_


def conn_count(client, ids_a: list[int], ids_b: list[int]) -> int:
    """ids_a → ids_b のシナプス数の合計（どちらかが空なら -1 = 確認できない）。"""
    from neuprint import NeuronCriteria as NC, fetch_adjacencies
    if not ids_a or not ids_b:
        return -1
    _, conn = _retry(fetch_adjacencies, NC(bodyId=ids_a), NC(bodyId=ids_b), client=client)
    return int(conn.weight.sum()) if len(conn) else 0


def reverse_in_source(ds: Dataset) -> dict[tuple[str, str], int]:
    """元データでのタイプの組ごとのシナプス数（逆向きを引くために使う）。"""
    types = {i: str(n.meta.get("type")) for i, n in ds.neurons.items()}
    rev: dict[tuple[str, str], int] = {}
    for (p, q), v in ds.directed_weights().items():
        k = (types[p], types[q])
        rev[k] = rev.get(k, 0) + v
    return rev


def stage_conn(cand: pd.DataFrame, rev: dict, dataset: str, top: int, min_syn: int) -> pd.DataFrame:
    from common.loaders import neuprint_client
    client = neuprint_client(dataset)
    cand = cand.head(top)
    ids, nsyn = side_ids(client, sorted(set(cand.t1) | set(cand.t2)))
    vis: dict[str, tuple[int, int]] = {}
    rows = []
    for r in cand.itertuples(index=False):
        sides = sorted({s for (t, s) in ids if t in (r.t1, r.t2)})
        det, fwd, back, n_sides = [], 0, 0, 0
        if not {t for (t, s) in ids} >= {r.t1, r.t2}:
            sides = []
        # もし結合があれば見えるはずか（上流に出力、下流に入力が見えているか）
        thin = []
        if not sides:
            thin = ["（タイプ名なし）"]
        else:
            for t, which in ((r.t1, 0), (r.t2, 1)):
                key = t
                if key not in vis:
                    vis[key] = visible_weight(client, [b for (tt, s), bs in ids.items() if tt == t for b in bs])
                if vis[key][which] < min_syn:
                    thin.append(f"{t} の{'出力' if which == 0 else '入力'}が見えない")
        if thin:
            rows.append({
                "t1": r.t1, "t2": r.t2, "dataset": dataset,
                "contact_um_src": r.contact_total_um, "expected_src": r.expected,
                "syn_src_fwd": int(r.synapses), "syn_src_rev": int(rev.get((r.t2, r.t1), 0)),
                "n_sides": 0, "syn_cross_fwd": -1, "syn_cross_rev": -1, "cross_total": -1,
                "detail": "（" + "・".join(thin) + "）",
                "total_syn": int(r.synapses) + int(rev.get((r.t2, r.t1), 0)),
                "comp_type": r.comp_type, "comp_rate": r.comp_rate,
                "availability": getattr(r, "availability", np.nan)})
            print(f"{r.t1:>14s} → {r.t2:<14s} {rows[-1]['detail']}")
            continue
        # 細胞体の左右で揃えて数えてはいけない。MBON01 のように細胞体と樹状突起が反対側にある細胞があり、
        # 揃えると本物の結合を見落とす。判定は「両半球をまとめた合計」で行い、左右別は参考として出す。
        all_a = [b for (t, s), bs in ids.items() if t == r.t1 for b in bs]
        all_b = [b for (t, s), bs in ids.items() if t == r.t2 for b in bs]
        fwd, back = conn_count(client, all_a, all_b), conn_count(client, all_b, all_a)
        n_sides = len(sides) if fwd >= 0 and back >= 0 else 0
        det = [f"合計 {max(fwd, 0)}/{max(back, 0)}"]
        fwd, back = max(fwd, 0), max(back, 0)
        src_fwd, src_rev = int(r.synapses), int(rev.get((r.t2, r.t1), 0))
        rows.append({
            "t1": r.t1, "t2": r.t2, "dataset": dataset,
            "contact_um_src": r.contact_total_um, "expected_src": r.expected,
            "syn_src_fwd": src_fwd, "syn_src_rev": src_rev,
            "n_sides": n_sides, "syn_cross_fwd": fwd if n_sides else -1,
            "syn_cross_rev": back if n_sides else -1,
            "cross_total": (fwd + back) if n_sides else -1,
            "detail": " ".join(det) if det else "（タイプ名なし）",
            "total_syn": src_fwd + src_rev + (fwd + back if n_sides else 0),
            "comp_type": r.comp_type, "comp_rate": r.comp_rate,
            "availability": getattr(r, "availability", np.nan),
        })
        print(f"{r.t1:>14s} → {r.t2:<14s} 元 {src_fwd}/{src_rev} | {dataset} {rows[-1]['detail']}"
              + ("  ★" if rows[-1]["total_syn"] == 0 and n_sides else ""))
    return pd.DataFrame(rows)


def stage_contact(cross: pd.DataFrame, dataset: str, th: Thresholds, top: int, out_dir: Path,
                  max_total: int, pairs: pd.DataFrame | None = None) -> pd.DataFrame:
    """そのデータセットで、代表ニューロンの接触長（アプローチ G）を計算する。"""
    from approach_g_output_sites.run import output_contacts
    from common.loaders import load_neuprint, neuprint_client

    client = neuprint_client(dataset)
    if pairs is not None:  # shortlist.csv など、外から組を指定された場合
        want = set(zip(pairs.t1, pairs.t2))
        sel = cross[[(a, b) in want for a, b in zip(cross.t1, cross.t2)]]
        sel = sel[cross.n_sides > 0] if "n_sides" in sel else sel
    else:
        sel = cross[(cross.total_syn <= max_total) & (cross.n_sides > 0)]
    sel = sel.sort_values("contact_um_src", ascending=False).head(top)
    rows = []
    for r in sel.itertuples(index=False):
        ids, _ = side_ids(client, [r.t1, r.t2])
        for s in sorted({k[1] for k in ids}):
            a, b = ids.get((r.t1, s), []), ids.get((r.t2, s), [])
            if not a or not b:
                continue
            ds2 = load_neuprint(a[:4] + b[:4], dataset=dataset)
            ds2.save(out_dir / f"{tag_of(dataset)}_{r.t1}_{r.t2}_{s}.pkl".replace("/", "_"))
            df = output_contacts(ds2, th)
            ta = {i for i, n in ds2.neurons.items() if str(n.meta.get("type")) == r.t1}
            tb = {i for i, n in ds2.neurons.items() if str(n.meta.get("type")) == r.t2}
            d = df[df.pre.isin(ta) & df.post.isin(tb)]
            rows.append({"t1": r.t1, "t2": r.t2, "dataset": dataset, "side": s,
                         "n1": len(ta), "n2": len(tb), "n_pairs_contact": len(d),
                         "contact_um": float(d.contact_um.sum())})
            print(f"{r.t1} → {r.t2} [{dataset} {s}] 接触 {rows[-1]['contact_um']:.1f} µm（{len(d)} ペア）")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="元データ（hemibrain）の pickle")
    ap.add_argument("--datasets", nargs="+", default=DEFAULT_DATASETS, help="確かめる先の neuPrint データセット")
    ap.add_argument("--stage", default="conn", choices=["conn", "contact"])
    ap.add_argument("--top", type=int, default=60)
    ap.add_argument("--min-expected", type=float, default=8.0, help="元データでの期待シナプス数の下限")
    ap.add_argument("--min-syn-annotated", type=int, default=100,
                    help="そのデータセットで見える結合の重みの下限（上流の出力・下流の入力それぞれ）。"
                         "これ未満なら『結合があっても見えない』として照合不能にする")
    ap.add_argument("--max-total", type=int, default=2,
                    help="contact: 両方向・全データセットの合計シナプス数の上限")
    ap.add_argument("--pairs", default=None,
                    help="contact: 調べる組を書いた CSV（列 t1, t2）。shortlist.csv を渡せる")
    ap.add_argument("--pairs-as-source", action="store_true",
                    help="--pairs の CSV をそのまま候補として使う（複数領域にまたがる組を調べるとき）")
    ap.add_argument("--tag", default=None, help="出力先 results/<tag>/")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()

    th = Thresholds().update(parse_overrides(args.set))
    res = ROOT / "results" / (args.tag or Path(args.data).stem)
    res.mkdir(parents=True, exist_ok=True)

    if args.stage == "conn":
        if args.pairs_as_source:
            cand = pd.read_csv(args.pairs)
        else:
            cand = pd.read_csv(ROOT / "results" / Path(args.data).stem / "candidates_ranked.csv")
        cand = cand[cand.synapses == 0]
        if args.pairs and not args.pairs_as_source:  # 調べる組を外から指定（他のデータセットで 0 だった組だけ、など）
            p = pd.read_csv(args.pairs)
            want = set(zip(p.t1, p.t2))
            cand = cand[[(a, b) in want for a, b in zip(cand.t1, cand.t2)]]
        else:
            cand = cand[cand.expected >= args.min_expected]
        cand = cand.sort_values("expected", ascending=False)
        rev = reverse_in_source(Dataset.load(args.data))
        for dsname in args.datasets:
            print(f"\n----- {dsname}")
            out = stage_conn(cand, rev, dsname, args.top, args.min_syn_annotated)
            f = res / f"cross_{tag_of(dsname)}.csv"
            out.to_csv(f, index=False, encoding="utf-8")
            ok = out[(out.total_syn == 0) & (out.n_sides > 0)]
            print(f"{len(out)} 組を確認（照合できたのは {int((out.n_sides > 0).sum())} 組）、"
                  f"合計 0 の組 {len(ok)} 組 -> {f}")
    else:
        d = ROOT / "data" / "cross"
        d.mkdir(parents=True, exist_ok=True)
        for dsname in args.datasets:
            f = res / f"cross_{tag_of(dsname)}.csv"
            if not f.exists():
                print(f"skip {dsname}: {f} がありません（先に stage conn）")
                continue
            print(f"\n----- {dsname}")
            pairs = pd.read_csv(args.pairs) if args.pairs else None
            out = stage_contact(pd.read_csv(f), dsname, th, args.top, d, args.max_total, pairs)
            g = res / f"cross_{tag_of(dsname)}_contact.csv"
            save_merged(out, g, ["t1", "t2", "side"])
            print(f"-> {g}")


if __name__ == "__main__":
    main()
