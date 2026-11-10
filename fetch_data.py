"""実データを取得して Dataset (pickle) として保存する。

  # hemibrain: 細胞タイプ or bodyId を指定（トークンは環境変数 NEUPRINT_APPLICATION_CREDENTIALS）
  python fetch_data.py neuprint --types "KCg-m" "MBON01" --out data/hemibrain_mb.pkl
  python fetch_data.py neuprint --ids 5813105172 1009378435 --meshes

  # FlyWire (公開版 783): root ID を指定（fafbseg の CAVE トークン設定が必要）
  python fetch_data.py flywire --ids 720575940621039145 720575940613583001 --out data/fw.pkl

保存後は各アプローチを
  python approach_a_skeleton_kdtree/run.py --data data/hemibrain_mb.pkl
のように実行するか、run_all.py でまとめて実行する。

※ 対象ニューロン数の目安: スケルトンなら数百〜数千、メッシュ (B) は数十〜数百。
   全ペア比較なので、まず ROI や細胞タイプで絞ってから使う。
"""
from __future__ import annotations

import argparse

from common.loaders import load_flywire, load_neuprint


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="source", required=True)

    np_ = sub.add_parser("neuprint")
    np_.add_argument("--ids", nargs="*", type=int, default=[])
    np_.add_argument("--types", nargs="*", default=[], help="細胞タイプ名（正規表現可）")
    np_.add_argument("--dataset", default="hemibrain:v1.2.1")
    np_.add_argument("--group", nargs="*", default=[], metavar="NAME=REGEX",
                     help="細胞タイプの正規表現でグループを定義（例: visual='LC.*' olfactory='.*_l?PN'）。"
                          "メタデータ group に保存され、type_pairs.py でグループ間の集計に使う")
    np_.add_argument("--rois", nargs="*", default=[], help="この ROI のいずれかにシナプスをもつニューロンに限る（例: LH(R) PLP(R)）")
    np_.add_argument("--min-roi-synapses", type=int, default=20,
                     help="--rois 内の入力+出力シナプスがこれ以上のニューロンに限る")
    np_.add_argument("--max-per-type", type=int, default=0, help="1 タイプあたりの最大ニューロン数（0 = 制限なし。多いタイプはランダムに間引く）")
    np_.add_argument("--dry-run", action="store_true", help="取得せず、該当するニューロン数だけ表示")
    np_.add_argument("--meshes", action="store_true")
    np_.add_argument("--out", default="data/neuprint.pkl")

    meta = sub.add_parser("neuprint-meta", help="既存の pickle に細胞名 (instance) などを追加する（スケルトンは取り直さない）")
    meta.add_argument("--data", required=True)
    meta.add_argument("--dataset", default="hemibrain:v1.2.1")

    fw = sub.add_parser("flywire")
    fw.add_argument("--ids", nargs="+", type=int, required=True)
    fw.add_argument("--materialization", type=int, default=783)
    fw.add_argument("--meshes", action="store_true")
    fw.add_argument("--types-csv", default=None,
                    help="細胞タイプの表（列 root_id と primary_type、または root_id と type）。アプローチ E 用")
    fw.add_argument("--out", default="data/flywire.pkl")

    a = ap.parse_args()
    if a.source == "neuprint-meta":
        from neuprint import NeuronCriteria as NC, fetch_neurons
        from common.loaders import neuprint_client
        from common.model import Dataset
        ds = Dataset.load(a.data)
        df, _ = fetch_neurons(NC(bodyId=list(ds.neurons)), client=neuprint_client(a.dataset))
        cols = [c for c in ("type", "instance", "status", "cropped") if c in df]
        for r in df[["bodyId"] + cols].itertuples(index=False):
            n = ds.neurons.get(int(r.bodyId))
            if n is not None:
                n.meta.update({c: getattr(r, c) for c in cols})
        ds.save(a.data)
        print(f"updated {len(df)} neurons ({', '.join(cols)}) -> {a.data}")
        return
    if a.source == "neuprint":
        import re
        ids = list(a.ids)
        groups = dict(g.split("=", 1) for g in a.group)
        patterns = list(a.types) + list(groups.values())
        if patterns:
            from neuprint import NeuronCriteria as NC, fetch_neurons
            from common.loaders import neuprint_client
            c = neuprint_client(a.dataset)
            df, roi_df = fetch_neurons(NC(type="|".join(f"({p})" for p in patterns), regex=True,
                                          rois=a.rois or None, roi_req="any"), client=c)
            if a.rois:
                r = roi_df[roi_df.roi.isin(a.rois)].assign(n=lambda d: d.pre + d.post).groupby("bodyId").n.sum()
                df = df[df.bodyId.map(r).fillna(0) >= a.min_roi_synapses]
            if a.max_per_type:
                df = df.sample(frac=1, random_state=0).groupby("type").head(a.max_per_type)
            for name, pat in groups.items():
                hit = df.type.fillna("").str.fullmatch(pat)
                print(f"group {name}: {hit.sum()} neurons, {df.type[hit].nunique()} types")
            print(f"total {len(df)} neurons, {df.type.nunique()} types")
            print(df.type.value_counts().head(40).to_string())
            if a.dry_run:
                return
            ids += df.bodyId.tolist()
        if not ids:
            ap.error("--ids / --types / --group のどれかを指定してください")
        ds = load_neuprint(sorted(set(ids)), dataset=a.dataset, with_meshes=a.meshes)
        for n in ds.neurons.values():
            t = n.meta.get("type") or ""
            n.meta["group"] = next((g for g, pat in groups.items() if re.fullmatch(pat, t)), None)
    else:
        types = None
        if a.types_csv:
            import pandas as pd
            t = pd.read_csv(a.types_csv)
            col = "primary_type" if "primary_type" in t else "type"
            types = dict(zip(t.root_id.astype("int64"), t[col]))
        ds = load_flywire(a.ids, materialization=a.materialization, with_meshes=a.meshes, types=types)

    ds.save(a.out)
    print(f"saved {len(ds.neurons)} neurons, {len(ds.edges)} edges -> {a.out}")


if __name__ == "__main__":
    main()
