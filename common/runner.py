"""各アプローチの run.py が共有する CLI 部分。"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path
from typing import Callable

import pandas as pd

from .config import Thresholds, parse_overrides
from .model import Dataset

ROOT = Path(__file__).resolve().parents[1]


def list_approaches() -> list[str]:
    return sorted(p.parent.name for p in ROOT.glob("approach_*/run.py"))


def load_approach(name: str):
    """approach_xxx/run.py をモジュールとして読み込む（フォルダごとに独立させるため importlib で）。"""
    mod_name = f"_approach_{name}"
    if mod_name in sys.modules:
        return sys.modules[mod_name]
    spec = importlib.util.spec_from_file_location(mod_name, ROOT / name / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


def cli(name: str, find_candidates: Callable[..., pd.DataFrame], extra_args=None) -> None:
    ap = argparse.ArgumentParser(description=name)
    ap.add_argument("--data", required=True, help="Dataset の pickle（fetch_data.py / run_benchmark.py で作成）")
    ap.add_argument("--out", default=None, help="結果 CSV（既定: results/<approach>.csv）")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE", help="閾値の上書き")
    ap.add_argument("--top", type=int, default=20)
    if extra_args:
        extra_args(ap)
    args = ap.parse_args()

    th = Thresholds().update(parse_overrides(args.set))
    ds = Dataset.load(args.data)
    kw = {k: v for k, v in vars(args).items() if k not in ("data", "out", "set", "top")}
    t0 = time.perf_counter()
    df = find_candidates(ds, th, **kw)
    dt = time.perf_counter() - t0

    out = Path(args.out or Path(__file__).resolve().parents[1] / "results" / f"{name}.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    flagged = int(df["flagged"].sum()) if "flagged" in df else 0
    print(f"[{name}] {len(df)} proximal pairs, {flagged} flagged, {dt:.1f}s -> {out}")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(df.head(args.top).to_string(index=False))


def save_merged(df, path, keys):
    """前回までの結果を消さずに保存する。

    組を指定して少しずつ測るとき、毎回上書きすると前に測った組が消えてしまう。
    keys が一致する行だけ入れ替え、それ以外は残す。
    """
    import pandas as pd
    path = Path(path)
    if path.exists() and len(df):
        old = pd.read_csv(path)
        if all(k in old.columns for k in keys):
            keep = old[~old.set_index(keys).index.isin(df.set_index(keys).index)]
            df = pd.concat([keep, df], ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8")
    return df
