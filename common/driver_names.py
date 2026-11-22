"""ドライバー表の「区画名」を hemibrain の型名に読み替える。

なぜ要るか:
  ドライバーの出典が 2 系統あって、細胞の呼び方が違う。

    Aso 2014 のキノコ体コレクション   `MBON-α1`  `PPL1-α'2α2`  `PAM-α1`   （区画名）
    Shuai 2024 / Meissner 2024        `MBON07`   `PPL105`      `PAM11`    （hemibrain の型名）

  照合は hemibrain の型名で走っているので、**区画名で書かれた系統は一度も候補に出てこない**。
  これが `MB058B`（PPL105 を単独で狙える唯一の系統）や
  `MB310C`（MBON07 を単独で狙える、しかも BDSC にある）を見落とした原因だった
  （`results/stocks_bdsc.md`）。

どう対応をつけるか:
  hemibrain の instance 名が両方を持っている。

    MBON07(a1)_R          -> MBON-α1
    PPL105(a'2a2)(PDL05)_L -> PPL1-α'2α2
    MBON04(B'2mp_bilateral)_R -> MBON-β'2mp_bilateral

  instance の最初の括弧が区画名で、ギリシャ文字が ASCII に置き換えてある
  （α→a、β→B、γ→y）。ただし `B'2mp_bilateral` の `a` のように、
  ギリシャ文字ではない `a` も混ざる。置換規則を決め打ちすると壊れるので、
  **両方の表記を同じ形に畳んでから突き合わせる**（canon）。
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

from common.model import Dataset
from common.runner import ROOT

CACHE = ROOT / "data" / "compartment_to_type.csv"

# ギリシャ文字と、hemibrain instance の ASCII 表記を同じ文字に畳む。
# アポストロフィは出典によって ′ ’ ´ が混在するので 1 種類に揃える。
# **落としてはいけない。** α'1 と α1、α'3 と α3 が別の細胞なので、
# 落とすと MBON15 と MBON07、PPL104 と PPL106 が衝突する。
_FOLD = str.maketrans({
    "α": "a", "β": "b", "γ": "g",
    "B": "b", "y": "g",
    "′": "'", "’": "'", "´": "'", "‘": "'",
    " ": "",
})


def canon(s: str) -> str:
    """区画名を、出典による表記ゆれを無視した形にする。"""
    return str(s).translate(_FOLD).lower()


def _family(t: str) -> str | None:
    if t.startswith("MBON"):
        return "MBON-"
    if t.startswith("PPL"):
        return "PPL1-"
    if t.startswith("PAM"):
        return "PAM-"
    return None


def compartment_of(type_name: str, instance: str) -> str | None:
    """hemibrain の instance 名から区画名を作る。MBON07(a1)_R -> MBON-a1"""
    m = re.match(r"[^(]+\(([^)]*)\)", str(instance or ""))
    fam = _family(type_name)
    return fam + m.group(1) if (m and fam) else None


def build(verbose: bool = False) -> dict[str, str]:
    """data/hemibrain_*.pkl の instance 名から {canon(区画名): 型名} を作る。

    同じ畳み方になる区画名が 2 つの型に付いたら、**どちらも捨てる**。
    取り違えるくらいなら拾わないほうがよい。
    """
    seen: dict[str, set[str]] = {}
    for f in sorted(glob.glob(str(ROOT / "data" / "hemibrain_*.pkl"))):
        try:
            ds = Dataset.load(f)
        except Exception:
            continue
        for n in ds.neurons.values():
            t, inst = n.meta.get("type"), n.meta.get("instance")
            if not t or not inst:
                continue
            c = compartment_of(t, inst)
            if c:
                seen.setdefault(canon(c), set()).add(t)
    out, dropped = {}, []
    for k, v in seen.items():
        (out.__setitem__(k, next(iter(v))) if len(v) == 1 else dropped.append((k, sorted(v))))
    if verbose and dropped:
        print(f"衝突したので捨てた区画名 {len(dropped)} 件:")
        for k, v in dropped:
            print(f"  {k} -> {v}")
    return out


def load(rebuild: bool = False) -> dict[str, str]:
    """対応表。無ければ作って data/compartment_to_type.csv に保存する。"""
    if CACHE.exists() and not rebuild:
        rows = [l.split(",", 1) for l in CACHE.read_text(encoding="utf-8").splitlines()[1:] if "," in l]
        return {a: b for a, b in rows}
    m = build()
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text("canon,type\n" + "".join(f"{k},{v}\n" for k, v in sorted(m.items())),
                     encoding="utf-8")
    return m


def to_type(name: str, table: dict[str, str] | None = None) -> str:
    """区画名なら hemibrain の型名に、そうでなければそのまま返す。"""
    name = str(name).strip()
    if not name.startswith(("MBON-", "PPL1-", "PAM-")):
        return name
    return (table if table is not None else load()).get(canon(name), name)


if __name__ == "__main__":
    m = build(verbose=True)
    load(rebuild=True)
    print(f"\n{len(m)} 件 -> {CACHE}")
    for k in sorted(m)[:12]:
        print(f"  {k:24s} -> {m[k]}")
