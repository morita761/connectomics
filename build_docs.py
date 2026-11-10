"""results/ の Markdown 報告書を HTML に変換し、サイトから読めるようにする。

  python build_docs.py                 # results/*.md -> results/docs/*.html と docs/index.html
  python build_docs.py --open          # 変換後に一覧をブラウザで開く

Markdown はブラウザがそのまま表示してくれないので（生のテキストになる）、
ビルド時に HTML に変換しておく。表・コードブロックを使っているので
markdown の tables / fenced_code 拡張を有効にする。

リンクの書き換え:
  - 報告書の中の `results/...` は docs/ から見ると 1 つ上なので `../...` にする
  - `data/...` は `../../data/...`
  - 変換対象の .md へのリンクは .html に差し替える（報告書どうしを行き来できる）
"""
from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

from common.runner import ROOT

RESULTS = ROOT / "results"
DOCS = RESULTS / "docs"

CSS = """
:root{--bg:#fff;--fg:#1a1a1a;--sub:#666;--line:#e3e3e3;--acc:#0b5fa5;--code:#f6f8fa}
@media (prefers-color-scheme:dark){:root{--bg:#15171a;--fg:#e8e8e8;--sub:#9aa0a6;--line:#2c2f33;--acc:#6cb4ff;--code:#1e2227}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font-family:system-ui,"Segoe UI","Hiragino Sans","Noto Sans JP",sans-serif;line-height:1.75}
.wrap{max-width:900px;margin:0 auto;padding:24px 16px 80px}
nav{font-size:.9rem;color:var(--sub);margin-bottom:18px}
nav a{color:var(--acc);text-decoration:none}nav a:hover{text-decoration:underline}
h1{font-size:1.6rem;line-height:1.35;border-bottom:2px solid var(--line);padding-bottom:.4em;margin-top:1.2em}
h2{font-size:1.25rem;margin-top:2em;border-bottom:1px solid var(--line);padding-bottom:.3em}
h3{font-size:1.05rem;margin-top:1.6em}
a{color:var(--acc)}
table{border-collapse:collapse;width:100%;margin:1.2em 0;font-size:.92rem;display:block;overflow-x:auto}
th,td{border:1px solid var(--line);padding:6px 10px;text-align:left;vertical-align:top}
th{background:var(--code);font-weight:600;white-space:nowrap}
code{background:var(--code);padding:.15em .4em;border-radius:4px;font-size:.9em}
pre{background:var(--code);padding:12px 14px;border-radius:6px;overflow-x:auto}
pre code{background:none;padding:0}
blockquote{margin:1.2em 0;padding:.6em 1em;border-left:4px solid var(--acc);background:var(--code);color:var(--fg)}
ul,ol{padding-left:1.4em}
hr{border:0;border-top:1px solid var(--line);margin:2em 0}
.cards{display:grid;gap:12px;grid-template-columns:repeat(auto-fill,minmax(260px,1fr));margin-top:20px}
.card{display:block;border:1px solid var(--line);border-radius:8px;padding:14px 16px;text-decoration:none;color:var(--fg);background:var(--bg)}
.card:hover{border-color:var(--acc)}
.card b{color:var(--acc)}
.card span{display:block;color:var(--sub);font-size:.85rem;margin-top:4px}
"""

PAGE = """<!doctype html>
<html lang="ja"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{css}</style></head>
<body><div class="wrap">
<nav><a href="index.html">← 報告書の一覧</a> ／ <a href="../index.html">結果一覧</a></nav>
{body}
</div></body></html>
"""


def rewrite_links(md: str, md_names: set[str]) -> str:
    """Markdown のリンク先を docs/ から辿れる形に直す。"""
    def fix(m):
        label, target = m.group(1), m.group(2)
        if re.match(r"^(https?:|mailto:|#)", target):
            return m.group(0)
        t = target
        if t.startswith("results/"):
            t = "../" + t[len("results/"):]
        elif t.startswith("data/"):
            t = "../../" + t
        # 同じ docs/ に変換される .md は .html へ
        stem = Path(t).name
        if stem.endswith(".md") and stem[:-3] in md_names:
            t = stem[:-3] + ".html"
        return f"[{label}]({t})"
    return re.sub(r"\[([^\]]*)\]\(([^)\s]+)\)", fix, md)


def convert(files: list[Path]) -> list[tuple[Path, str, str]]:
    import markdown as md_lib
    DOCS.mkdir(parents=True, exist_ok=True)
    names = {f.stem for f in files}
    made = []
    for f in files:
        text = rewrite_links(f.read_text(encoding="utf-8"), names)
        body = md_lib.markdown(text, extensions=["tables", "fenced_code", "toc", "sane_lists"])
        m = re.search(r"<h1[^>]*>(.*?)</h1>", body, re.S)
        title = re.sub(r"<[^>]+>", "", m.group(1)).strip() if m else f.stem
        # 概要として最初の段落を取る（一覧カードの説明に使う）
        p = re.search(r"<p>(.*?)</p>", body, re.S)
        lead = re.sub(r"<[^>]+>", "", p.group(1)).strip()[:110] if p else ""
        out = DOCS / f"{f.stem}.html"
        out.write_text(PAGE.format(title=html.escape(title), css=CSS, body=body), encoding="utf-8")
        made.append((out, title, lead))
    return made


def index(made: list[tuple[Path, str, str]]) -> Path:
    cards = "\n".join(
        f'<a class="card" href="{p.name}"><b>{html.escape(t)}</b>'
        f'<span>{html.escape(lead)}</span></a>'
        for p, t, lead in sorted(made, key=lambda x: x[1]))
    body = ("<h1>報告書の一覧</h1>\n<p>results/ の Markdown 報告書を HTML にしたもの。"
            f"全 {len(made)} 件。</p>\n<div class=\"cards\">\n{cards}\n</div>")
    f = DOCS / "index.html"
    f.write_text(PAGE.format(title="報告書の一覧", css=CSS, body=body)
                 .replace('<nav><a href="index.html">← 報告書の一覧</a> ／ ', "<nav>"), encoding="utf-8")
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--open", action="store_true", help="変換後に一覧をブラウザで開く")
    args = ap.parse_args()

    files = sorted(f for f in RESULTS.glob("*.md") if f.name != "index.md")
    if not files:
        print("results/*.md が無い")
        return
    made = convert(files)
    for p, t, _ in made:
        print(f"  {t} -> {p.relative_to(ROOT)}")
    f = index(made)
    print(f"\n{len(made)} 件 -> {f.relative_to(ROOT)}")
    if args.open:
        import webbrowser
        webbrowser.open(f.as_uri())


if __name__ == "__main__":
    main()
