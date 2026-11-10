"""GitHub Pages で公開する分だけを docs/ に集める。

  python build_docs.py && python publish_site.py
  python publish_site.py --with-datasets   # データセット別ページも入れる（3D は巨大なので除く）

なぜ全部を公開しないのか:
  data/ が 5.4 GB、results/ の 3D（results/<データ名>/3d/）が合計 3.5 GB ある。
  GitHub は 1 ファイル 100 MB・リポジトリ 1 GB が目安なので、そのままでは載らない。
  報告書・最終候補・候補ペアの 3D だけなら約 100 MB で収まる。

集めるもの:
  results/docs/        報告書（Markdown を HTML にしたもの）
  results/shortlist.html  最終候補
  results/zero_3d/, results/cross_3d/  候補ペアの 3D
  results/*.md, *.csv のうち報告書が参照するもの
  docs/index.html      入口（ここで作る）
  docs/reports/        報告書（results/docs/ をコピーしたもの）

GitHub Pages の「Deploy from a branch」はフォルダを / (root) か /docs しか選べないので、
出力先は docs/ に固定している。
"""
from __future__ import annotations

import argparse
import html
import shutil
from pathlib import Path

from common.runner import ROOT

RESULTS = ROOT / "results"
SITE = ROOT / "docs"   # GitHub Pages が選べるフォルダは / か /docs だけ
LIMIT_MB = 100  # GitHub の 1 ファイル上限


def copy(src: Path, dst: Path) -> int:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
        return sum(f.stat().st_size for f in dst.rglob("*") if f.is_file())
    shutil.copy2(src, dst)
    return dst.stat().st_size


INDEX = """<!doctype html>
<html lang="ja"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>近接するがシナプスを作らないニューロンの探索</title>
<style>
:root{{--bg:#fff;--fg:#1a1a1a;--sub:#666;--line:#e3e3e3;--acc:#0b5fa5;--code:#f6f8fa}}
@media (prefers-color-scheme:dark){{:root{{--bg:#15171a;--fg:#e8e8e8;--sub:#9aa0a6;--line:#2c2f33;--acc:#6cb4ff;--code:#1e2227}}}}
body{{margin:0;background:var(--bg);color:var(--fg);font-family:system-ui,"Segoe UI","Hiragino Sans","Noto Sans JP",sans-serif;line-height:1.7}}
.wrap{{max-width:900px;margin:0 auto;padding:28px 16px 80px}}
h1{{font-size:1.7rem;margin-bottom:.2em}}
.lead{{color:var(--sub);margin-top:0}}
.cards{{display:grid;gap:12px;grid-template-columns:repeat(auto-fill,minmax(260px,1fr));margin-top:24px}}
.card{{display:block;border:1px solid var(--line);border-radius:8px;padding:16px;text-decoration:none;color:var(--fg)}}
.card:hover{{border-color:var(--acc)}}
.card b{{color:var(--acc);font-size:1.05rem}}
.card span{{display:block;color:var(--sub);font-size:.88rem;margin-top:6px}}
.hero{{grid-column:1/-1;background:var(--code)}}
ul{{color:var(--sub);font-size:.92rem}}
code{{background:var(--code);padding:.15em .4em;border-radius:4px}}
</style></head><body><div class="wrap">
<h1>近接するがシナプスを作らないニューロンの探索</h1>
<p class="lead">ショウジョウバエのコネクトーム 3 個体（hemibrain・FlyWire・male-CNS）を突き合わせて、
<b>膜は接しているのにシナプスが 1 つも無いニューロンのペア</b>を探す。
人工接着分子で結合を強制的に作らせる実験の標的を選ぶのが目的。
条件は 5 つ — 両方向でシナプス 0 ／ 複数個体で 0 ／ EM で膜が接している ／
GAL4 と LexA の 2 系統が使える ／ 上流がその場所で出力している。</p>
<div class="cards">
{cards}
</div>
<h2>公開していないもの</h2>
<ul>
<li>骨格・メッシュのデータ（<code>data/</code>、約 5.4 GB）</li>
<li>データセット別ページの全ペア 3D（<code>results/&lt;データ名&gt;/3d/</code>、約 3.5 GB）</li>
</ul>
<p class="lead">どちらも <code>fetch_data.py</code> と <code>run_all.py</code> で作り直せる。</p>
</div></body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-datasets", action="store_true",
                    help="データセット別ページも入れる（3D は除くのでリンクは切れる）")
    ap.add_argument("--clean", action="store_true", help="docs/ を作り直す")
    args = ap.parse_args()

    if args.clean and SITE.exists():
        shutil.rmtree(SITE)
    SITE.mkdir(exist_ok=True)

    docs = RESULTS / "docs"
    if not docs.exists():
        print("results/docs/ がない。先に python build_docs.py を実行してください。")
        return

    total = 0
    total += copy(docs, SITE / "reports")
    for d in ("zero_3d", "cross_3d"):
        if (RESULTS / d).exists():
            total += copy(RESULTS / d, SITE / d)
    for f in ("shortlist.html", "shortlist_final.csv", "mb_final19.csv",
              "gmr_lexa_match.csv", "hemidriver_lexa.csv", "em_gap/summary.csv"):
        src = RESULTS / f
        if src.exists():
            total += copy(src, SITE / f)
    # 最終候補ページから報告書へのリンクは results/docs/ を指しているので、
    # 公開先の docs/reports/ に読み替える
    sl = SITE / "shortlist.html"
    if sl.exists():
        t = sl.read_text(encoding="utf-8")
        # 行データは JSON として埋め込まれていて引用符が \" になっているので、両方の形を直す
        t = t.replace('href="docs/', 'href="reports/').replace('href=\\"docs/', 'href=\\"reports/')
        sl.write_text(t, encoding="utf-8")
    if (ROOT / "data" / "lexa_lines_nonmb.csv").exists():
        total += copy(ROOT / "data" / "lexa_lines_nonmb.csv", SITE / "lexa_lines_nonmb.csv")
    if args.with_datasets:
        for d in sorted(p for p in RESULTS.iterdir() if p.is_dir() and (p / "index.html").exists()):
            total += copy(d / "index.html", SITE / d.name / "index.html")

    n_docs = len(list((SITE / "reports").glob("*.html"))) - 1
    n3d = len(list((SITE / "zero_3d").glob("*.html"))) + len(list((SITE / "cross_3d").glob("*.html")))
    cards = [
        f'<a class="card hero" href="reports/index.html"><b>◇ 報告書（{n_docs} 件）</b>'
        f'<span>候補の検討・ドライバー・EM・文献のまとめ</span></a>',
        '<a class="card" href="reports/candidate_MBON04_PPL105.html"><b>本命: MBON04 ↔ PPL105</b>'
        '<span>キノコ体案の第一候補。両方向が候補になる</span></a>',
        '<a class="card" href="reports/candidate_MBON04_FB4N.html"><b>本命: MBON04 → FB4N</b>'
        '<span>全脳探索の第一候補</span></a>',
        '<a class="card" href="shortlist.html"><b>★ 最終候補の一覧</b>'
        '<span>3 個体で確かめた組</span></a>',
        f'<a class="card" href="reports/em_links_MBON04_PPL105.html"><b>EM と 3D のリンク</b>'
        f'<span>neuroglancer で膜の接触を見る ／ 3D {n3d} 件</span></a>',
        '<a class="card" href="reports/lexa_catalog.html"><b>LexA ドライバーの棚卸し</b>'
        '<span>キノコ体以外・神経とグリア</span></a>',
    ]
    (SITE / ".nojekyll").write_text("", encoding="utf-8")  # Jekyll の処理を止める
    (SITE / "index.html").write_text(INDEX.format(cards="\n".join(cards)), encoding="utf-8")

    big = [f for f in SITE.rglob("*") if f.is_file() and f.stat().st_size > LIMIT_MB * 1e6]
    print(f"docs/ に {sum(1 for f in SITE.rglob('*') if f.is_file())} ファイル / "
          f"{total / 1e6:.0f} MB")
    if big:
        print(f"警告: {LIMIT_MB} MB を超えるファイルがある（GitHub に push できない）:")
        for f in big:
            print(f"  {f.stat().st_size / 1e6:.0f} MB  {f.relative_to(SITE)}")
    print(f"-> {SITE / 'index.html'}")
    print("""
GitHub Pages に載せる手順:
  git init                        # まだ Git リポジトリでない場合
  git add .gitignore docs         # data/ と巨大な 3D は .gitignore で除外済み
  git commit -m "公開用サイト"
  git branch -M main
  git remote add origin https://github.com/<ユーザー名>/<リポジトリ名>.git
  git push -u origin main
そのあと GitHub の Settings → Pages で
  Source: Deploy from a branch / Branch: main / フォルダ: /docs
を選ぶと https://<ユーザー名>.github.io/<リポジトリ名>/ で公開される。
（選べるのは / か /docs だけ。別名のフォルダを公開したいときは GitHub Actions の
  ワークフローを使うか、gh-pages ブランチの root に置く）""")


if __name__ == "__main__":
    main()
