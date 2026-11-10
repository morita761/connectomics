#!/usr/bin/env bash
# よく使う処理をまとめたスクリプト（Makefile の代わり。Git Bash / Linux / macOS 用）
#
#   ./tasks.sh server                        結果サイトを作ってローカルサーバーで開く（http://127.0.0.1:8000/）
#   ./tasks.sh site                          結果サイトを作るだけ
#   ./tasks.sh benchmark                     合成データで全アプローチを比較
#   ./tasks.sh analyze data/hemibrain_mb.pkl
#                                            全アプローチ → 一致ペアの集計 → タイプの組の集計 → サイト作成
#
# 環境変数で設定を変えられる:  PORT=8080 N3D=50 ./tasks.sh server
set -euo pipefail
cd "$(dirname "$0")"
if [ -x .venv/Scripts/python.exe ]; then PY=.venv/Scripts/python.exe; else PY=.venv/bin/python; fi
PORT="${PORT:-8000}"
N3D="${N3D:-30}"

run() { echo "> python $*"; "$PY" "$@"; }

case "${1:-help}" in
  server)    run serve.py --port "$PORT" --n3d "$N3D" ;;
  site)      run build_site.py --n3d "$N3D" ;;
  benchmark) run run_benchmark.py ;;
  analyze)
    data="${2:-data/hemibrain_mb.pkl}"
    stem="$(basename "$data" .pkl)"
    run run_all.py --data "$data" --compartments navis
    run compare_methods.py --results "results/$stem" --data "$data" --top-pct 1 --min-methods 4
    run type_pairs.py --data "$data"
    run build_site.py --n3d "$N3D"
    ;;
  *) sed -n '3,10p' "$0" | sed 's/^# \{0,1\}//' ;;
esac
