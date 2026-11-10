#!/usr/bin/env bash
# 1 つの領域について、3 個体での絞り込みをまとめて行う。
#   ./run_screen.sh smp            # 既に run_regions.sh で取得・集計した領域に対して
# 順番: 比較相手を付けて順位付け → FlyWire で 0 の組だけ残す → male-CNS で確認 → 接触領域を調べる
set -euo pipefail
cd "$(dirname "$0")"
if [ -x .venv/Scripts/python.exe ]; then PY=.venv/Scripts/python.exe; else PY=.venv/bin/python; fi
export PYTHONIOENCODING=utf-8

k=${1:?領域名を指定してください}
data="data/hemibrain_$k.pkl"
res="results/hemibrain_$k"

echo "########## $k 順位付け $(date +%H:%M)"
"$PY" rank_candidates.py --data "$data" --top 0 | tail -1

echo "########## $k FlyWire $(date +%H:%M)"
"$PY" flywire_check.py --data "$data" --min-expected 3 --top 2000 | tail -1
"$PY" - "$res" "results/zero_$k.csv" <<'PY'
import sys
import pandas as pd
d = pd.read_csv(f"{sys.argv[1]}/cross_flywire_783.csv")
z = d[(d.n_sides > 0) & (d.total_syn == 0)]
z.to_csv(sys.argv[2], index=False)
print("hemibrain + FlyWire でゼロ:", len(z))
PY

echo "########## $k male-CNS $(date +%H:%M)"
"$PY" validate_cross.py --data "$data" --top 300 --datasets male-cns:v1.0 \
  --pairs "results/zero_$k.csv" 2>&1 | grep -avE "it/s|s/it" | tail -2
"$PY" - "$res" "results/zero3_$k.csv" <<'PY'
import sys
import pandas as pd
d = pd.read_csv(f"{sys.argv[1]}/cross_male_cns_v1_0.csv")
z = d[(d.n_sides > 0) & (d.total_syn == 0)]
z.to_csv(sys.argv[2], index=False)
print("3 個体すべてでゼロ:", len(z))
PY

echo "########## $k 接触領域 $(date +%H:%M)"
"$PY" inspect_pairs.py --data "$data" --any-region --pairs "results/zero3_$k.csv" | tail -1
echo "########## $k done $(date +%H:%M)"
