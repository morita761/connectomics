#!/usr/bin/env bash
# キノコ体以外の領域で「近いのに結合しない」タイプの組を探す（取得 → G → タイプの組 → ドライバー照合）。
#   ./run_regions.sh            # 3 領域すべて
#   ./run_regions.sh lh cx      # 一部だけ
# 領域とグループ（上流→下流の階層）は下の case で定義。データがあれば取得は飛ばす。
set -euo pipefail
cd "$(dirname "$0")"
if [ -x .venv/Scripts/python.exe ]; then PY=.venv/Scripts/python.exe; else PY=.venv/bin/python; fi
export PYTHONIOENCODING=utf-8

run_region() {
  local key=$1 data="data/hemibrain_$1.pkl"
  local -a rois groups levels fetch_opts
  case "$key" in
    lh)  # 外側角: 投射ニューロン → 外側角ニューロン → 前大脳（SLP・SMP など）
      rois=("LH(R)"); fetch_opts=(--min-roi-synapses 200 --max-per-type 3)
      groups=("PN=.*PN.*" "LHLN=LHLN.*" "LHN=LH.*" "MBON=MBON.*" "DAN=PAM.*|PPL.*" "other=.*")
      levels=("PN=0" "LHLN=1" "LHN=1" "other=2") ;;
    cx)  # 中心複合体: TuBu → リング → EB の柱状 → FB → 出力
      rois=(EB FB PB "NO"); fetch_opts=(--min-roi-synapses 100 --max-per-type 5)
      groups=("TuBu=TuBu.*" "ER=ER.*|ExR.*" "EBcol=EPG.*|PEG|PEN.*|EL|Delta7.*"
              "FB=FB.*|FS.*|FC.*|FR.*|PFN.*|PFG.*|hDelta.*|vDelta.*" "out=PFL.*|PFR.*|LAL.*" "other=.*")
      levels=("TuBu=0" "ER=1" "EBcol=2" "FB=3" "out=4") ;;
    vis) # 視覚の投射先: 視覚投射ニューロン → AOTU・PVLP の下流
      rois=("AOTU(R)" "PVLP(R)"); fetch_opts=(--min-roi-synapses 200 --max-per-type 3)
      groups=("VPN=LC.*|LPLC.*|LLPC.*|LT.*|MC.*|MeTu.*|aMe.*|LPC.*" "TuBu=TuBu.*" "other=.*")
      levels=("VPN=0" "TuBu=1" "other=1") ;;
    al)  # 触角葉: 嗅覚受容神経 → 投射ニューロン → 局所介在ニューロン
      rois=("AL(R)"); fetch_opts=(--min-roi-synapses 200 --max-per-type 3)
      groups=("ORN=ORN.*" "PN=.*PN[0-9]*$|.*_[a-z]*PN.*" "LN=.*LN.*|lLN.*|vLN.*" "other=.*")
      levels=("ORN=0" "PN=1" "LN=1" "other=2") ;;
    smp) # 上前大脳: MBON の出力先・下行ニューロン・神経ペプチド細胞
      rois=("SMP(R)" "SIP(R)" "SLP(R)"); fetch_opts=(--min-roi-synapses 800 --max-per-type 2)
      groups=("MBON=MBON.*" "DAN=PAM.*|PPL.*" "LHN=LH.*" "DN=DN.*" "CRE=CRE.*" "other=.*")
      levels=("MBON=0" "DAN=0" "LHN=0" "CRE=1" "other=1" "DN=2") ;;
    lal) # 側方付属葉・食道下帯: 前運動。下行ニューロンの手前
      rois=("LAL(R)" "VES(R)" "GOR(R)"); fetch_opts=(--min-roi-synapses 300 --max-per-type 3)
      groups=("CX=PFL.*|PFR.*|FB.*|ER.*|ExR.*|EPG.*|PEN.*|PEG.*" "LAL=LAL.*" "DN=DN.*" "other=.*")
      levels=("CX=0" "LAL=1" "other=1" "DN=2") ;;
    wed) # 楔形領域・前帯状部: 聴覚・機械感覚（求愛歌など）
      rois=("WED(R)" "SAD" "CAN(R)"); fetch_opts=(--min-roi-synapses 400 --max-per-type 2)
      groups=("JO=JO.*" "AN=AN.*|WED.*" "DN=DN.*" "other=.*")
      levels=("JO=0" "AN=1" "other=1" "DN=2") ;;
    vlp) # 腹外側前大脳: 多感覚の収束。視覚投射ニューロンの下流でもある
      rois=("AVLP(R)" "PLP(R)"); fetch_opts=(--min-roi-synapses 700 --max-per-type 2)
      groups=("VPN=LC.*|LPLC.*|LT.*|MC.*|LLPC.*" "AVLP=AVLP.*" "PLP=PLP.*" "DN=DN.*" "other=.*")
      levels=("VPN=0" "AVLP=1" "PLP=1" "other=1" "DN=2") ;;
    ps)  # 後斜面・食道下帯: 下行ニューロンへの入力。運動指令の手前
      rois=("SPS(R)" "IPS(R)" "EPA(R)" "GNG"); fetch_opts=(--min-roi-synapses 400 --max-per-type 2)
      groups=("VPN=LC.*|LPLC.*|LPT.*|H[1-9].*|VS|HS.*" "CX=PFL.*|PFR.*|LAL.*" "DN=DN.*" "other=.*")
      levels=("VPN=0" "CX=0" "other=1" "DN=2") ;;
    cre) # クランプ・脚間部まわり: キノコ体と中心複合体をつなぐ中継
      rois=("CRE(R)" "ATL(R)" "IB" "SCL(R)" "ICL(R)"); fetch_opts=(--min-roi-synapses 700 --max-per-type 2)
      groups=("MBON=MBON.*" "DAN=PAM.*|PPL.*" "CX=FB.*|ER.*|ExR.*|PFL.*|PFR.*|EPG.*" "DN=DN.*" "other=.*")
      levels=("MBON=0" "DAN=0" "CX=1" "other=1" "DN=2") ;;
    *) echo "unknown region $key"; return 1 ;;
  esac
  if [ ! -f "$data" ]; then
    "$PY" fetch_data.py neuprint --group "${groups[@]}" --rois "${rois[@]}" "${fetch_opts[@]}" --out "$data"
  fi
  "$PY" run_all.py --data "$data" --approaches approach_g_output_sites
  "$PY" type_pairs.py --data "$data" --levels "${levels[@]}" --memory-groups none --top 20
  "$PY" match_drivers.py --data "$data" --top 30
}

for r in "${@:-lh cx vis}"; do
  for k in $r; do echo "===== $k"; run_region "$k"; done
done
