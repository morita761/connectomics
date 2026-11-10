# connectomics — 接しているのにシナプスを作らない神経ペアを探す

ショウジョウバエのコネクトーム 3 個体（hemibrain・FlyWire・male-CNS）を突き合わせて、
**膜は接しているのにシナプスが 1 つも無いニューロンのペア**を探す。
人工接着分子で結合を強制的に作らせる実験の標的を選ぶのが目的。

「近接していれば結合する」という Peters' rule からの逸脱は、
結合の有無を決めているのが場所ではなく分子であることの候補になる。
そこに新しいシナプスを 1 本足せたとき、回路と行動が何を返すかを見たい。

## 探索の条件

標的に求めるものは 5 つ。どれか 1 つでも欠ければ実験にならない。

1. **両方向でシナプスが 0** — 片方向だけでは基準値が立たない
2. **複数個体で 0** — 1 個体の欠落は再構成の取りこぼしかもしれない
3. **膜が実際に接している** — 骨格の近さではなく、EM のセグメンテーションで確かめる
4. **発現系が 2 つ使える** — 上流と下流に別の分子を入れるので GAL4 と LexA が要る
5. **上流がその場所で出力している** — 通り道が重なっただけの接触を外す

この 5 つで、全脳 10 領域の 285 組とキノコ体の 392 組を絞り込んだ。
結果と候補の検討は [`results/`](results/) に、まとめは
[報告書の一覧](results/docs/index.html)（`python build_docs.py` で生成）にある。

## 全体設計

```
            ┌──────────────┐   ┌───────────────┐   ┌──────────────────┐
データ取得   │ neuPrint      │   │ FlyWire        │   │ 合成データ（正解つき）│
            │ (hemibrain)   │   │ (fafbseg)      │   │ common/synthetic  │
            └──────┬───────┘   └──────┬────────┘   └────────┬─────────┘
                   └──────────┬───────┴─────────────────────┘
                              ▼
               common/model.py: Dataset（スケルトン・シナプス・メッシュ・結合表, nm 単位）
                              ▼
     ┌───────────┬───────────┬────────┴──┬────────────┬────────────┐
     ▼           ▼           ▼           ▼            ▼            ▼
 approach_a  approach_b  approach_c  approach_d   approach_e   approach_f
 スケルトン    膜表面の     軸索→樹状   Peters' rule  細胞タイプ   接触パッチ
 KD-tree     接触面積     突起の有向   統計モデル     の一貫性     ごとの解析
     └───────────┴───────────┴─────┬─────┴────────────┴────────────┘
                               ▼
            各アプローチ共通の出力: ペア / 接触量 / シナプス数 / score / flagged
                               ▼
         run_benchmark.py（合成データで AP・適合率・再現率を比較）
         run_all.py（実データで全アプローチを実行し、順位を統合）
```

**拡張のしかた**: `approach_xxx/run.py` に `find_candidates(ds, th, **kw) -> DataFrame` を書くだけで、
`run_benchmark.py` と `run_all.py` が自動で比較対象に含める。
返す表には `score`（大きいほど「近いのに結合が弱い」）と `flagged`（閾値を満たすか）の列を入れる。

## アプローチ

| | フォルダ | 切り口 | 近接の量 | 方向 |
|---|---|---|---|---|
| A | `approach_a_skeleton_kdtree` | スケルトン中心線の点群を KD-tree で近傍探索（ベースライン） | 接触ケーブル長 µm | 無向 |
| B | `approach_b_surface_contact` | 膜表面（メッシュ頂点 or 半径から作るチューブ）どうしの接触 | 接触面積 µm² | 無向 |
| C | `approach_c_axon_dendrite` | 軸索/樹状突起に分け、A の軸索 × B の樹状突起 の近接だけを見る | 軸索–樹状突起 接触長 | 有向 |
| D | `approach_d_peters_null_model` | 接触量から期待シナプス数を回帰し（負の二項）、期待より有意に少ないペアを FDR 補正で検出 | 期待値との差 (p, q) | 有向 |
| E | `approach_e_celltype_consistency` | 同じ細胞タイプの組の**他のペア**の結合密度から期待値を出す（leave-one-out）。「仲間はみんな結合しているのに、このペアだけ結合しない」を検出 | タイプ単位の期待値との差 | 有向 |
| F | `approach_f_contact_patches` | 接触を空間的なまとまり（パッチ）に分け、**同じペアの他のパッチ**と比べてシナプスが不自然に少ないパッチを検出。1 行 = 1 パッチ | パッチ単位の期待値との差 | 有向 |

詳しい考え方と長所・短所は各 `run.py` の冒頭に書いてある。

## 閾値（`common/config.py`、すべて変更可）

| パラメータ | 既定値 | 根拠 |
|---|---|---|
| `weak_synapse_max` | 4 | FlyWire (Dorkenwald et al. 2024) は 5 シナプス以上を既定の「結合」とする |
| `skeleton_contact_um` | 1.0 µm | 中心線どうしの距離。突起の半径 (0.1–0.5 µm) のぶんを見込む。哺乳類の Peters' rule 研究では約 2 µm |
| `surface_contact_um` | 0.2 µm | 膜どうし。シナプス間隙は約 20 nm、メッシュの簡略化誤差を見込む |
| `min_contact_cable_um` | 5 µm | 1 点だけの偶然の接触を除く |
| `min_contact_area_um2` | 2 µm² | 同上（B 用） |
| `fdr_alpha` | 0.05 | D の Benjamini–Hochberg FDR |
| `min_expected_synapses` | 3 | D, E: 「本来はこれくらい結合しているはず」の下限 |
| `min_type_pairs` | 3 | E: タイプ単位の期待値を使うのに必要な、同じタイプの組の他のペア数（足りなければ全体平均） |
| `patch_link_um` | 2 µm | F: 接触点がこの距離でつながっていれば同じパッチ |
| `patch_synapse_radius_um` | 1 µm | F: パッチからこの距離以内のシナプスをそのパッチのものとみなす |
| `min_patch_contact_um` / `min_patch_expected` | 3 µm / 2 | F: 評価するパッチの最小サイズ / 期待シナプス数の下限 |

コマンドラインからは `--set skeleton_contact_um=2.0 weak_synapse_max=2` のように上書きできる。

## 使い方

```bash
# 環境（作成済み）
.venv/Scripts/python -m pip install -r requirements.txt

# 1) 合成データで全アプローチを比較
.venv/Scripts/python run_benchmark.py

# 2) 実データを取得（トークンは環境変数 NEUPRINT_APPLICATION_CREDENTIALS）
.venv/Scripts/python fetch_data.py neuprint --types "KCg-m" "MBON01" --out data/hemibrain_mb.pkl

# 3) 全アプローチを実行して順位を統合
.venv/Scripts/python run_all.py --data data/hemibrain_mb.pkl

# 個別に実行
.venv/Scripts/python approach_c_axon_dendrite/run.py --data data/hemibrain_mb.pkl --compartments navis

# 結果の一覧サイト（一致数・上位何 % で絞り込み → 3D 表示へ）
.\tasks.ps1 server          # PowerShell（Git Bash なら ./tasks.sh server）。中身は .venv/Scripts/python serve.py
.venv/Scripts/python fetch_data.py neuprint-meta --data data/hemibrain_mb.pkl   # 細胞名 (instance) を追加

# ドライバーの照合（初回だけ表を取り込む）
.venv/Scripts/python match_drivers.py --import-elife path/to/elife-94168-supp1-v1.xlsx
.venv/Scripts/python match_drivers.py --data data/hemibrain_mb_all.pkl

# 複数のアプローチで一致したペアを探す
.venv/Scripts/python compare_methods.py --results results/hemibrain_mb --data data/hemibrain_mb.pkl --top-pct 1 --min-methods 4

# 候補ペアを 3D で見る（results/<データ名>/3d/ に HTML を出力。ブラウザで開くと回転・拡大できる）
.venv/Scripts/python view_pair.py --data data/hemibrain_mb.pkl --pair 784204967 673509195
.venv/Scripts/python view_pair.py --data data/hemibrain_mb.pkl --csv results/hemibrain_mb/approach_c_axon_dendrite.csv --top 5
```

## 合成データでの比較結果

合成データは、細胞タイプ（5 細胞 × 16 タイプ）と、タイプの組ごとのシナプスの作りやすさの差を含む
（`common/synthetic.py`）。正解は 2 種類: 結合を弱めたペア（各 seed 36–44 組）と、
普通に結合しているペアの中でシナプスを置かなかったパッチ（各 seed 11–19 個）。

### ペア単位（80 ニューロン × 3 seed、軸索/樹状突起は navis で推定）

| アプローチ | AP | flagged の適合率 / 再現率 | 実行時間 |
|---|---|---|---|
| D 統計モデル | **0.57** | 0.56 / 0.13 | 3.2 s |
| E 細胞タイプの一貫性 | 0.56 | 0.57 / 0.17 | 3.2 s |
| C 軸索→樹状突起 | 0.52 | 0.17 / 0.92 | 3.8 s |
| F 接触パッチ | 0.25 | 0.32 / 0.18 | 4.5 s |
| B 表面接触 | 0.13 | 0.08 / 0.93 | 3.0 s |
| A スケルトン KD-tree | 0.13 | 0.06 / 0.98 | 0.4 s |

160 ニューロン（2 seed）では E 0.62 > C 0.60 > D 0.52、E の flagged の適合率は 0.98。
正解の軸索/樹状突起ラベルを使うと C 0.73 / E 0.72 まで上がる。

### パッチ単位（F のみ）

AP 0.38、flagged の適合率 0.38 / 再現率 0.54。

読み取れること:
- **方向と軸索/樹状突起の区別が最も効く**（A, B → C で大きく改善）。分割の精度も性能を左右する。
- **C は再現率重視、D・E は適合率重視**。C は怪しいペアをほぼ全部含むが、flagged の 8 割以上は
  単に結合しにくいタイプの組。D・E は flagged が少ないが半分以上が正解（160 細胞の E では 98%）。
  実データでは「C で候補を広く拾い、E で絞る」使い方が良さそう。
- **E はデータが大きいほど有利**。タイプの組ごとの期待値は、同じ組のペアが多いほど安定するため
  （80 細胞では D とほぼ同等、160 細胞で D を上回る）。
- **F は別の問いに答えるもの**。ペア単位では弱いが（ペア全体が 0 のときは比べる「他のパッチ」がない）、
  ペア単位の方法では原理的に見つけられない「結合しているペアの中の避けられた部位」を拾える。

## hemibrain での実データ試験（2026-11-08）

対象: hemibrain v1.2.1 の KCg-m 590 個 + MBON01 2 個（`fetch_data.py neuprint --types "KCg-m" "MBON01"`）。
スケルトンのノードは 163 万個、ペア内のシナプスは 207,585 個。

| アプローチ | 実行時間 | 近接ペア | flagged |
|---|---|---|---|
| A スケルトン | 77 s | 166,095（可能な全ペアの 95%） | 119,183 |
| B 表面接触 | 806 s | 152,045 | 71,720 |
| C 軸索→樹状突起 | 107 s | 37,587（有向） | 4,804 |
| D 統計モデル | 106 s | 37,587 | 1 |
| E 細胞タイプ | 106 s | 37,587 | 0 |
| F 接触パッチ | 161 s | 48,126 パッチ | 0 |

### 分かったこと

1. **KC 同士のシナプスは「軸索→樹状突起」ではない。** navis で分けた区画で数えると、
   KC→KC のシナプスの 69% は軸索→軸索、18% は樹状突起→樹状突起で、軸索→樹状突起は 4% だけ。
   C〜F は軸索→樹状突起の接触しか見ないので、KC 同士の結合のほとんどを対象外にしてしまう。
   一方 KC→MBON01 は 99% が軸索→樹状突起で、前提どおり。
   → **細胞の種類によって、どの区画どうしの接触を「シナプスの機会」とみなすかを変える必要がある。**
2. **束になって走る軸索では、中心線の近さはほぼ意味を持たない。** A では KC のほぼ全ペアが「近接」になる。
   それでも接触長とシナプス数には単調な関係はある（Spearman 0.41。接触 100 µm 以上では結合なしは 7%、5 µm 以下では 80%）。
3. **KC→MBON01 は、2 個の MBON01 で全く違う。**
   - 612371421: 樹状突起に触れている KCg-m 588 個すべてと結合している（結合なしは 0、0.22 シナプス/µm）。
     KC→MBON がコンパートメント内でほぼ全対全という既知の性質と合う。
   - 673509195: 286 個の KC と接触するが、269 個とは結合なし（5 µm 以上触れていて結合なしは 34 個、0.02 シナプス/µm）。
     3D で見ると（`view_pair.py`）、KC との接触は MBON の樹状突起の茂みではなく、そこから伸びる細い突起の上にあった。
     同じ KC (784204967) は 612371421 とは樹状突起の茂みの中で 170 µm 接し、30 シナプスを作っている。
4. **D・E・F が何も検出しなかった理由。** KC 同士の有向ペアは期待シナプス数が小さく（中央値 約 1）、
   37,587 ペアで FDR 補正すると有意にならない。E は 2 個の MBON01 を同じタイプとしてまとめるので、
   上の 2 個の違いが平均されて消える。D の唯一の検出（MBON01→MBON01）は、出力/入力総数の非常に大きい細胞に
   回帰式を外挿したための見かけの結果。

### 複数のアプローチで一致したペア（`compare_methods.py`）

各アプローチのスコア上位 1% に、6 つすべてで入ったペアは 3 組。いずれも KCg-m どうしで、両方向ともシナプス 0。

| ペア | A 接触長 | B 接触面積 | 各アプローチでの順位 (A/B/C/D/E/F) |
|---|---|---|---|
| 692818802 – 5813062631 | 149 µm | 44.6 µm² | 28 / 6 / 11 / 199 / 18 / 271 |
| 1142024290 – 5812979197 | 184 µm | 41.6 µm² | 5 / 11 / 46 / 285 / 82 / 237 |
| 540730639 – 571424748 | 127 µm | 24.5 µm² | 97 / 229 / 42 / 291 / 75 / 187 |

5 つで一致したペアは 6 組（A・C・E で 1 位の 722821616 – 847302850 など）、4 つ以上で一致したペアは計 39 組。
一覧は `results/hemibrain_mb/consensus_pairs.csv`、3D 表示は `results/hemibrain_mb/3d/`。

- 6 つで一致した 6 個の KC は、どれもこのデータ内でほかの KC 約 200 個と結合し（シナプス 700–800 個）、MBON01 にも 28–46 シナプスを出している。
  つまり「シナプスをほとんど作らない細胞」ではなく、**普通に結合している細胞どうしが、この相手とだけ結合していない**。
- 722821616 は KC の結合相手が 113 個と少なく、MBON01 とも結合していないので、再構築が不完全な可能性がある（このペアは割り引いて見る）。
- **一致の解釈の注意**: C・D・E・F は同じ軸索–樹状突起の接触計算を使うので、互いに独立ではない（C と E の上位 1% は 302/361 が重なる）。
  実質的に独立なのは「幾何だけの A・B」と「区画を考える C〜F」の 2 系統。上の 3 組は、その両方の上位に入ったペア。
  A と B（上位 1% で 964 組が重なる）は似た結果になり、A・B と C〜F の重なりは小さい（上位 1% で 5–19 組）。

### この結果を受けた変更・次の課題

- E に `--type-key instance` を追加した（左右 _L/_R などを区別したグループで比べる）。
  ローダーも instance を保存するようにしたので、**データを取り直すと使える**（今の pickle には入っていない）。
- C〜F に「どの区画の組を接触とみなすか」の選択肢を追加する（KC 同士なら軸索–軸索、樹状突起–樹状突起も含める）。
- 軸索/樹状突起の分割を 1 回だけ計算して使い回す（今は C〜F がそれぞれ計算し直すので、1 つ約 100 秒かかる）。
- D・E・F の FDR 補正は、細胞タイプの組ごとに分けて行う（全ペアまとめてだと検出力が落ちる）。
- 結果を ROI（γ1–γ5、ペダンクル、カリックスなど）に対応づけて、接触と結合がどこで食い違うかを見る。

## hemibrain 視覚系 × 嗅覚系（2026-11-08）

対象: 視覚系の投射ニューロン（LC・LPLC・LT・LLPC・LPC・MC、93 タイプ 357 個）と、嗅覚系の投射ニューロン
（触角葉 PN、VP 糸球体を除く、163 タイプ 295 個）。LH・PLP・SCL・SLP・CA のいずれかに 20 シナプス以上もつもの、
1 タイプ最大 8 個。

```bash
.venv/Scripts/python fetch_data.py neuprint --group "visual=(LC|LPLC|LT|LLPC|LPC|MC)[0-9].*" "olfactory=(?!VP)[A-Z][A-Za-z0-9]*_[a-z0-9]*PN[A-Za-z0-9]*|M_[a-z0-9]*PN.*" --rois "LH(R)" "PLP(R)" "SCL(R)" "SLP(R)" "CA(R)" --min-roi-synapses 20 --max-per-type 8 --out data/hemibrain_vis_olf.pkl
.venv/Scripts/python run_all.py --data data/hemibrain_vis_olf.pkl --compartments navis
.venv/Scripts/python type_pairs.py --data data/hemibrain_vis_olf.pkl --cross visual olfactory
```

### 分かったこと

1. **視覚系と嗅覚系は、ほとんど接触しない。** 5 µm 以上接触するニューロンペアは 184 組だけ
   （同じ系どうしは数万組）。結合も少ない（視覚→嗅覚 25 組・48 シナプス、嗅覚→視覚 63 組・322 シナプス）。
2. **ただし、接触した場合はよく結合する。** 接触 1 µm あたり 0.124 シナプスで、同じ系どうし（約 0.03）の 4 倍。
   同じ系どうしの接触の多くは束になって並走する軸索で、シナプスを作らないため。
   → 期待値は、データ全体ではなく**視覚系×嗅覚系のペアの密度を基準に**計算するようにした（`type_pairs.py`）。
3. **接触のわりに結合しないタイプの組（無向、期待値 5 以上・観測 0）は 3 組。**

| 視覚系 | 嗅覚系 | 細胞数 | 接触合計 | シナプス | 期待値 | p (Poisson) | 接触している部位 |
|---|---|---|---|---|---|---|---|
| LT35 | M_smPN6t2 | 1 × 1 | 84 µm（68 µm²） | 0 | 11.4 | 1e-5 | PN の**細胞体線維** × LT35 の樹状突起 |
| LT52 | M_l2PNl21 | 8 × 1 | 80 µm（19 µm²） | 0 | 10.8 | 2e-5 | PN の**軸索** × LT52 の軸索・細胞体線維・樹状突起（11 µm） |
| LT36 | M_smPNm1 | 2 × 1 | 39 µm（8 µm²） | 0 | 5.2 | 0.006 | **樹状突起どうし** |

   - LT52 × M_l2PNl21 は、サンプルした LT52 の 8 個すべてがこの PN の軸索に接しているのに、どれも結合していない。
     ただし LT52 の樹状突起との接触は 11 µm だけ。
   - LT35 × M_smPN6t2 の接触は、PN の細胞体線維（シナプスをふつう作らない部分）の上。通り道が重なっているだけの可能性が高い。
   - 116 組を調べているので、多重比較を考えると p が十分小さいのは上の 2 組だけ。
4. **軸索→樹状突起（C、有向）では、期待値 5 以上で結合しない組はなかった**（最大でも M_l2PN10t19 → LT36 の期待値 4.3、p = 0.01）。
   つまり「視覚系の軸索が嗅覚系の樹状突起に広く接するのに結合しない」（またはその逆）という、はっきりした組は見つからなかった。
5. hemibrain は脳の右半分が中心で、視葉の大部分は入っていない。今回見ているのは中心脳の中の接触だけ。
   instance 名の `_L` は細胞体が左側にあり、右側に投射している細胞。

### 同じ系どうしの組（同じデータ）

`type_pairs.py` は、グループの組（嗅覚×嗅覚、視覚×視覚、視覚×嗅覚）ごとに、それぞれの基準の密度で期待値を計算する。
flagged（期待値 5 以上・観測が期待値の 0.2 倍以下）の数:

| 組 | 無向（A） | 軸索→樹状突起（C） |
|---|---|---|
| 嗅覚 × 嗅覚 | 881 / 8,278 | 291 / 1,054 |
| 視覚 × 視覚 | 142 / 2,519 | 173 / 1,573 |
| 視覚 × 嗅覚 | 3 / 116 | 0 / 39 |

軸索→樹状突起で接触が長いのに結合 0 の例: M_vPNml64 → M_vPNml53（252 µm）、VM5v_adPN → VM2_adPN（181 µm）、
LC14 → LT81（500 µm）、LT52 → LC37b（444 µm）。無向の上位は、多糸球体 PN（M_vPNml…）どうしや隣り合う LC タイプどうしなど、
同じ神経路を並走する組が多い（束になった軸索の可能性）。

### 次の課題
- 1 タイプ最大 8 個の制限を外す・ROI を広げる（AOTU・PVLP・LO など）。
- male-cns:v1.0（neuPrint で利用可）なら視葉も含めて全体を見られる。
- 接触部位の区画（軸索・樹状突起・細胞体線維）を type_pairs の表に自動で入れる。
- 3 組のタイプのドライバーは未登録（サイトの「調べる」リンクから確認して drivers.csv に追記）。

## hemibrain キノコ体まわりの全タイプ（2026-11-08）

対象: キノコ体 (MB(R)) に 50 シナプス以上もつニューロン（227 タイプ、1 タイプ最大 10 個、計 586 個）。
グループ: PN 121・KC 102・DAN 205・MBON 54・MBIN（APL・DPM・OA）8・LHN 34・その他 62。

```bash
.venv/Scripts/python fetch_data.py neuprint --group "PN=.*PN.*" "KC=KC.*" "DAN=PAM.*|PPL.*" "MBON=MBON.*" "MBIN=APL.*|DPM.*|OA-.*" "LHN=LH.*" "other=.*" --rois "MB(R)" --min-roi-synapses 50 --max-per-type 10 --out data/hemibrain_mb_all.pkl
.venv/Scripts/python run_all.py --data data/hemibrain_mb_all.pkl --approaches approach_a_skeleton_kdtree approach_g_output_sites
.venv/Scripts/python type_pairs.py --data data/hemibrain_mb_all.pkl
.venv/Scripts/python inspect_pairs.py --data data/hemibrain_mb_all.pkl --kind output --top 40
```

### 新しく追加したもの
- **アプローチ G（出力部位→相手）**: 上流の T-bar から 3 µm 以内のケーブル（= 出力部位）が、下流のどこかに 1 µm 以内で接する長さ。
  受け手の区画を仮定しないので、DAN → KC の軸索のような結合も扱える。hemibrain_mb では全シナプスの 99.99% が
  「接しているペア」に含まれた（C は 27%）。合成データでは不利（合成データのシナプス前部位は結合した場所にしかないため）。
- **type_pairs.py**: 有向の集計を「グループの順序付きの組」ごとの基準密度で行い、上流→下流（feedforward）と
  記憶系（KC・MBON・DAN・MBIN を含む）の印と優先度を付ける。階層は PN 0 → KC/DAN/MBIN 1 → MBON/LHN 2。
- **inspect_pairs.py**: 代表ペアの接触がどの ROI（キノコ体の区画など）で起きているか、そこで上流が出力し下流が入力を受けているかを調べる。
- C（navis による軸索/樹状突起の分割）は、APL などシナプス数万の細胞で 100 GB 級のメモリを使ったため、
  このデータでは中止した。巨大な細胞では navis を使わないようにした（`MAX_NAVIS_SYNAPSES`）。

### 結果
上流→下流・記憶系・接触のわりに結合しない（G、期待値 5 以上・観測が期待値の 0.2 倍以下）タイプの組は 243 組。
接触部位と比較対象を確かめた上位候補:

| 上流 → 下流 | 接触（G） | シナプス / 期待値 | 接触している領域 | 同じ下流への比較 |
|---|---|---|---|---|
| **PAM01_b → MBON29** | 2,049 µm | 5 / 161 | γ5 93% | PAM08_b → MBON29 は 0.161 シナプス/µm、PAM01_b は 0.002（約 80 倍の差）。同じ γ5 の PAM01_a は 0.041 |
| **PAM15_b → MBON04** | 249 µm | 0 / 20 | β'2 89% | PAM05 0.131、PAM08_c 0.176 に対し 0 |
| **PAM07 → MBON01** | 223 µm | 0 / 18 | γ4 74%、γ5 26% | PAM02 0.095 に対し 0 |
| **APL → MBON33** | 840 µm | 2 / 36 | γ2 62%、γ3 23% | PPL103 0.235、APL は 0.002 |
| DA4m_adPN → KCg-t | 33 µm | 0 / 20 | カリックス 100% | PN → KC の基準 0.61 シナプス/µm |
| VP1d+VP4_l2PN1 → KCg-t | 21 µm | 0 / 13 | カリックス 100% | 同上 |

- いずれも、接触している領域で上流は出力を出し（例: PAM01_b は γ5 に 108 個の T-bar）、
  下流も入力を受けている（MBON29 は γ5 で 2,781 個）。通り道が重なっているだけではない。
- 特に **DAN → MBON** は「同じ区画で、同じ MBON に対して DAN のサブタイプによって結合の有無がはっきり分かれる」パターン。
- 除外・保留したもの: 温度・湿度系 PN（VP1m+_lvPN など）→ KC は、PN がその領域でほとんど出力していない（カリックスで 15 個）。
  PN → OA-VUMa2 は触角葉での接触で、記憶系の中心ではない。PPL101・PAM13 → MBON は接触の大半が CRE・SMP（キノコ体の外）。
- 3D: `results/hemibrain_mb_all/3d/`（上の 6 組）。一覧はサイトの hemibrain_mb_all ページ（接触の種類「出力部位→相手」、「上流→下流のみ」「記憶系」で絞り込み）。
- ドライバー（2026-11-08 調査。[Shuai et al. eLife 94168](https://elifesciences.org/articles/94168) の補足表 1・2 と
  [Rubin & Aso eLife 90523](https://elifesciences.org/articles/90523) による。詳細は `drivers.csv`）:

| 候補 | 上流（GAL4） | 下流（LexA） | 評価 |
|---|---|---|---|
| **PAM07 → MBON01** | **MB197B**（split-GAL4、PAM-γ1γ2>γ4、1–2 個、疎） | **SL86047**（split-LexA、MBON01 のみ、特異性 2） | **上流 GAL4・下流 LexA がそろう**。比較対照に PAM02（MB106C）→ MBON01（0.095 シナプス/µm） |
| PAM01_b → MBON29 | MB313B（PAM-γ5、4–6 個）／MB315C。サブタイプ（a/b）の区別は不明 | MBON29 は split-GAL4 が 2 系統（線名未確認）。LexA なし | サブタイプの対応が課題 |
| PAM15_b → MBON04 | 見つからず | SL86063（split-LexA、MBON01・03・04 を一緒に標識） | 上流のドライバーなし |
| APL → MBON33 | VT043924-GAL4（APL-Gal4） | SS75200・SS77383（split-GAL4）。LexA なし | 下流の LexA なし |

  補足表 2 には、split-GAL4 の DBD 側を LexA に置き換えるための LexADBD 系統（R27G01・R15B01・R17C11 など）がある。
  論文では変換した 34 系統中 22 系統で同じ発現パターンが得られたとされる。MB313B（R58E02-AD + R17C11-DBD）は
  17C11-LexADBD があるので split-LexA にできる可能性がある。

### ドライバーの網羅的な照合（`match_drivers.py`）

eLife 94168 の補足表 1（585 系統、うち split-LexA 20 系統）を `data/driver_lines.csv` に取り込み、
G の flagged 2,462 組すべてと照合した。上流に GAL4 系・下流に LexA 系がある組は 79 組、そのうち上流→下流は 16 組。
16 組について接触領域と「同じ下流タイプへの、同じグループの他の上流タイプの密度」を確かめた
（`results/hemibrain_mb_all/driver_candidates_checked.csv`）。

| 上流 → 下流 | 接触 | シナプス / 期待値 | 接触領域 | 上流の出力 / 下流の入力（主領域） | 比較（同じ下流へ） | 上流 GAL4 | 下流 LexA |
|---|---|---|---|---|---|---|---|
| **PAM07 → MBON01** | 223 µm | 0 / 18 | γ4 74%、γ5 26% | 140 / 1,505 | PAM06_b 0.113 | MB197B | **SL86047**（MBON01 のみ） |
| **PPL103 → MBON11** | 199 µm | 1 / 16 | γ2 73% | 1,974 / 99 | PPL101 0.181 | MB302B（特異性 2） | **SL86058**（MBON11 のみ） |
| **PAM08_e → MBON01** | 1,222 µm | 11 / 96 | γ5 55%、γ4 42% | 21 / 23,798 | PAM06_b 0.113 | MB518C（+PAM11、発現弱い） | **SL86047** |
| **PAM08_e → MBON04** | 467 µm | 4 / 37 | γ4 90% | 97 / 2,139 | **PAM08_c 0.176**（同じ PAM08 の別サブタイプ） | MB518C | SL86063（MBON01・03・04） |
| PPL101 → MBON04 | 714 µm | 4 / 56 | SMP 63%、CRE 36% | 259 / 3,113 | PAM08_c 0.176 | MB503B | SL86063 |

- **PPL103 → MBON11**: PPL103 は γ2 に T-bar が 1,974 個あるのに、MBON11 にはほぼ出力しない。ただし MBON11 の γ2 での入力は 99 個と少ない（区画の境界付近の可能性）。
- **PAM08_e → MBON04**: 同じ γ4 で、PAM08_c は MBON04 に密に結合するのに PAM08_e は結合しない。PAM01_a/b と同じ「サブタイプで分かれる」型。ただし MB518C は PAM08 のサブタイプを区別しない。
- 除外: 接触の大半がキノコ体の外（CRE・SMP・SIP）で上流の出力が少ないもの（PAM13 → MBON12/04/01、PPL104/105 → MBON04）、
  接触が今回の ROI 一覧の外（左半球など）にあるもの（PAM12・PAM07 → MBON04、PAM11・OA-VPM3 → MBON03、PPL102 → MBON04）。
- 3D: `results/hemibrain_mb_all/3d/` に追加（PPL103 × MBON11、PAM08_e × MBON04・MBON01、PPL101 × MBON04）。

### 行動実験の文献（`lit_search.py`）

Europe PMC の API で、細胞名（LM 名）・系統名 AND Drosophila を検索し（各 100 本まで）、オープンアクセス論文は本文、
それ以外は要旨から「細胞名・操作・行動」がそろう文を拾った（`results/literature/<タイプ>.csv`・`summary.md`）。
主要論文は本文でも確認した。

| タイプ（LM 名） | 行動の知見 | 使われた系統 | 単独タイプでの行動実験 |
|---|---|---|---|
| MBON01（MBON-γ5β'2a） | 活性化で忌避（Aso 2014、Owald 2015）、阻害で匂いの忌避が誘引に変わる（Owald 2015）、覚醒促進（睡眠）、甘味記憶・食物探索 | MB011B（MBON01・03・04 をまとめて）、VT1211、R21D02 | **ほぼ無し**（MBON01 だけの SS01308・SL86047 を行動に使った論文は見つからず） |
| MBON11（MBON-γ1pedc>α/β、MVP2） | 活性化で誘引（MB112C; Aso 2014、Tsao 2018）、阻害で 2 時間の嫌悪記憶が低下（MB112C・MB085C）、空腹時の食欲記憶の発現に必要（Perisse 2016） | MB112C、MB085C、MB262B | **あり**（MB112C は MBON11 にほぼ特異） |
| MBON04（MBON-β'2mp_bilateral） | MBON01 と同じ（MB011B でまとめて操作）、覚醒促進 | MB011B | 無し |
| PPL103（PPL1-γ2α'1、MV1） | 匂いと対にした活性化で嫌悪記憶（1 日・4 日、MB099C; Aso & Rubin 2016）、記憶の固定に必要（MB296B） | MB099C、MB296B | あり |
| PAM07（PAM-γ4<γ1γ2） | 該当する行動の記述なし（12 本） | — | 無し |
| PAM08（PAM-γ4） | 活性化で睡眠が減る（MB312B） | MB312B | あり（睡眠） |

- **注意 1**: split-LexA（SL 系統）は解剖（CsChrimson-mVenus の発現）でしか確かめられておらず、行動実験の報告はない。
  変換した 34 系統のうち 12 系統は、発現が弱いか不安定だった（Shuai et al. 2024）。行動に使う前に、LexAop の発現強度を確かめる必要がある。
- **注意 2**: MBON11 の樹状突起は γ1・pedc にあるので、γ2 の接触（PPL103 → MBON11）は樹状突起の本体から外れている可能性が高い。

## キノコ体以外の領域（2026-11-08、`run_regions.sh`）

外側角（LH）・中心複合体（CX）・視覚の投射先（AOTU・PVLP）で、同じ流れ（取得 → G → タイプの組 → ドライバー照合）を回した。
グループと上流→下流の階層は `run_regions.sh` に定義。候補は `rank_candidates.py` で
「同じ下流タイプに、同じグループの他の上流タイプがどれだけ密に結合しているか（comp_rate）」と並べ、
上位を `inspect_pairs.py --ranked --any-region` で領域ごとに確かめた（`results/hemibrain_<領域>/candidates_ranked.csv`・`inspect_output.csv`）。

```bash
./run_regions.sh lh cx vis
.venv/Scripts/python rank_candidates.py --data data/hemibrain_lh.pkl
.venv/Scripts/python inspect_pairs.py --data data/hemibrain_lh.pkl --ranked --top 30
```

| 領域 | ニューロン | 上流→下流で flagged | 主な型 |
|---|---|---|---|
| LH（LH に 200 シナプス以上、1 タイプ 3 個まで） | 1,502 | 1,781 | PN → LHN。接触は LH の中（100%）、両側ともそこで出力・入力している |
| CX（EB・FB・PB・NO に 100 以上、5 個まで） | 1,018 | 361 | FB の層状ニューロン → PFL・PFR、ExR7 → PEN・PEG（EB） |
| 視覚（AOTU・PVLP に 200 以上、3 個まで） | 979 | 509 | LC・LT → PVLP・AVLP の下流。LC は視覚糸球体（PVLP 100%）の中での接触 |

主な候補:

| 上流 → 下流 | 接触 | シナプス / 期待値 | 領域 | 比較（同じ下流へ、同じグループの別タイプ） |
|---|---|---|---|---|
| **ExR7 → PEN_a / PEN_b / PEG** | 2,648 / 3,138 / 1,808 µm（20 ペア） | 5 / 271、5 / 321、5 / 185 | EB 96–100% | ER3a_c → PEN_a 0.418、ER6 → PEN_b 0.183 シナプス/µm |
| FB4Z・FB5R・FB5S → PFL2 / PFL3 | 600–1,700 µm（25 ペア） | 1–5 / 40–110 | FB 100% | FB3B → PFL2 0.521、PFNp_a → PFL3 0.314 |
| VA1v_adPN → LHAV2m1 | 96 µm（5 ペア） | 0 / 11.5 | LH 100% | DP1l_adPN 0.752 |
| VA1d_adPN → LHPD2a2 | 84 µm（8 ペア） | 0 / 10.1 | LH 100% | VA4_lPN 0.363 |
| DM2_lPN → LHAV4a7 | 170 µm（6 ペア） | 0 / 20.2 | LH 100% | VM7d_adPN 0.410 |
| VC3m_lvPN → LHAV4g13 / LHPV4a2 | 146 / 116 µm | 0 / 17.4、0 / 13.8 | LH 100% | VP2_adPN 0.380、VL2a_adPN 0.386 |
| LC16 → AVLP466 | 99 µm（6 ペア） | 0 / 12.7 | PVLP 100% | LT67 0.355 |
| LC6 → PVLP058 | 65 µm（8 ペア） | 0 / 8.4 | PVLP 100% | LT87 0.224 |
| LC21 → PVLP096、LC31 → PVLP130、LC18 → PVLP094 | 60–114 µm | 0 / 8–15 | PVLP 87–100% | 0.33–0.50 |
| LT52 → TuTuA | 595 µm（9 ペア） | 13 / 77 | AOTU 100% | LT55 0.250 |

- LH の PN → LHN は 1 組あたりの接触が 100 µm 前後と小さいが、細胞ペアの多く（5–9 ペア）で一貫して 0。
  VA1v・VA1d はフェロモン（性行動）の糸球体。
- CX の FB は 1 つの ROI として扱ったので、FB の層（FBl1–9）までは見ていない。層の境界で触れているだけの可能性がある。
- 視覚の LC は 1 タイプ数十個あるうちの 3 個だけを見ているので、タイプ全体での確認が必要。LT（大きな接線ニューロン）の組は、
  主な接触が AVLP など投射先の外にあるものが多く、優先度は低い。

### キノコ体以外のドライバー照合

3 つの split-GAL4 集を `match_drivers.py` に取り込んだ（`data/driver_lines_{cx,lh,lc}.csv`）:
[Wolff et al. eLife 104764](https://elifesciences.org/articles/104764)（CX、593 系統・234 タイプ、hemibrain と同じ名前）、
[Dolan et al. eLife 43079](https://elifesciences.org/articles/43079)（LH、225 系統・81 タイプ。名前に LH を付けて hemibrain に合わせ、
hemibrain のサブタイプ（_a など）は区別しない）、[Wu et al. eLife 21022](https://elifesciences.org/articles/21022)（LC・LPLC、54 系統）。
PN 側は eLife 94168 の PN 系統を使う。

```bash
.venv/Scripts/python match_drivers.py --import-cx elife-104764-fig2-data1-v1.xlsx --import-lh elife-43079-supp2-v1.xlsx --import-lc elife-21022-supp1-v2.docx
```

**3 つの集とも split-LexA は無い**ので、上流 GAL4・下流 LexA がそろう組は 0 のまま（両側 GAL4 は LH 187・CX 320・視覚 6 組）。
eLife 94168 の LexADBD 系統（16 種類のエンハンサー）と、各系統の DBD 側を突き合わせたが、候補の下流タイプで split-LexA に変換できるものは無かった
（変換できそうなのは FR1・FC1・FB5F/Q・WED-PN1 の系統だけ）。

両側に特異的な GAL4 系統がある候補:

| 上流 → 下流 | 接触 / シナプス / 期待値 | 領域 | 比較 | 上流 GAL4 | 下流 GAL4 |
|---|---|---|---|---|---|
| **ExR7 → PEN_a / PEN_b / PEG** | 1,800–3,100 µm / 5 / 185–321 | EB | ER3a_c → PEN_a 0.418 | SS65157・SS82313 | PEN_a: SS62862、PEN_b: SS54295、PEG: SS02191 |
| **VL1_vPN → LHAV1a3** | 127 µm / 0 / 15.1（3 ペア） | LH 98% | VC4_adPN 0.60 | **SS86025**（1 個、特異性 2） | **LH1543**（Ideal） |
| **DL3_lPN → LHAV4b4** | 91 µm / 0 / 10.8（8 ペア） | LH 100% | DC4_adPN 0.35 | SS01147（特異性 2） | LH240（Good） |
| VL2a_vPN → LHPV2a1_c / _d | 97・149 µm / 0 / 11.5・17.8 | LH 100% | VA1v_vPN 0.44 | SS01850（特異性 2） | LH86（Ideal。LHPV2a1+LHPV2b1、サブタイプ未区別） |
| DM5_lPN → LHAV4g8 | 77 µm / 0 / 9.2 | LH 98% | DL5_adPN 0.36 | SS49306（特異性 2） | SS01836（特異性 1、他も標識） |
| FB4Z → PFL2 | 1,553 µm / 1 / 103 | FB | FB3B 0.52 | SS63390・SS95748 | SS82329・SS82331（他の CX タイプも） |
| LC16 → PVLP006 | 95 µm / 2 / 12.3 | PVLP 100% | LT82b 0.22 | OL0046B（行動実験で後退歩行） | SS48866（PVLP004+006） |

- 除外: VP1m+_lvPN → LHPV12a1（接触の 45% がカリックスで、上流はそこで 23 個しか出力しない）。

## 候補を広げた再探索（2026-11-09）

FlyWire のタイプ単位の表が手元にある（`data/flywire_type_edges.csv.gz`）ので、
**全候補を安く一次スクリーニングできる**ようになった。そこで手順を変えた。

```bash
# 1) FlyWire で全候補を一次スクリーニング（期待値の下限を 8 → 3 に下げた。ローカルなので安い）
python flywire_check.py --data data/hemibrain_*.pkl --min-expected 3 --top 2000
# 2) hemibrain と FlyWire の両方で 0 だった組だけ male-cns に問い合わせる（高い処理を絞る）
python validate_cross.py --data data/hemibrain_lh.pkl --pairs results/zero_hemibrain_flywire.csv
# 3) 3 個体すべてで 0 の組の接触と領域を確認
python validate_cross.py --data data/hemibrain_lh.pkl --stage contact --pairs results/zero_all3.csv
python inspect_pairs.py --data data/hemibrain_lh.pkl --pairs results/shortlist_zero.csv
```

| 段階 | 組数 |
|---|---|
| hemibrain でシナプス 0 | 1,427 |
| FlyWire でも 0（照合できた 1,089 組のうち） | 140 |
| male-cns でも 0 = **3 個体すべて 0** | **58** |
| 別個体でも接触が再現（10 µm 以上） | **36** |
| 上流も下流もその領域で働いている | **31** |
| 両側にドライバーがある | 3 |

以前は「3 個体すべて 0 かつ接触再現」が 2 組だけだったが、**36 組**になった。

### 最終候補（`results/shortlist_final.csv`）

| 上流 → 下流 | 領域 | hemibrain 接触 | male-cns 接触 | 上流出力 / 下流入力 | 上流ドライバー | 下流ドライバー |
|---|---|---|---|---|---|---|
| **DL4_adPN → LHPV6c1** | LH 85% | 41 µm | **27–37 µm** | 917 / 925 | SS10228（+VA1d_adPN、特異性 1） | **LH2220**（特異性 1） |
| **M_vPNml87 → LHAV2m1** | LH 100% | 96 µm | **41–48 µm** | 203 / 1,083 | SS00587・SS00586（特異性 1） | なし |
| **V_l2PN → LHAV3h1** | LH 89% | 58 µm | **47 µm**（左右ほぼ同じ） | 789 / 489 | SS01336・SS02076（特異性 1） | なし |
| DC2_adPN → LHAV5a8 | LH | 73 µm | 29–36 µm | 674 / 263 | なし | なし |
| M_imPNl92 → LHAV2k8 | LH | 63 µm | 31–36 µm | 341 / 917 | SS01097・SS01092 | なし |
| M_l2PNl22 → LHAV3e1 | LH | 63 µm | 31–41 µm | 153 / 355 | SS32193（特異性 2） | なし |
| LT79 → AVLP372 | AVLP | 54 µm | 44–50 µm | 357 / 723 | なし | なし |
| M_vPNml65 → LHCENT4 | LH | 52 µm | 28–42 µm | 137 / 1,183 | SS10947（特異性 2） | なし |

- **両側にドライバーがそろう中では DL4_adPN → LHPV6c1 が最有力**。3 個体すべてで両方向 0 個、
  接触は別個体でも 27–37 µm（元 41 µm）と最もよく再現し、上流 917・下流 925 と両者が外側角で活発。
- **M_vPNml87 → LHAV2m1** は接触が一番大きく（96 µm、別個体で 41–48 µm）、接触が外側角に 100% 収まる。下流のドライバーが無い。
- **V_l2PN → LHAV3h1** は左右でほぼ同じ接触長（47.4・47.7 µm）という再現性。

### 除外したもの（ドライバーはあるが条件を満たさない）

| 組 | 理由 |
|---|---|
| **ER3a_d → FB4L** | 接触は LAL 98% だが、**上流がそこで 6 個しか出力していない**。放出できない場所での接触 |
| **PPL105 → MBON04** | 接触がキノコ体の外（SIP 95%）。上流の出力も 139 個と少ない |
| FB6L → PFL2・PFL3 | 3 個体ともシナプス 0 だが、別個体で接触がほぼ消える（0–14 µm） |
| PFNv → PFL3 | 接触 938 µm で最大だが、male-cns では結合している |

### 触角葉（AL）の探索

触角葉（566 ニューロン・295 タイプ）も同じ手順で回した（`results/hemibrain_al/`）。
FlyWire で 0 の組は 19 組だが、**male-cns で照合できたのは 4 組だけ**（AL の PN のタイプ名が揃わない）。
3 個体すべて 0 は DA4m_adPN → ALBN1（116 µm）の 1 組。

接触の大きい組（hemibrain と FlyWire で 0。male-cns 未照合）:
VA3_adPN → M_l2PNl20（483 µm）、VA5_lPN → M_l2PNl20（351 µm）、VL1_ilPN → M_lv2PN9t49（224 µm）。

**ORN（嗅覚受容神経）→ PN の組**も出た（ORN_VC3l → M_vPNml87 など 9 組）。ORN は glomerulus ごとに
Or-GAL4 があり **LexA 版が存在するものもある**ので、ドライバーの面では有望。ただし接触は 20–39 µm と小さい。

## ドライバーと読み出しの調査（2026-11-09）

候補の細胞について、ドライバー系統の中身と、行動・生理の読み出しがあるかを調べた
（`lit_search.py` の TARGETS に追加。結果は `results/literature/`）。

### 行動の読み出しがある下流タイプ

| 下流 | 読み出し | 出典 | 系統 |
|---|---|---|---|
| **LHAV1a1** | **活性化で忌避**（PI −0.23）。**抑制でジオスミン（忌避臭）への産卵反応が消える** | Dolan 2019、Huoviala 2018 | LH728（Good）、LH1983（Ideal） |
| **PFL3** | 目標指向の操舵の出力。頭方向の信号を操舵指令に変換する | Westeinde 2024 | SS82335、SS60291 |
| **PFL2** | 同（P2X2 で活性化した実験がある） | Westeinde 2024 | SS82329、SS82331 |
| LHPD2a1 | 学習した匂いの想起に必要 | Dolan 2018 | LH989・LH991（Good） |
| LHPV5a1・LHPV6a3・LHAD1b2 | 活性化で誘引（PI +0.38・+0.42・+0.53） | Dolan 2019 | LH291・LH2193・LH1538 |
| **LHAV1a3**（最有力候補の下流） | 光遺伝学で**検査済みだが有意差なし**（PI 0.22） | Dolan 2019 | **LH1543（Ideal、AD 76E07 / DBD 66B12）** |

LHAV1a3 は「系統として行動実験に使える実績はあるが、単独で活性化しても行動が変わらない」。
強制シナプス形成の結果を読むには、カルシウムイメージングなどの生理の読み出しが必要。

### 下流のドライバーが無いことの確認

候補の下流タイプにドライバーが見つからないのは、名前の対応づけの失敗ではなく**本当に公開されていない**ためだった。
Dolan 2019 の補足表に載っている外側角のタイプは **66 種類だけ**で、
LHAV2m1・LHAV3h1・LHAV5a8・LHAV2k8・LHAV3e1・LHCENT4・LHPD5d1・LHAD1c3 などはどれも含まれていない。
候補の下流で収録されているのは **LHPV6c1（LH2220）と LHAV1a3（LH1543）の 2 つだけ**。

### 最有力候補の上流ドライバーの確認（DL4_adPN → LHPV6c1）

SS10228 は DL4_adPN と **VA1d_adPN** を一緒に標識する。そこで一緒に標識される細胞も調べた:

| 組 | hemibrain | FlyWire | 判断 |
|---|---|---|---|
| DL4_adPN → LHPV6c1 | 接触 41 µm・シナプス 0 | 0 | 本命 |
| VA1d_adPN → LHPV6c1 | 接触 6.7 µm・シナプス 0 | 0 | **一緒に標識されても結合しない**（問題なし） |
| **VA1v_adPN → LHPV6c1** | 接触 15 µm・**シナプス 5**（密度 0.33） | **32** | **結合する**。SS10228 の注記に「VA1v_adPN の可能性あり」とあり、**要確認** |

VA1d と VA1v の区別が付いていない可能性があるので、実験前に SS10228 の標識像を確認する必要がある。
なお LHAV2m1 は M_vPNml87（96 µm）と M_vPNml86（66 µm）の**どちらからも 0 個**で、
2 つの独立した上流が結合しないという点では筋が通っている。

### 2 つの発現系（GAL4 と LexA）がそろうか

上流と下流に別々の人工分子を入れるには、発現系が 2 つ要る。取り込んだ全系統 1,457 のうち
**split-LexA は 20 系統だけ**（すべて eLife 94168 のキノコ体向け）。
3 個体で確かめた 106 組のうち、**上流 GAL4 + 下流 LexA がそろうのは 1 組だけ**だった。

| 組 | 上流 | 下流 | 判定 |
|---|---|---|---|
| **PPL105 → MBON04** | MB059B〈GAL4〉（特異性 1） | **SL86063〈LexA〉**（MBON01・03・04、特異性 2） | **◎ 2 つの発現系がそろう** |
| DL4_adPN → LHPV6c1 | SS10228〈GAL4〉 | LH2220〈GAL4〉 | ○ どちらも GAL4。2 つ目の系が要る |

LexA がキノコ体に偏っているため、外側角・中心複合体・視覚の候補はすべて「両方 GAL4」か「片側のみ」。

eLife 94168 の LexADBD 系統（エンハンサー 16 種類）と候補の全系統の DBD 側を突き合わせたが、
**新たに LexA 化できる候補は無かった**。
参考: LHAV1a3（LH1543）と LHAV1a1（LH728・LH1983）は AD が同じ **76E07**。

### 系統の注意点

| 系統 | 注意 |
|---|---|
| SS71548（FB5AA） | 表に「FB5AA かどうか確実でない」と記載 |
| SS75533・SS69168（FB8F_a） | FB8F_a と _b の**両方**を標識 |
| SS82335（PFL3）・SS60291 | 「PFL2 も否定できない」「LAL の枝から PFL2 か PFL3 かどちらとも」 |
| SS86025（VL1_vPN） | LM 1 細胞、特異性 2、主要ドライバー。最もきれい |
| LH1543（LHAV1a3） | Ideal 評価。LHAV1a3 は出力ニューロン、アセチルコリン＋GABA、約 5.7 細胞 |

## 3 個体での再現確認（2026-11-09、`validate_cross.py` + `flywire_check.py` → `shortlist.py`）

使えた別個体は次の 2 つ。これに hemibrain を加えた **3 個体**で照合する。

| データ | 個体 | 使えるか |
|---|---|---|
| hemibrain:v1.2.1 | 雌・右半球 | 元データ |
| **male-cns:v1.0** | 雄・全神経系 | **使える**（v0.9 は同じ個体の古い版） |
| **flywire:783** | 雌・全脳 | **使える**（Zenodo 10676866 の結合表 852 MB ＋ GitHub の注釈 32 MB） |
| optic-lobe:v1.1 | 別個体・視覚系 | **使えない**。視覚葉の中のシナプスしか公開されておらず、中枢脳の結合は 0 件に見える |
| manc:v1.2.3 | 雄・腹部神経索 | 対象外（脳がない） |
| mushroombody | 別個体・キノコ体 | 使わない。17 タイプだけで名前の体系が違い（MBON-06・PPL1-05 など）、対応づけが推測になる |

```bash
# FlyWire は 1 回だけ集計してタイプ単位の小さい表にする（data/flywire_type_edges.csv.gz、19 MB）
python flywire_check.py --build --connections <...>.feather --annotations <...>.tsv
python flywire_check.py --data data/hemibrain_cx.pkl data/hemibrain_lh.pkl
python validate_cross.py --data data/hemibrain_cx.pkl --datasets male-cns:v1.0 --top 150
python shortlist.py --max-total 2   # → 最終候補
python validate_cross.py --data data/hemibrain_cx.pkl --stage contact --pairs results/shortlist.csv
```

### 途中で見つけた 2 つの落とし穴

1. **公開範囲の違い**: optic-lobe では中枢脳のシナプスが公開されておらず、「シナプス 0」が見かけ上の 0 になる。
   そこで「もし結合があれば見えるはずか」（上流の出力・下流の入力がそのデータセットで見えているか）を確かめ、
   見えないものは照合不能として total に加算しない（`--min-syn-annotated`）。
2. **細胞体の左右と樹状突起の左右が逆**: FlyWire では KCg-m（左）→ MBON01（右、領域は MB_ML_L）のように、
   細胞体が樹状突起と反対側にある細胞がある。細胞体の左右で揃えて数えると本物の結合を見落とす
   （KCg-m → MBON01 が 0 に見えた）。判定は**両半球をまとめた合計**で行うことにした。修正後は 6,974 個。

### 結果

候補 375 組のうち、**3 個体すべてで照合できたのは 26 組**、そのうち別個体でも接触が再現したのは 17 組。

| 上流 → 下流 | 領域 | hemibrain 接触 | 3 個体の合計シナプス | male-cns 接触 | 上流ドライバー | 下流ドライバー |
|---|---|---|---|---|---|---|
| **VL1_vPN → LHAV1a3** | LH | 127 µm | **2** | 76–104 µm | **SS86025**（特異性 2） | **LH1543**（Ideal） |
| **ExR3 → FB8F_a** | FB | 260 µm | **2** | 31–41 µm | SS82798・SS82773 | SS75533・SS69168（_a と _b 両方の可能性） |
| **ExR3 → FB5AA** | FB | 217 µm | **2** | 50–87 µm | SS82798・SS82773 | SS71548 |
| **M_vPNml87 → LHAV2m1** | LH | 96 µm | **0** | 41–48 µm | SS00587・SS00586（特異性 1） | なし |
| **M_vPNml86 → LHPD3a5** | LH | 93 µm | **0** | 16–47 µm | なし | なし |
| M_l2PNl22 → LHAD2b1 | LH | 120 µm | 1 | 30–75 µm | SS32193（特異性 2） | なし |
| M_vPNml79 → LHAV3g1 | LH | 89 µm | 1 | 50–53 µm | なし | なし |
| ExR6 → LAL015 | EB・LAL | 61 µm | 1 | 33–40 µm | SS53617・SS58754 | なし |
| VA5_lPN → LHPV4a11 | LH | 127 µm | 2 | 21–31 µm | なし | なし |

- **3 個体すべてでシナプスが 1 個も無く、かつ接触が再現したのは 2 組だけ**（M_vPNml87 → LHAV2m1、M_vPNml86 → LHPD3a5）。
- **ドライバーが両側そろう中では VL1_vPN → LHAV1a3 が最有力**。上流は 1 細胞だけを標識する系統、下流は Ideal 評価。
  接触の再現も一番良い（hemibrain 127 µm に対し male-cns で 76–104 µm）。
- 3D（male-cns）: `results/cross_3d/`。

### 個体差が見えた組

| 組 | hemibrain | FlyWire | male-cns | 解釈 |
|---|---|---|---|---|
| **PAM07 → MBON01** | 0 | 0 | 順 7・逆 4 | 3 個体のうち male-cns だけ結合。個体差か雌雄差（male-cns のみ雄）。強制形成には不適 |
| FB6L → PFL2・PFL3 | 0 | 0 | 0 | シナプスは 3 個体とも 0 だが、**male-cns では接触がほぼ消える**（0–14 µm）。hemibrain の接触が偶然の可能性 |
| KCab-p → MBON13 | 0 | 0 | 0 | 同じく接触が再現しない（1–4 µm） |

- **注意**: Schlegel et al. 2024 によると hemibrain のタイプの約 3 分の 1 は FlyWire で確実に同定できない。
  照合できなかった組（LH で 33 組、キノコ体で 14 組など）は「結合が無い」ではなく「分からない」。
- **注意**: FlyWire の注釈で `hemibrain_type` が "AVLP386,AVLP388" のように複数あるときは、どちらにも数えている
  （多めに数える方が「シナプスが無い」と言い切るうえで安全なため）。

## 以前の 2 個体での確認（`validate_cross.py` → `shortlist.py`）

人工分子でシナプスを作らせる実験に向けて、判定基準を変えた。期待値との差ではなく、
**「両方向・両個体でシナプスがほぼ 0」かつ「膜が接している」** ことを条件にする
（接着分子が働くには向かい合っている必要があり、新しくできたシナプスを 0 を基準に数えたい）。

別個体は neuPrint の **male-cns:v1.0**（左右の半球がある。タイプ名は hemibrain と同じ）。

```bash
# 段階 1: 結合だけ確認（スケルトン不要。速い）
.venv/Scripts/python validate_cross.py --data data/hemibrain_cx.pkl --top 150
# 段階 2: 残った組の接触を確認（スケルトンを取得）
.venv/Scripts/python validate_cross.py --data data/hemibrain_cx.pkl --stage contact --top 10 --max-total 2
# 4 領域をまとめて最終候補に
.venv/Scripts/python shortlist.py --max-total 2
```

確認した 375 組のうち、148 組は male-cns にタイプ名が無く照合できなかった。
合計シナプス 2 個以下は 70 組。接触まで確認した 20 組が `results/shortlist.csv`・`shortlist.md`。

### 最終候補（別個体でも接触が再現し、シナプスがほぼ 0）

| 上流 → 下流 | 領域 | hemibrain 接触 | シナプス合計 | male-cns 接触 R / L | 上流 GAL4 | 下流 GAL4 |
|---|---|---|---|---|---|---|
| **ExR3 → FB8F_a** | FB | 260 µm | **1** | 41 / 31 µm | SS82798・SS82773 | SS75533・SS69168（FB8F_a と _b 両方の可能性） |
| **ExR3 → hDeltaD** | FB | 168 µm | **1** | 26 / 22 µm | SS82798・SS82773 | SS65202・SS00045 |
| **PFGs → PFL1** | FB・LAL | 385 µm | 2 | 56 / 90 µm | SS52547・SS77032 | SS100561・SS47735（PFLC も標識） |
| ExR3 → FB5AA | FB | 217 µm | 2 | 87 / 50 µm | SS82798・SS82773 | SS71548 |
| ExR5 → FB5K | FB | 174 µm | 2 | 29 / 64 µm | SS51128・SS82294 | SS75848・SS70334 |
| ER6 → PFGs | EB・FB | 242 µm | 2 | 11 / 21 µm | SS58833・SS58778 | SS52547・SS77032 |
| **VA5_lPN → LHPV4a11** | LH | 127 µm | **0** | 31 / 21 µm | なし | なし |
| M_l2PNm16 → LHCENT8 | LH | 137 µm | 1 | **116 / 124 µm** | なし | なし |
| VL1_vPN → LHAV2b3 | LH | 183 µm | 2 | 53 / 140 µm | SS86025（特異性 2） | なし |
| V_l2PN → LHPD4c1 | LH | 137 µm | 2 | 102 / 107 µm | SS01336（特異性 1） | なし |

- **ドライバーがそろう中では ExR3 → FB8F_a が最有力**。両個体・両半球で接触があり、シナプスは hemibrain の逆向き 1 個だけ。
- **接触の再現が一番良いのは M_l2PNm16 → LHCENT8**（別個体でも 116・124 µm）。ただし両側ともドライバーが無い。
- **VA5_lPN → LHPV4a11 は 4 回の計測すべてで 0 個**だが、ドライバーが無い。
- 3D（male-cns）: `results/cross_3d/`。

### 落ちた候補（重要）

| 組 | 理由 |
|---|---|
| **PAM07 → MBON01** | hemibrain では両方向 0 だが、male-cns の右半球で順 7・逆 4 個。**個体によって結合する**ので、強制シナプス形成には不適 |
| FB6L → PFL2 / PFL3 | シナプスは両個体で 0 だが、**別個体では接触がほぼ消える**（0〜14 µm）。接触が偶然だった可能性 |
| ER6 → PFR_a、WEDPN4 → mALB1 | 同じく別個体で接触なし |
| ExR7 → PEN_a・hDeltaK、ExR1 → PEG、PFNv → PFL2 | male-cns では結合している |
| PAM15_b → MBON04 など 148 組 | male-cns にそのタイプ名が無く、照合できなかった |

- **注意**: CX のドライバー表にある記号（†・§・角かっこ）の意味は、論文の図の説明にも補足ファイルにも記載がなく未確認。
  補足ファイル 1 は研究室の作業メモに近い内容だった。系統を使う前に著者に確認が必要。
- **注意**: どの候補も LexA が無い。上流と下流に別の分子を入れるには、片側を LexA に変換する必要がある。

## 実験への落とし込み（細胞名とドライバー）

- 結果サイトの各ペアに、細胞タイプ名（instance があればそれも）と、`drivers.csv` に登録したドライバーを表示する。
  未登録のタイプには、調べるためのリンクを出す。`drivers.csv` に行を足せばサイトに反映される（`status` は「確認済み」か「要確認」）。
- 現在の登録内容: MBON01 → MB011B（確認済み。ただし MBON03・MBON04 も含む）、SS01308・MB210B（要確認）、KCg-m → MB009B（要確認）。
- **注意: ドライバーは細胞タイプ単位。** hemibrain の KC で見つかった一致ペアは KCg-m どうしなので、
  そのペアだけを遺伝学的に操作することはできない。実験につなげるには、**異なる細胞タイプのペア**
  （例: タイプ X の軸索がタイプ Y の樹状突起に広く接するのに結合しない）を探すのが現実的。
  そのためには、いろいろなタイプを含むデータを取得し、ペアの結果をタイプの組ごとにまとめる。

## 次に試すアイデア

- **G: データセット間の再現性** — hemibrain と FlyWire の両方で同じタイプのペアが「近いのに結合しない」かを確認（再構築エラーとの区別）
- **C + E の組み合わせ** — C の候補を E の p 値で並べ替える 2 段階の方法
- **F の改良** — パッチの境界をスケルトンの枝の単位で決める、隣り合うパッチの密度の空間的な相関を考慮する
- 合成データの改良: 軸索/樹状突起が混ざった局所介在ニューロン、タイプの組ごとに結合を避けるケース

## 注意

- E には細胞タイプが必要。hemibrain は neuPrint の type を自動で使う。FlyWire は `fetch_data.py flywire --types-csv` で
  Codex の細胞タイプ表を渡す。
- F には個々のシナプスの位置と相手が必要（ローダーが `Dataset.synapses` に入れる）。
- neuPrint のローダーは hemibrain で動作確認済み。FlyWire のローダーは**まだ実データで動作確認していない**（fafbseg の返り値の列名はバージョンで変わるので、初回は中身を確認すること）。
- 全ペア比較なので、まず ROI や細胞タイプで数百〜数千ニューロンに絞ってから使う。

## ドライバー表の拡張（下行ニューロンと中枢神経系全体）

```bash
.venv/Scripts/python.exe match_drivers.py   --import-dn  <Namiki 2018 eLife 34272 Supplementary file 3 (.xlsx)>   --import-cns <Meissner 2024 eLife 98405 Figure 1-source data 1 (.xlsx)>
```

`data/driver_lines_dn.csv`（134 系統 / 62 型）と `data/driver_lines_cns.csv`（2,454 系統 / 1,926 型）
を作る。これで照合できる型は 4 コレクションの頃の約 500 から **2,266 型**になった。

- Namiki の表は `bg`（背景の標識）を反転して特異性に読み替える。
- Meissner の表は中枢神経系全体。`Quality` は 1 が最良なので `3 - Quality` を特異性とし、
  `VNC only` の系統（脳に発現しない）と幼虫の系統は落とす。

LexA はこれだけ集めても **20 系統 / 12 型**しかなく、2 つの発現系がそろうかどうかが
候補選びの最大の制約になっている。

## 領域ごとの絞り込みをまとめて実行

```bash
./run_regions.sh smp          # 取得 → 接触の計算 → タイプの組 → ドライバー照合
./run_screen.sh  smp          # 順位付け → FlyWire → male-CNS → 接触領域
```

`run_screen.sh` は領域名を引数で受ける。`sed` でスクリプトを書き換えて使い回すと
`results/hemibrain_cre/` のような語中の領域名が置換されずに事故るので、必ずこちらを使う。

調べた領域（hemibrain）:

| 鍵 | 領域 | 細胞数 |
|---|---|---|
| `mb_all` | キノコ体 | — |
| `lh` | 外側角 | — |
| `cx` | 中心複合体 | — |
| `vis` | AOTU・PVLP | — |
| `al` | 触角葉 | — |
| `smp` | SMP・SIP・SLP | 1094 |
| `lal` | LAL・VES・GOR | 807 |
| `wed` | WED・SAD・CAN | 378 |
| `ps` | SPS・IPS・EPA・GNG | 548 |
| `cre` | CRE・ATL・IB・SCL・ICL | 717 |
| `vlp` | AVLP・PLP | 875 |

## LexA を軸にした再探索

2 つの発現系（上流と下流に別々の分子を入れる）がそろうかどうかが最大の制約。
集めた 4,045 系統のうち LexA は **20 系統・12 型**しかない:

```
MBON01 MBON03 MBON04 MBON07 MBON11 MBON12 PAM11 PPL101 PPL106 SLP400 SMP108 SMP348
```

そこで、この 12 型を含む組だけを全領域の候補から集めて絞り直した。
領域をまたぐので `--pairs-as-source` と `--tag` を使う:

```bash
.venv/Scripts/python.exe flywire_check.py --pairs results/lexa_zero_hemibrain.csv --tag lexa   --min-expected 3 --data data/hemibrain_smp.pkl data/hemibrain_cre.pkl data/hemibrain_mb_all.pkl
.venv/Scripts/python.exe validate_cross.py --data data/hemibrain_smp.pkl --tag lexa   --pairs results/lexa/zero2_cand.csv --pairs-as-source --min-expected 3 --datasets male-cns:v1.0
```

結果は 245 組 → 98 組（hemibrain でゼロ）→ 16 組（FlyWire でもゼロ）→ **7 組**（male-CNS でもゼロ）。
7 組すべてが別個体の両半球で接触を再現した。ただしドライバーがそろうのは 2 組だけで、
残り 5 組は下流（CRE022・CRE050・SIP067・LAL175・SMP053）の系統が 1 つも無い。

## 本命候補の詳細

`results/candidate_MBON04_FB4N.md` に MBON04 → FB4N の調査をまとめた。
要点は「同じ β'2 区画の MBON26・MBON27 は FB4N と結合し、MBON04 自身も同じ 4 層の
FB4C と結合しているのに、この 1 組だけ 108.9 µm 接触して 0」という特異性と、
FB4N が既にグルタミン酸入力を 23.8% 受けている（受容体がある）こと。

## 接触部位に何が挟まっているかを EM で確かめる

```bash
.venv/Scripts/python.exe em_mesh_gap.py   # メッシュ（膜）どうしの最短距離
.venv/Scripts/python.exe em_gap_check.py --pairs results/mb_em_targets.csv  # 間にある物体を数える
```

`em_gap_check.py` は cloud-volume で hemibrain のセグメンテーション
（`gs://neuroglancer-janelia-flyem-hemibrain/v1.2/segmentation`、8 nm 等方）を
切り出し、A のボクセルと B のボクセルが直接隣接しているか、隣接していなければ
間にどの物体 ID があるかを数える。骨格ベースの「1 µm 以内」では分からない
「グリアや別のニューロンが挟まっていないか」をこれで判定できる。

neuroglancer のリンクは `results/em_links_MBON04_FB4N.md` に生成してある。

## キノコ体の中だけで強制シナプス形成を行う案

`results/plan_mushroom_body.md` にまとめた。LexA がある 12 型がすべてキノコ体の細胞なので、
探索をキノコ体に閉じると「2 つの発現系がそろう組」が全脳 2 組 → 12 組に増える。

絞り込みで効いたのは次の 2 つ。どちらも全脳の探索には無かった工程:

- **型単位の数え直し** — キノコ体のデータは細胞を間引いて取得しているので、
  間引いた細胞でゼロでも型全体ではゼロでない組が混ざる（43 組中 6 組）。
- **ドライバーの混同の判定** — 上流のドライバーが標識する型の集合 U と下流の集合 D の、
  **U × D すべての組**が 0 でなければならない。U どうし・D どうしの既存結合は
  前側と後側の分子が揃わないので読み出しに混ざらない。22 組中、本当に使えないのは 3 組だけ。

残った 4 組のうち上位 2 組を詰めた:

| 候補 | 報告書 | 強み |
|---|---|---|
| **MBON04 ↔ PPL105** | `results/candidate_MBON04_PPL105.md` | 接触 93 µm、別個体で 33〜38 µm（全体で最高）、両方向が候補になる、同じクラスターの PPL107 が対照になる |
| **MBON07 → PPL106** | `results/refs_top_pairs.md` 2 節 | ドライバーが両方とも単一タイプ・特異性 2（MB319B + SL86062）。U × D が 1 組だけ |

EM と 3D の neuroglancer リンクは `results/em_links_MBON04_PPL105.md`（`em_links.py` で再生成できる）。

## LexA ドライバーの棚卸し（キノコ体以外・神経とグリア）

`results/lexa_catalog.md` と `data/lexa_lines_nonmb.csv`。

「LexA は 20 系統・12 型しかない」は**細胞型特異な split-LexA に限った話**で、
範囲を広げると GMR-LexA 約 1,400 系統・VT-LexA・遺伝子ノックイン LexA がある。
グリアでは repo-LexA（全グリア）・alrm-LexA（アストロサイト）・wrapper-LexA（コルテックス）・
R83E12-LexA（BDSC 54288、エンシーシング）・R85G01-LexA（ペリニューリアル）・
R86E01-LexA（視葉アストロサイト）が特定できた。
Shuai 2024 が split-GAL4 → split-LexA の変換を 34 系統で試して 22 系統成功しているので、
**このボトルネックは原理的な制約ではない**。

## GMR-LexA の照合

`results/gmr_lexa_report.md`（`python gmr_lexa_match.py` で再現）。

「上流にも下流にも split-GAL4 はあるが LexA が無い」ために落ちた 32 組を、
GMR-LexA（Rubin のエンハンサー断片 LexA）で救えるかを調べた。**結論は救えない。**

- GMR-LexA の実在を確認できたのは **150 系統**（FlyBase の allele `Ecol\lexA[GMRxxxxx]`、VFB 経由）。
- そのうち標識細胞が分かるのは **286 の発現パターン / 43 件が hemibrain の型に対応**。
- 落ちた 32 組の細胞型 49 種類のうち、注釈があるのは **0 件**。
  LAL098・DNb05・SMP503 のようなコネクトーム由来の型名は、Gen1 の注釈体系に載っていない。
- キノコ体では代替が見つかった（MBON04 に **R15B01-LexA** = SL86063 と同じ構成、
  PAM11 に R71C03 / R72D01）が、どれも既に split-LexA がある細胞。

次の一手は **split-GAL4 → split-LexA の変換**（Shuai 2024 で 22/34 成功）か、
split-GAL4 のヘミドライバー断片に LexA 版があるかの確認。

## ヘミドライバー断片の LexA 版

`results/hemidriver_lexa_report.md`（`python hemidriver_lexa.py` で再現）。

split-GAL4 は 2 つのエンハンサー断片の積（例 `SS49882 = VT045661-GAL4.DBD ∩ VT039475-p65.AD`）。
**片方の断片に LexA 版があれば、特異性は落ちるが発現系は 2 つになる。**
VirtualFlyBrain の neo4j から 60 系統のヘミドライバー構成を取り、
FlyBase の lexA allele 190 件と突き合わせた。

| | |
|---|---|
| 片方の断片に LexA 版がある split-GAL4 | **7 系統** |
| それで救える組 | **6 組**（上位の条件を満たすものが **5 組**） |

救えた 5 組（★は接触が別個体でも再現し、上流もその領域で出力している）:
MBON30 → oviIN（VT026347-LexA）、MBON21 → FB4X（GMR11C07-LexA）、
ExR2 → LAL098（VT045661-LexA）、AVLP478 → DNp06（VT019018-LexA）、
DL4_adPN → LHPV6c1（VT033006-LexA）。

**ExR2 → LAL098 は別個体での接触 47〜77 µm で、これまでのどの候補より良い。**
ただし断片 1 つが標識する細胞の全体（U / D の中身）は未確定なので、
混同の判定は FlyLight の画像か NeuronBridge の逆引きで取り直す必要がある。

## 最終候補のページ（results/shortlist.html）

`python shortlist.py --max-total 0 && python build_docs.py && python build_site.py` で作る。
3 個体すべてで両方向のシナプスが 0 だった 306 組を 1 つの表にして、判定・接触・膜の密着・
ドライバー・3D・報告書へのリンクを並べる。列はクリックで並べ替えられる。

| 列 | 中身 |
|---|---|
| 判定 | ◎ 発現系が 2 つそろう／○ 両方 GAL4／△ 片側のみ／要注意 上流がそこで出力していない |
| 接触 元 → 別個体 | hemibrain の接触と、別個体の左右それぞれでの接触 |
| 上流出力 / 下流入力 | 主な接触領域にある上流のシナプス前部位と下流の後部位の数 |
| **膜の密着** | EM のセグメンテーションで数えた直接隣接ボクセルと、間に何も無い表面点の割合（`em_gap_check.py`） |
| 3D・報告書 | 骨格の 3D（hemibrain と別個体）、その組について書いた報告書と EM リンク集 |

3D のファイル名は作ったときの向きで決まるので、`t1_t2` と `t2_t1` の両方を探してリンクする。
別個体の 3D は `python cross_3d.py --data data/male_cns_top.pkl --pair MBON04 PPL105 …` で作る。

**注意**: `shortlist.py` は各領域の `results/` を読み直して `shortlist_final.csv` を作る。
新しく接触や接触領域を測ったら、**必ず `shortlist.py` を流し直す**こと。
流し忘れると、測ったはずの組が「接触未確認」のまま表に出る。

## サイトと公開

- `python build_docs.py` — `results/*.md` を HTML に変換（`results/docs/`）。
  表・コードブロック・報告書どうしのリンクをそのまま出す。`build_site.py` の入口にも出る。
- `python publish_site.py` — GitHub Pages に載せる分だけを **`docs/`** に集める（約 103 MB）。
  GitHub Pages の「Deploy from a branch」はフォルダを `/ (root)` か `/docs` しか選べないので、
  出力先は `docs/` に固定している（報告書は `docs/reports/`）。
  `data/`（5.4 GB）とデータセット別の全ペア 3D（3.5 GB）は `.gitignore` で除外。
  公開設定は Settings → Pages → Branch: `main` / フォルダ: `/docs`。

## 上位候補の行動実験の文献

`results/refs_top_pairs.md`。`lit_search.py` で Europe PMC の本文を検索した結果から、
候補 4 組に出てくる 9 つの細胞について「操作 + 行動」の前例を整理した。
MBON04 が 10 本で最も恵まれており、MBON31 と FB4N は 0 本。
