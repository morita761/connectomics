# キノコ体以外の LexA ドライバー（神経とグリア）

2026-11-10 調査。構造化した一覧は `data/lexa_lines_nonmb.csv`（34 行）。

## 0. 結論を先に

これまで「LexA 系統は 20 本・12 タイプしかなく、すべてキノコ体」としていたが、
これは **「1 つの細胞型だけを狙える split-LexA」に限った話**だった。範囲を広げると状況が変わる。

| 種類 | 規模 | 特異性 | この実験に使えるか |
|---|---|---|---|
| **split-LexA**（細胞型特異） | 20 本 / 12 タイプ | ◎ 1 細胞型 | **そのまま使える**。ただし全部キノコ体 |
| **GMR-LexA**（Rubin, エンハンサー断片） | **約 1,400 本** | △ 断片ごとに数個〜数十の細胞型 | 混同（U × D）の確認が必須 |
| **VT-LexA**（Vienna Tiles） | 数千本のうち LexA 版 | △ 同上 | 同上 |
| **遺伝子ノックイン LexA** | 数十遺伝子 | ✕ クラス単位（全グリア・全ドーパミンなど） | 単独では広すぎる |
| **StanEx エンハンサートラップ** | 約 100 本＋後続数百 | △ 脳以外が中心 | 脳の候補には向かない |

**キノコ体の外に出るなら、「細胞型特異な LexA が無い」のではなく
「細胞型特異性が落ちる LexA しか無い」というのが正確な状況。**
GMR-LexA 約 1,400 本は、断片が標識する細胞型の集合 U（または D）を FlyLight の画像で確定できれば、
U × D がすべて 0 かどうかを今のパイプラインで判定できる。これは次の探索の現実的な出口になる。

## 1. グリアの LexA

グリアはコネクトームに登録されていないので「シナプス 0 の相手」としては使えないが、
**接触部位にグリアを介入させる／除く操作の道具**としては使える。

| 系統 | 構築 | 標的 | 備考 | 出典 |
|---|---|---|---|---|
| **repo-LexA** | P{repo-LexA-GAD}, M{repo-LG.B}, M{repo-LHG.B} | **全グリア** | repo は全グリアのマーカー（中央線グリアを除く）。FlyBase に pan-glial と明記 | FlyBase [FBal0221468](https://flybase.org/reports/FBal0221468.html) |
| **repo-T2A-LexA-GAD** | ノックイン | 全グリア | GAL4 版で出る腸・肛門板の余計な発現が LexA 版では出ない | Zirin 2024 [eLife 94073](https://elifesciences.org/reviewed-preprints/94073v1) |
| **alrm-LexA** | P{alrm(-FRT)lexA::GAD} (FBtp0098149), P{alrm-LexA::GAD.S}, M{alrm-LG.B}, M{alrm-LHG.B} | **アストロサイト様グリア** | alrm はグリアの中でアストロサイトだけを標識する（repo-GAL80 で抑制される） | FlyBase [FBtp0098149](https://flybase.org/reports/FBtp0098149), [FBto0000153](https://flybase.org/reports/FBto0000153) |
| alrm(FRT.GAL4)lexA::GAD | P{alrm(FRT.GAL4)lexA::GAD} (FBtp0098148) | アストロサイト | FRT で挟んだ GAL4 を FLP で切り出すと LexA になる。モザイク解析向き | FlyBase [FBtp0098148](https://flybase.org/reports/FBtp0098148) |
| **wrapper-LexA** | P{wrapper-lexA::GAD.932i} | コルテックスグリア・中央線グリア | wrapper は細胞体を包むグリアのマーカー | FlyBase [FBto0000153](https://flybase.org/reports/FBto0000153) |
| **R86E01-LexA** | P{R86E01-lexA::GADfl} | アストロサイト（視葉）、ラミナの上皮グリア | エンハンサー断片。GMR コレクション由来 | glia atlas 2023 [PLOS Biol](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002328) |
| **R85G01-LexA** | P{R85G01-lexA::GADfl} (attP40) | ペリニューリアルグリア | メデュラ・ロブラ・ロブラプレートで形態が均一 | 同上 |
| **R83E12-LexA** | BDSC **54288** | エンシーシンググリア | 「ensheathing glia に特異的に見える」と記載。**在庫番号が判明している唯一のグリア LexA** | Nat Commun 2021 [βH-Spectrin](https://www.nature.com/articles/s41467-021-26462-x) |

**注意点**: alrm-LexA は FlyBase 上で紐づく系統が 0 なので、構築はあるが
配布系統が無い可能性がある。実際に取り寄せられるのは R83E12-LexA (BDSC 54288) と repo-LexA が確実。
アストロサイト LexA が必要なら、2026 年の Cell Reports Methods の論文が
「LexA>LexAop では効率が足りず split-GAL4>UAS に置き換えた」と報告しているので、
効率の点でも事前検討が必要。

## 2. 神経の LexA（キノコ体以外）

### 2-1. 細胞型・サブタイプが絞れるもの

| 系統 | 標的 | 細胞数 | 出典 |
|---|---|---|---|
| **Or82a-LexA** | 単一の嗅覚受容ニューロン型（糸球体 1 つ） | 数十 | FlyBase [FBto0000153](https://flybase.org/reports/FBto0000153) |
| **Rh3 / Rh4 / Rh5 / Rh6-LexA** | 光受容細胞 R7/R8 のサブタイプ | — | 同上 |
| **SIFa-LexA** | SIFamide ニューロン | **4 個** | 同上 |
| **Ilp2-LexA** | インスリン産生細胞（IPC） | 約 14 個 | 同上, Zirin 2024 |
| **GH146-LexA** | 嗅覚投射ニューロン | 約 90 個 | FlyBase [FBal0244422](https://flybase.org/reports/FBal0244422.html), BDSC 606260 |
| **eve(CQ2)-LexA** | 腹髄の CQ ニューロン | — | FlyBase FBto0000153 |
| **Gr43a / Gr66a split-LexA** | 糖受容 / 苦味受容ニューロン | — | Shuai 2024 [eLife 94168](https://elifesciences.org/articles/94168) |

**Gr43a / Gr66a split-LexA は、キノコ体以外で唯一の split-LexA**。
Shuai 2024 は split-GAL4 を split-LexA に変換する試みを 34 系統で行い 22 系統が成功しており、
**同じ変換を他の split-GAL4 に適用できる**ことを示している。これが最も有望な拡張手段。

### 2-2. クラス単位（単独では広すぎる）

| 系統 | 標的 |
|---|---|
| **elav-LexA** | 全神経 |
| **ple-LexA** | ドーパミンニューロン全般（ple = pale = チロシン水酸化酵素）。**PPL105 を含む** |
| **ChAT-LexA** | コリン作動性ニューロン全般 |
| **Orco-LexA** (= Or83b-LexA) | 嗅覚受容ニューロン全般 |
| **cry-LexA** | 概日時計ニューロン |
| **NPF-LexA** | NPF 作動性ニューロン |
| **ppk-LexA** | 多樹状突起性感覚ニューロン（末梢） |
| **tubP-LexA / Act5C-LexA** | ユビキタス（対照・効果器の試験用） |

Chemoconnectome（Deng 2019 Neuron）は伝達物質・神経ペプチド・受容体の遺伝子を
網羅的にノックインしており、**10 遺伝子で LexA 版**がある。GAL4 版とほぼ同じパターンになる。

## 3. 大規模コレクション

| コレクション | 規模 | 活性化因子 | 挿入部位 | 出典 |
|---|---|---|---|---|
| **GMR-LexA**（Rubin） | **約 1,400 本** | lexA::GAD / lexA-p65 | ほぼ attP40 (25C7, 2L) | Pfeiffer 2010 Genetics 186:735; [BDSC Janelia lexA](https://bdsc.indiana.edu/stocks/lexa/lexa_janelia.html) |
| **VT-LexA**（Vienna Tiles） | VT 18,000 本のうち LexA 版 | lexA::GAD | — | Tirian & Dickson 2017 [bioRxiv 198648](https://www.biorxiv.org/content/10.1101/198648v1.full) |
| **StanEx LexA** エンハンサートラップ | 約 100 本（2016）＋数百本（2019・2023） | LHG | P 因子の転移、挿入をマッピング済み | Kockel 2016 G3 6:3017 [PMC5068927](https://pmc.ncbi.nlm.nih.gov/articles/PMC5068927); Kockel 2019 G3 9:2097 |
| **T2A-LexA-GAD ノックイン** | da, elav, **repo**, nub, btl, hh, dpp, CG9458, CG17560 ほか | LexA-GAD | 内在遺伝子の 5' コード領域直後 | Zirin 2024 [eLife 94073](https://elifesciences.org/reviewed-preprints/94073v1) |
| **Chemoconnectome (CCT)** | 伝達物質・ペプチド・受容体。10 遺伝子で LexA 版 | LexA | ノックイン | Deng 2019 [Neuron 101:876](https://www.cell.com/neuron/fulltext/S0896-6273(19)30072-8) |

**GMR-LexA の設計上の利点**: LexA は attP40（2L）、GAL4 は attP2（3L）に入っているので、
同じ個体で組み合わせたときに transvection（染色体間の干渉）が起きにくい。
GAL4 + LexA を 1 匹に入れる前提で作られたコレクションである。

**GMR/VT の限界**: どちらも脳のエンハンサー断片なので、非神経の標的には向かない。
また 1 断片が数個〜数十の細胞型を標識するため、本実験の
「U × D のすべての組が 0」という条件を満たすかは系統ごとに画像で確認するしかない。

## 4. 発現系を 3 つにする選択肢

LexA が足りないという制約そのものを外す手もある。

| 手段 | 内容 | 出典 |
|---|---|---|
| **UAS-LexA-GAD** | GAL4 の下流で LexA を発現させる。既存の split-GAL4 を「LexAop を駆動する系統」に変えられる。ただし GAL4 はその細胞に占有されるので、もう一方の細胞には第 3 の系が必要 | FlyBase FBto0000153 |
| **T2A-QF2-T2A-LexA-GAD 二重ドライバー** | 1 つのノックインで QF2 と LexA-GAD を同時発現。GAL4 を空けたまま 2 系統を使える | Zirin 2024 eLife 94073（hh・dpp で実証） |
| **split-GAL4 → split-LexA 変換** | DNA 結合ドメインを GAL4 から LexA に差し替える。34 系統試して 22 系統成功 | Shuai 2024 eLife 94168 |
| **split-QF2** | GAL4 / LexA に加えて QF/QUAS を使う | Zirin 2024 ほか |

**この実験に最も効きそうなのは「split-GAL4 → split-LexA 変換」**。
成功率 22/34（65%）という実績があり、既にある細胞型特異な split-GAL4
（本プロジェクトで集めた 4,045 系統 / 2,266 タイプ）のどれでも候補にできる。
つまり **LexA のボトルネックは「作れば解ける」種類の制約**であって、
現存系統だけで考えると候補をキノコ体に閉じ込めてしまう。

## 5. 次にやること（提案）

1. GMR-LexA 約 1,400 本について、FlyLight の発現画像から標識細胞型の集合を取り、
   既存の候補ペアの下流側に当てはまるものを探す。
2. 全脳探索で「両方に split-GAL4 があるが LexA が無い」ために落ちた 14 組
   （`results/shortlist_final.csv`）を、split-LexA への変換を前提に再評価する。
3. グリアを介入させる実験を組むなら、まず R83E12-LexA (BDSC 54288, エンシーシンググリア) と
   repo-LexA を確保する。接触部位（SIP など）にどのグリアが来るかは
   EM の切り出しで確認できる（`em_gap_check.py` の「間にあった物体」）。

## 参考文献

- Pfeiffer BD, Ngo TT, Hibbard KL, Murphy C, Jenett A, Truman JW, Rubin GM (2010)
  Refinement of tools for targeted gene expression in Drosophila. *Genetics* 186:735–755.
  doi:10.1534/genetics.110.119917 — GMR-LexA コレクションの元論文。LexA は GAL4 より活性が弱いと報告。
- Tirian L, Dickson BJ (2017) The VT GAL4, LexA, and split-GAL4 driver line collections
  for targeted expression in the Drosophila nervous system. *bioRxiv* 198648.
  doi:10.1101/198648
- Kockel L, et al. (2016) A Drosophila LexA enhancer-trap resource for developmental biology
  and neuroendocrine research. *G3* 6:3017–3026. doi:10.1534/g3.116.031229
- Kockel L, et al. (2019) An interscholastic network to generate LexA enhancer trap lines
  in Drosophila. *G3* 9:2097–2106. doi:10.1534/g3.119.400105
- Kockel L, et al. (2023) Generation of LexA enhancer-trap lines in Drosophila by an
  international scholastic network. *G3*. PMC10468311
- Zirin J, Jusiak B, Lopes R, Ewen-Campen B, Bosch JA, Risbeck A, Forman C, Villalta C,
  Hu Y, Perrimon N (2024) Expanding the Drosophila toolkit for dual control of gene
  expression. *eLife* 94073. doi:10.7554/eLife.94073 — repo を含む T2A-LexA-GAD ノックイン。
- Deng B, et al. (2019) Chemoconnectomics: mapping chemical transmission in Drosophila.
  *Neuron* 101:876–893. doi:10.1016/j.neuron.2019.01.045
- Shuai Y, et al. (2024) Driver lines for studying associative learning in Drosophila.
  *eLife* 94168. doi:10.7554/eLife.94168 — split-LexA 20 本（キノコ体）＋ Gr43a/Gr66a。
  split-GAL4 → split-LexA 変換の成功率 22/34。
- Kremer MC, Jung C, Batelli S, Rubin GM, Gaul U (2017) The glia of the adult Drosophila
  nervous system. *Glia* 65:606–638. doi:10.1002/glia.23115
- Drosophila glial cell atlas (2023) *PLOS Biology* 21:e3002328.
  doi:10.1371/journal.pbio.3002328 — R86E01-LexA（アストロサイト）、R85G01-LexA（ペリニューリアル）。
- Nature Communications (2021) Drosophila βHeavy-Spectrin is required in polarized
  ensheathing glia. doi:10.1038/s41467-021-26462-x — R83E12-LexA (BDSC 54288)。
- FlyBase lexA::GAD 実験ツール一覧: https://flybase.org/reports/FBto0000153
  （構築 200 件超、挿入 619 件。LexA 系統を網羅的に探すならここが起点）
- BDSC Janelia lexA 在庫: https://bdsc.indiana.edu/stocks/lexa/lexa_janelia.html
  （対照系統 77691 = 制御配列なしの lexA-p65、attP40）
