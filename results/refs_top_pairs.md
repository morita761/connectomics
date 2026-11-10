# 上位候補ペアの細胞について、行動実験を行っている論文

2026-11-10 作成。`lit_search.py` で Europe PMC を検索し、
オープンアクセス論文は**本文全体**、それ以外は要旨を対象に
「細胞名・系統名」「操作（光遺伝学・活動阻害など）」「行動（学習・選好など）」の
3 つがそろう文を探した結果。生データは `results/literature/<タイプ>.csv`、
根拠の文は `results/literature/summary.md`。

```
python lit_search.py PPL105 PPL106 MBON07 PAM11 MBON33 MBON04 MBON11 MBON31 FB4N --max-papers 200
```

**注意**: 語の一致で拾っているので、最終的な判断は論文本体で確認すること。

## 0. 細胞ごとの収穫

| 細胞 | 候補ペアでの役割 | 検索された論文 | 行動実験の根拠がある論文 |
|---|---|---|---|
| **MBON04** (MBON-β'2mp_bilateral) | 第一候補の上流／下流 | 64 本 | **10 本** |
| **PPL105** (PPL1-α'2α2) | 第一候補の下流／上流 | 54 本 | **4 本** |
| **PPL106** (PPL1-α3) | 第2候補の下流 | 49 本 | **4 本** |
| **MBON07** (MBON-α1) | 第2候補の上流 | 36 本 | **5 本** |
| **PAM11** (PAM-α1) | 第3候補の下流 | 46 本 | **3 本** |
| MBON11 (MBON-γ1pedc) | 却下された組の上流 | 107 本 | 5 本 |
| MBON33 (MBON-γ2γ3) | 却下された組 | 6 本 | 1 本 |
| **MBON31** (MBON-α'1a) | 第3候補の上流 | 15 本 | **0 本** |
| **FB4N** | 全脳本命の下流 | 3 本 | **0 本** |

**読み出しの観点では MBON04 ↔ PPL105 が最も恵まれている。**
両方の細胞に行動実験の前例があり、同じ論文（Tsao 2018、Aso 2016）が両方を扱っている。

**MBON31 → PAM11 は上流側（MBON31）に行動実験の前例が無い。** 15 本のうち
行動の根拠がある論文は 0 本で、最も関連するのは Felsenberg 2017 Nature で 3 回言及されるだけ。
下流の PAM11 には Ichinose 2015 eLife という決定的な論文があるので、
読み出しは PAM11 側（報酬性長期記憶）に置くしかない。

**FB4N（全脳本命の下流）は 3 本しかなく、すべて形態の記載。** 行動の前例が無い点は変わらない。

---

## 1. MBON04 ↔ PPL105（第一候補）

### MBON04 = MBON-β'2mp_bilateral

| 年 | 著者 | 論文 | 誌 | DOI | 何をしたか |
|---|---|---|---|---|---|
| 2014 | Aso Y, et al. | Mushroom body output neurons encode valence and guide memory-based action selection in Drosophila | eLife | [10.7554/eLife.04580](https://doi.org/10.7554/eLife.04580) | **MB011B** で MBON-β'2mp_bilateral を含む 5 つのグルタミン酸作動性 MBON を CsChrimson で活性化し、忌避反応を定量。価数のマッピングの基礎論文 |
| 2015 | Owald D, Felsenberg J, Talbot CB, Das G, Perisse E, Huetteroth W, Waddell S | Activity of defined mushroom body output neurons underlies learned olfactory behavior in Drosophila | Neuron | [10.1016/j.neuron.2015.03.025](https://doi.org/10.1016/j.neuron.2015.03.025) | MBON-β'2mp を VT1211 / R21D02 で阻害し、報酬・忌避記憶の両方で成績低下を示した |
| 2015 | Sitaraman D, Aso Y, Rubin GM, Nitabach MN | Control of sleep by dopaminergic inputs to the Drosophila mushroom body | Front Neural Circuits | [10.3389/fncir.2015.00073](https://doi.org/10.3389/fncir.2015.00073) | MBON-γ5β'2a/β'2mp/β'2mp_bilateral が覚醒を促し、生来の忌避を媒介することを光遺伝学で示した |
| 2018 | Tsao CH, Chen CC, Lin CH, Yang HY, Lin S | Drosophila mushroom bodies integrate hunger and satiety signals to control innate food-seeking behavior | eLife | [10.7554/eLife.35264](https://doi.org/10.7554/eLife.35264) | **MBON04 と PPL105 の両方を扱う**。β'2 の MBON と PPL1-α'2α2 を操作して酵母の探索行動を測定 |
| 2018 | Al-Anzi B, Zinn K | Identification and characterization of mushroom body neurons that regulate fat storage in Drosophila | Neural Dev | [10.1186/s13064-018-0116-7](https://doi.org/10.1186/s13064-018-0116-7) | MBON-β'2mp を含む MBON の阻害・過活性化で脂肪蓄積と摂食行動を測定 |
| 2022 | Hobin M, Dorfman K, Adel M, Rivera-Rodriguez EJ, Kuklin EA, Ma D, Griffith LC | The Drosophila microRNA bantam regulates excitability in adult mushroom body output neurons to promote early night sleep | iScience | [10.1016/j.isci.2022.104874](https://doi.org/10.1016/j.isci.2022.104874) | **MB011B** で γ5β'2a/β'2mp/β'2mp_bilateral の興奮性を操作し睡眠を測定。根拠文 12 本で最多 |
| 2023 | Zhou F, et al. | Optimized design and in vivo application of optogenetically functionalized Drosophila dopamine receptors | Nat Commun | [10.1038/s41467-023-43970-0](https://doi.org/10.1038/s41467-023-43970-0) | **MB011B で optoDop1R2 を発現**させ、摂食行動（sip 数）を測定。**MBON04 にドーパミン受容体を導入した前例** |

**Zhou 2023 は特に重要**: MBON04 を含む MBON に光操作できるドーパミン受容体を導入し、
行動（摂食）で読み出した例である。PPL105 → MBON04 の向きで新しいドーパミンシナプスを作る本実験と、
操作の対象と読み出しがそのまま重なる。

### PPL105 = PPL1-α'2α2（= MB-V1）

| 年 | 著者 | 論文 | 誌 | DOI | 何をしたか |
|---|---|---|---|---|---|
| 2012 | Aso Y, Herb A, Ogueta M, Siwanowicz I, Templier T, Friedrich AB, Ito K, Scholz H | Three dopamine pathways induce aversive odor memories with different stability | PLoS Genet | [10.1371/journal.pgen.1002768](https://doi.org/10.1371/journal.pgen.1002768) | **MB-V1（= PPL1-α'2α2）を含む 3 経路を切り分け**、安定性の異なる嫌悪記憶を誘導することを示した。根拠文 9 本 |
| 2016 | Aso Y, Rubin GM | Dopaminergic neurons write and update memories with cell-type-specific rules | eLife | [10.7554/eLife.16135](https://doi.org/10.7554/eLife.16135) | **MB099C**（PPL103+104+105+106）で CsChrimson 駆動し、嫌悪記憶の形成を区画ごとに比較。**PPL105 と PPL106 の両方を扱う** |
| 2018 | Tsao CH, et al. | （上記） | eLife | [10.7554/eLife.35264](https://doi.org/10.7554/eLife.35264) | PPL1-α'2α2 の人工活性化が酵母匂いの探索を促進。Aso 2016 の嫌悪記憶とは逆向きの効果 |
| 2021 | Jacob PF, Vargas-Gutierrez P, Okray Z, Vietti-Michelina S, Felsenberg J, Waddell S | Prior experience conditionally inhibits the expression of new learning in Drosophila | Curr Biol | [10.1016/j.cub.2021.05.056](https://doi.org/10.1016/j.cub.2021.05.056) | PPL1-α'2α2 の阻害は pre-exposure 学習には不要と判定（除外の根拠として有用） |

### 関連（PPL105 の光遺伝学・イメージング）

- Toshima N, et al. Individual dopaminergic neurons induce unique, yet overlapping combinations
  of memory components. *eLife* — 個々の DAN が誘導する記憶成分の分解。
- Huang C, et al. (2018) Long-term optical brain imaging in live adult fruit flies.
  *Nat Commun* [10.1038/s41467-018-02873-1](https://doi.org/10.1038/s41467-018-02873-1)
  — PPL1 を長期イメージングした例。新生シナプスの経時観察に使える手法。

---

## 2. MBON07 → PPL106（第2候補）

### MBON07 = MBON-α1

| 年 | 著者 | 論文 | 誌 | DOI | 何をしたか |
|---|---|---|---|---|---|
| 2015 | Ichinose T, Aso Y, Yamagata N, Abe A, Rubin GM, Tanimoto H | Reward signal in a recurrent circuit drives appetitive long-term memory formation | eLife | [10.7554/eLife.10719](https://doi.org/10.7554/eLife.10719) | **MB323B で MBON-α1 を操作**。活性化を匂いと対提示して 24 時間記憶を測定、試験時のみの阻害で条件付け接近が大きく損なわれる。**MBON07 と PAM11 の両方を扱う決定的な論文**。根拠文 8 本 |
| 2018 | Tsao CH, et al. | （上記） | eLife | [10.7554/eLife.35264](https://doi.org/10.7554/eLife.35264) | **MB310C**（MBON-α1 特異的）で GABA-A 受容体をノックダウンし、酵母探索行動の減少を示した |
| 2023 | Mohamed A, Malekou I, Sim T, O'Kane CJ, Maait Y, Scullion B, Masuda-Nakagawa LM | Mushroom body output neurons MBON-a1/a2 define an odor intensity channel that regulates behavioral odor discrimination learning in larval Drosophila | Front Physiol | [10.3389/fphys.2023.1111244](https://doi.org/10.3389/fphys.2023.1111244) | 幼虫の MBON-a1/a2 を CsChrimson で活性化し、匂い識別学習が損なわれることを示した。言及 117 回で最多 |
| 2024 | Davidson AM, Hige T | Roles of feedback and feed-forward networks of dopamine subsystems | Learn Mem | [10.1101/lm.053807.123](https://doi.org/10.1101/lm.053807.123) | MBON-α1 の抑圧が SMP353/354 → SMP108 → PAM を介して二次条件づけに関わる回路を整理した総説 |
| 2024 | Nanami T, Yamada D, Someya M, Hige T, Kazama H, Kohno T | A lightweight data-driven spiking neuronal network model of Drosophila olfactory nervous system | Front Neurosci | [10.3389/fnins.2024.1384336](https://doi.org/10.3389/fnins.2024.1384336) | MBON-α1 の活性化は上風走行を誘導しない（MBON-α3 とは異なる）という対照的な知見 |

### PPL106 = PPL1-α3

| 年 | 著者 | 論文 | 誌 | DOI | 何をしたか |
|---|---|---|---|---|---|
| 2016 | Aso Y, Rubin GM | （上記） | eLife | [10.7554/eLife.16135](https://doi.org/10.7554/eLife.16135) | MB099C で PPL1-α3 を含む DAN を光駆動。根拠文 5 本 |
| 2017 | Takemura SY, et al. | A connectome of a learning and memory center in the adult Drosophila brain | eLife | [10.7554/eLife.26975](https://doi.org/10.7554/eLife.26975) | **PPL1-α3 単独の光遺伝学的活性化が、匂い提示なしで接近反応を生じる**ことを示した |
| 2018 | Tsao CH, et al. | （上記） | eLife | [10.7554/eLife.35264](https://doi.org/10.7554/eLife.35264) | PPL1-α3 の活性化が酵母匂い探索を促進 |
| 2021 | Jacob PF, et al. | （上記） | Curr Biol | [10.1016/j.cub.2021.05.056](https://doi.org/10.1016/j.cub.2021.05.056) | **PPL1-α3 の阻害は pre-exposure 由来の忌避を接近に反転させる**。PPL1-α3 だけの阻害では即時の電撃記憶は保たれる。根拠文 3 本 |

**MBON07 → PPL106 の読み出しは明快**: Ichinose 2015 の報酬性長期記憶（α1 区画のループ）と
Jacob 2021 の pre-exposure 学習の反転が、どちらも 1 細胞型の操作で出ている。
さらにこの組はドライバーが最もきれい（下記）。

### この組のドライバーは候補中で最良

| 役割 | 系統 | 発現系 | 標識する型 | 特異性 |
|---|---|---|---|---|
| 上流 MBON07 | **MB319B** / MB319C / MB323B | split-GAL4 | **MBON07 のみ** | 2（理想） |
| 下流 PPL106 | **SL86062** / SL86050 | split-LexA | **PPL106 のみ** | 2（理想） |

**両方が単一タイプの系統で特異性 2。** U × D が 1 組だけになるので、混同の心配が原理的に無い。
`results/mb_final19.csv` では別の組み合わせ（SL86055 + SS01246）で記録しているが、
この組み合わせのほうが明らかに良い。MBON04 ↔ PPL105 は接触の長さと別個体での再現で勝るが、
**ドライバーの質ではこの組が勝つ**。

---

## 3. MBON31 → PAM11（第3候補）

### PAM11 = PAM-α1

| 年 | 著者 | 論文 | 誌 | DOI | 何をしたか |
|---|---|---|---|---|---|
| 2015 | Ichinose T, et al. | Reward signal in a recurrent circuit drives appetitive long-term memory formation | eLife | [10.7554/eLife.10719](https://doi.org/10.7554/eLife.10719) | **MB299B で PAM-α1 を阻害**。条件づけ直後 1 時間の阻害で報酬性長期記憶が損なわれ、22 時間後の阻害では損なわれない |
| 2017 | Takemura SY, et al. | （上記） | eLife | [10.7554/eLife.26975](https://doi.org/10.7554/eLife.26975) | **MB043C × CsChrimson** で PAM-α1 を光駆動し、報酬性記憶の形成・想起時の効果を測定 |
| 2022 | Yamagata N, Imanishi Y, Wu H, Kondo S, Sano H, Tanimoto H | Nutrient responding peptide hormone CCHamide-2 consolidates appetitive memory | Front Behav Neurosci | [10.3389/fnbeh.2022.986064](https://doi.org/10.3389/fnbeh.2022.986064) | **MB299B で PAM-α1 特異的に CCHa2-R をノックダウン**し、長期記憶特異の障害を示した |
| 2021 | Yamagata N, et al. | Presynaptic inhibition of dopamine neurons controls optimistic bias | eLife | [10.7554/eLife.64907](https://doi.org/10.7554/eLife.64907) | PAM の前シナプス抑制。言及 11 回（行動の根拠文は未検出） |

### MBON31 = MBON-α'1a — 行動実験の前例が無い

15 本中、行動の根拠がある論文は 0 本。最も近いのは:

- Felsenberg J, Barnstedt O, Cognigni P, Lin S, Waddell S (2017)
  Re-evaluation of learned information in Drosophila. *Nature* 544:240–244.
  [10.1038/nature21716](https://doi.org/10.1038/nature21716) — 言及 3 回。
  記憶の再評価（extinction）の回路。MBON → DAN のフィードバックを扱う論文群の中核。
- Rubin GM, Aso Y (2024) New genetic tools for mushroom body output neurons in Drosophila.
  *eLife* [10.7554/eLife.90523](https://doi.org/10.7554/eLife.90523)
  — MBON の新しい系統（MBON31 の系統 SS52050 / SS80896 の出所）。形態と系統の記載。

**この組を進めるなら、読み出しは PAM11 側（報酬性長期記憶）に置く。**
MBON31 → PAM11 の新生シナプスが PAM-α1 の活動を変えれば、
Ichinose 2015 の長期記憶課題で差が出ると予想できる。

---

## 4. 却下された組の細胞（参考）

### MBON11 = MBON-γ1pedc>α/β（= MVP2）— 文献は最も多い

107 本。行動の根拠がある論文 5 本。

- Aso Y, et al. (2014) *eLife* [10.7554/eLife.04580](https://doi.org/10.7554/eLife.04580)
  — MB085C / MB112C / MB262B で操作。言及 66 回。
- Tsao CH, et al. (2018) *eLife* [10.7554/eLife.35264](https://doi.org/10.7554/eLife.35264)
  — MB112C。GABA 作動性 MBON-γ1pedc>αβ が β'2 の MBON と MBON-α1 を抑制して酵母探索を促す。
- Chouhan NS, Griffith LC, Haynes P, Sehgal A (2021) Availability of food determines the need
  for sleep in memory consolidation. *Nature* [10.1038/s41586-021-03954-w](https://doi.org/10.1038/s41586-021-03954-w)
  — MB112C。

MBON11 は道具も前例も最も充実しているが、接触が別個体で再現せず（4〜12 µm）、
上流としての出力も γ1 で 95 しかないため候補から外れた。

### MBON33 = MBON-γ2γ3

- Yamada D, Bushey D, Li F, et al. (2023) Hierarchical architecture of dopaminergic circuits
  enables second-order conditioning in Drosophila. *eLife*
  [10.7554/eLife.79042](https://doi.org/10.7554/eLife.79042) — SS67221 で SMP108 を光駆動。

---

## 5. 陽性対照として使う記憶の競合実験（設計済み）

本実験の陽性対照（物理的なシナプス操作ではなく、既知の光遺伝学で同じ読み出しが動くことの確認）。

- Felsenberg J, Barnstedt O, Cognigni P, Lin S, Waddell S (2017)
  Re-evaluation of learned information in Drosophila. *Nature* 544:240–244.
  [10.1038/nature21716](https://doi.org/10.1038/nature21716)
- Felsenberg J, Jacob PF, Walker T, Barnstedt O, Edmondson-Stait AJ, Pleijzier MW,
  Otto N, Schlegel P, Sharifi N, Perisse E, Smith CS, Lauritzen JS, Costa M, Jefferis GSXE,
  Bock DD, Waddell S (2018) Integration of parallel opposing memories underlies memory
  extinction. *Cell* 175:709–722. [10.1016/j.cell.2018.08.021](https://doi.org/10.1016/j.cell.2018.08.021)
- Perisse E, Owald D, Barnstedt O, Talbot CB, Huetteroth W, Waddell S (2016)
  Aversive learning and appetitive motivation toggle feed-forward inhibition in the
  Drosophila mushroom body. *Neuron* 90:1086–1099.
  [10.1016/j.neuron.2016.04.034](https://doi.org/10.1016/j.neuron.2016.04.034)

条件 A（嫌悪 DAN のみ、PI < 0）／条件 B（報酬 DAN のみ、PI > 0）／
条件 C（両方同時、PI ≈ 0）を、CsChrimson で PAM と PPL1 を駆動して確認する設計。
読み出しは Preference Index と、移動速度・旋回。

---

## 6. 再現方法

```
python lit_search.py PPL105 PPL106 MBON07 PAM11 MBON33 MBON04 MBON11 MBON31 FB4N --max-papers 200
```

細胞型と論文上の名前・系統名の対応は `lit_search.py` の `TARGETS` にある。
新しい細胞を追加するときは、neuPrint の instance 名
（例 `PPL105(a'2a2)` → `PPL1-α'2α2`）と `data/driver_lines.csv` の
Descriptive Cell Type 列から別名を作る。取得した本文は `data/literature_cache/` に残る。
