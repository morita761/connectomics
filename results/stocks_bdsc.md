# 使う予定の系統は BDSC にあるか

2026-11-23。[bdsc.indiana.edu](https://bdsc.indiana.edu/) の在庫検索で 1 本ずつ引いた。
対象は [plan_three_pairs.md](plan_three_pairs.md) の 3 組と、
[candidate_MBON04_FB4N.md](candidate_MBON04_FB4N.md)・
[techdemo_MBON31_PAM11.md](techdemo_MBON31_PAM11.md) で使う系統。

## 0. 結論

**どの組も、いま決めている構成のままでは BDSC だけでは揃わない。**
**足りないのは必ず split-LexA の側である。**

| 組 | 上流 | 下流 | BDSC だけで揃うか |
|---|---|---|---|
| **MBON07 → PPL106**（本命） | MB319B **✕** | SL86050 **✕** | **揃わない** |
| MBON02 → PPL106 | SS01337 **✕** | SL86050 **✕** | **揃わない** |
| MBON04 → PPL105 | SL86063 **✕** | **MB058B ◎ 68278** | **揃わない** |
| MBON04 → FB4N | SL86063 **✕** | **SS66004 ◎ 89507** | **揃わない** |
| MBON31 → PAM11 | **SS52050 ◎ 89438** | SL86059 **✕** | **揃わない** |

**split-LexA は BDSC に 1 本も無い。** 検索語 `split-LexA` でも 0 件。
出典の Shuai 2024 本文に、どこから入手するかが書いてある。

> Fly strains can be requested from Janelia or the Bloomington stock center.
> **Split-LexA expression data and fly strains are also available upon request from Aso lab.**
> — [eLife 94168](https://elifesciences.org/articles/94168)

**→ split-LexA は Aso 研（Janelia）に直接依頼する。** これが最初にやることになる。

**ただし朗報がある。** 名前の不一致で見落としていた **BDSC 在庫の代替系統**が見つかった（4 節）。
**本命の上流 MB319B は、BDSC にある MB310C で置き換えられる。**

## 1. 在庫表

◎ = BDSC にあり（番号）、✕ = 検索して 0 件。

### 使う予定のもの

| 系統 | 発現系 | 標識 | 出典 | BDSC |
|---|---|---|---|---|
| **MB319B** | split-GAL4 | MBON07 | 94168 | **✕** |
| MB319C / MB323B | split-GAL4 | MBON07 | 94168 | ✕ / ✕ |
| **SL86050** | **split-LexA** | PPL106 | 94168 | **✕** |
| SL86062 | split-LexA | PPL106 | 94168 | ✕ |
| **SS01337** | split-GAL4 | MBON02 | 94168 | **✕** |
| SS01142 / MB517B | split-GAL4 | MBON02 | 94168 | ✕ / ✕ |
| **SL86063** | **split-LexA** | MBON01+03+04 | 94168 | **✕** |
| **MB058B** | split-GAL4 | PPL1-α'2α2（=PPL105） | Aso 2014 | **◎ 68278** |
| MB059B / MB320B / SS48817 | split-GAL4 | PPL105 ほか | 94168 | ✕ / ✕ / ✕ |
| **SS66004** | split-GAL4 | FB4N | 98405 | **◎ 89507** |
| SS95185 | split-GAL4 | FB4N | 104764 | ✕ |
| **SS52050** | split-GAL4 | MBON31 | 98405 | **◎ 89438** |
| SS80896 | split-GAL4 | MBON31 | 98405 | ◎ 88662 |
| **SL86059** | **split-LexA** | PAM11 | 94168 | **✕** |
| SS67728 / SS67755 | split-GAL4 | MBON04 | 98405 | ◎ 88330 / ◎ 88263 |

### 参照用・対照用

| 系統 | 標識 | BDSC |
|---|---|---|
| MB011B | γ5β'2a + β'2mp + β'2mp_bilateral | ◎ 68294 |
| MB099C | PPL103 + PPL105 | ◎ 68290 |
| MB320C | PPL101（γ1pedc） | ◎ 68253 |
| MB060B | α'2α2 + α'3 + α3 + γ2α'1 | ◎ 68279 |
| MB310C | **MBON-α1（= MBON07）** | **◎ 68313** |

## 2. 規則性 — 出典で決まる

調べていて気づいたが、**在庫の有無は系統の名前ではなく出典で決まる。**

| 出典 | BDSC | 番号帯 |
|---|---|---|
| **Aso 2014**（キノコ体 split-GAL4 コレクション） | **ある** | 682xx〜683xx |
| **Meissner 2024**（eLife 98405、中枢神経系コレクション） | **ある** | 882xx〜895xx |
| **Shuai 2024/2025**（eLife 94168、連合学習コレクション） | **無い** | — |
| Wolff 2024（eLife 104764） | **部分的** | SS66004 はあり、SS95185 は無し |

**本計画がいちばん頼っているのが Shuai 2024 のコレクションで、そこが丸ごと BDSC に無い。**
MBON・DAN を 1 細胞型で狙える系統が多いのはこのコレクションなので、仕方のない面はある。

## 3. split-LexA をどう手に入れるか

**3 組すべてで split-LexA が要る。** どの組も「上流と下流に別の分子を入れる」ので、
GAL4 と LexA の両方が必要で、LexA 側はこのコレクションにしかない。

| 要る split-LexA | どの組で | 代替 |
|---|---|---|
| **SL86050**（PPL106） | MBON07 →・MBON02 → の両方 | **無い**（5 節） |
| **SL86063**（MBON04+01+03） | MBON04 → PPL105・MBON04 → FB4N | 自作なら改善の余地あり（5 節） |
| SL86059（PAM11） | MBON31 → PAM11 | 無い |

**Aso 研への依頼が最初の一歩になる。**
同じ依頼で 3 本まとめて頼めるので、手間は 1 回分で済む。

## 4. 名前の不一致で見落としていた BDSC 在庫

**`data/driver_lines*.csv` には、区画名で書かれた系統が 50 本ある**
（`MBON-α1`、`PPL1-α'2α2` のように）。**照合は hemibrain の型名（`MBON07`、`PPL105`）で
走っているので、この 50 本は一度も候補に出てこなかった。**
`candidate_MBON04_PPL105.md` で MB058B が漏れていたのと同じ原因である。

**しかもこの 50 本はすべて Aso 2014 のコレクション — つまり BDSC にある側である。**
見えていなかったのが、よりによって注文できる系統だった。

拾い直すと、こうなる。

| 目標 | 見落としていた系統 | 標識 | 特異性 | BDSC |
|---|---|---|---|---|
| **MBON07** | **MB310C** | **MBON-α1 のみ** | **2** | **◎ 68313** |
| MBON02 | **MB399B** | **MBON-β2β'2a のみ** | 1 | **◎ 68369** |
| PPL105 | MB058B | PPL1-α'2α2 のみ | 2 | ◎ 68278（既出） |
| PAM11 | **MB043C** | **PAM-α1 のみ** | **2** | **◎ 68363** |
| MBON31 | MB050B | MBON-α'1 ＋ MBON-α2sc | 2 | ◎ 68365 |
| **PPL106** | **MB060B** | α'2α2 ＋ α'3 ＋ **α3** ＋ γ2α'1（4 型） | 2 | ◎ 68279 |

**MB310C は本命の上流をそのまま置き換えられる。**
MBON-α1 だけを標識する quality 1 の系統で、**Tsao 2018 が MBON-α1 特異的として使った実績**もある
（`results/candidate_MBON07_PPL106.md` 6 節）。

**MB399B も MBON02 の上流に使える**（特異性 1 なので SS01337 と同格）。

**PPL106 だけは代わりが無い。** MB060B は 4 型を標識し、
そのうち PPL1-γ2α'1 と PPL1-α'2α2 は MBON07 のすぐ近くを通っている
（それぞれ 62 µm、255 µm — `plan_three_pairs.md` 8.5 節）。**使えない。**

## 5. 組ごとの入手計画

### MBON07 → PPL106（本命）

| 役割 | 第一希望 | BDSC で代替 |
|---|---|---|
| 上流 MBON07〈GAL4〉 | MB319B | **MB310C ◎ 68313** |
| 下流 PPL106〈LexA〉 | **SL86050** | **無い。Aso 研に依頼** |

**上流は BDSC で足りる。下流の split-LexA だけが律速。**

### MBON02 → PPL106

| 役割 | 第一希望 | BDSC で代替 |
|---|---|---|
| 上流 MBON02〈GAL4〉 | SS01337 | **MB399B ◎ 68369** |
| 下流 PPL106〈LexA〉 | SL86050 | 無い（本命と共用） |

**本命と同じ SL86050 を使うので、依頼は 1 本で 2 組ぶん。**

### MBON04 → PPL105 ／ MBON04 → FB4N

| 役割 | 第一希望 | BDSC で代替 |
|---|---|---|
| 上流 MBON04〈LexA〉 | **SL86063**（+MBON01,03） | 無い。Aso 研に依頼 |
| 下流 PPL105〈GAL4〉 | **MB058B ◎ 68278** | — |
| 下流 FB4N〈GAL4〉 | **SS66004 ◎ 89507** | — |

**下流は両方 BDSC にある。上流の split-LexA だけが足りない。**

> **自作する手もある。** Shuai 2024 は
> 「split-GAL4 の DBD を LexA の DBD に置き換えれば split-LexA になる」と書いている。
> **MBON04 だけを標識する split-GAL4（SS67728 ◎ 88330 / SS67755 ◎ 88263）を変換すれば、
> MBON01・MBON03 を含まない「きれいな MBON04 split-LexA」が作れる。**
> これは入手の回避策であると同時に、**MBON04 → FB4N で見つかった
> 「MBON01 が FB4N と 4.3 µm 接している」問題をまるごと消す**
> （`candidate_MBON04_FB4N.md` 4 節）。やる価値はある。

### MBON31 → PAM11

| 役割 | 第一希望 | BDSC で代替 |
|---|---|---|
| 上流 MBON31〈GAL4〉 | **SS52050 ◎ 89438** | SS80896 ◎ 88662 |
| 下流 PAM11〈LexA〉 | **SL86059** | 無い。Aso 研に依頼 |

## 6. まとめ — 次にやること

1. **Aso 研（Janelia）に split-LexA を依頼する。** `SL86050`・`SL86063`・`SL86059` の 3 本。
   論文の記載どおり「upon request」なので、まとめて 1 通で足りる。
2. **並行して BDSC に注文できるものを出す。**
   **68313（MB310C、本命の上流）**・68278（MB058B）・89507（SS66004）・89438（SS52050）。
   68369（MB399B）と 88330/88263（SS67728/SS67755）も要るなら一緒に。
3. **MBON04 split-LexA の自作を検討する。** SS67728 から DBD を差し替える。
   入手の保険であると同時に、混同がひとつ消える。
4. **解析側を直す。** 区画名で書かれた 50 本が照合から漏れている。
   `match_drivers.py` に区画名と hemibrain 型名の対応表を足せば、
   今回のような見落としは無くなる。

## 使ったデータ

- 在庫: [BDSC Stock Search](https://bdsc.indiana.edu/Home/Search)。系統名そのもので検索した。
  見つからなかったものは「0 件」と表示されたことを確認している。
  **BDSC が別名で持っている可能性までは排除できない。**
- split-LexA の入手先: [Shuai et al., eLife 94168](https://elifesciences.org/articles/94168) 本文。
- 系統の中身・出典: `data/driver_lines.csv`、`data/driver_lines_cns.csv`。
- 近傍の接触: [plan_three_pairs.md](plan_three_pairs.md) 8.5 節。
