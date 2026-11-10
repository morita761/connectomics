# GMR-LexA を候補ペアに照合した結果

2026-11-10。再現は `python gmr_lexa_match.py`。
一覧は `results/gmr_lexa_vfb.csv`（VFB の注釈）と `results/gmr_lexa_match.csv`（論文の共起）。

## 0. 結論

> **【2026-11-10 追記】この報告書の結論は、系統まるごとの GMR-LexA についてのもの。
> その後ヘミドライバー断片まで降りて調べたところ、9 組中 5 組が救えることが分かった。
> → `results/hemidriver_lexa_report.md` を先に読むこと。**

**今ある GMR-LexA では、LexA が無くて落ちた 32 組のどれも救えない。**
キノコ体の細胞には GMR-LexA の代替が複数見つかったが、それらは既に split-LexA がある細胞ばかりで、
新しい候補は増えなかった。

| 問い | 答え |
|---|---|
| GMR-LexA は実在するか | **150 系統**を FlyBase / VirtualFlyBrain で確認（約 1,400 系統の一部） |
| そのうち標識細胞が分かっているものは | **286 の発現パターン**に注釈あり。型に読み替えられたのは MBON・PAM の一部 |
| 落ちた 32 組の細胞を標識する GMR-LexA はあるか | **0 件**。VFB の注釈にも論文本文にも出てこない |
| キノコ体の候補に使えるものはあるか | **ある**（下表）。ただし既存の split-LexA の代替にとどまる |

## 1. 何を照合したか

これまでの探索で、**「上流にも下流にも split-GAL4 はあるのに、どちらにも LexA が無いので
発現系が 2 つにならない」という理由だけで落ちた組が 32 組**あった。
GMR-LexA（Rubin のエンハンサー断片 LexA、約 1,400 系統）が使えれば、この壁は外れるはずだった。

照合の手順:

1. **GMR-LexA が実在する系統名を確定する** — FlyBase の allele 名 `Ecol\lexA[GMRxxxxx]` を
   VirtualFlyBrain の solr から取得。**150 系統**（`data/gmr_lexa_lines.txt`）。
2. **その系統が何を標識するかを調べる** — VFB の neo4j に、
   `P{GMRxxxxx-GAL4/lexA} expression pattern` に `part_of` でぶら下がる
   Gen1 MCFO の単一細胞を問い合わせ、その細胞型（FlyBase 解剖学名）を取る。
   LexA 版は GAL4 版と**同じエンハンサー断片**なので、GAL4 版の注釈も証拠として使える。
3. **論文本文でも裏を取る** — `data/literature_cache/` の全文を文単位で走査し、
   「GMR-LexA が実在する R 系統名」と「細胞型の名前」が同じ文に出る組を拾う。

## 2. VFB の注釈で分かったこと

286 の発現パターンに 714 件の細胞注釈があり、そのうち hemibrain の型名に読み替えられた 43 件は次のとおり。
**太字は LexA 版が実在する系統。**

| 型 | 標識する GMR 系統（LexA 版があるもの） | 既存の split-LexA |
|---|---|---|
| MBON01 | **R14C08, R15B01, R21D02** | SL86047, SL86046 |
| MBON02 | R21D02（GAL4 のみ注釈） | — |
| MBON03 | **R14C08, R15B01, R21D02** | SL86063 経由 |
| **MBON04** | **R14C08, R15B01** | **SL86063** |
| MBON05 | **R24E12** | — |
| MBON11 | **R12G04** | SL86058 |
| MBON19 | R71D08 | — |
| PAM01 | **R48B04, R58E02** | — |
| PAM02 | **R48B04** | — |
| PAM05 | **R14C08, R48B04** | — |
| PAM06 | **R14C08, R48B04** | — |
| PAM08 | **R48B04** | — |
| PAM10 | **R11C07** | — |
| **PAM11** | **R71C03, R72D01** | **SL86065, SL86059, SL86053** |
| PAM13 | **R24E12**, R58E02 | — |
| PAM14 | **R24E12**, R58E02 | — |

### MBON04 の代替: R15B01-LexA

**R15B01-LexA は MBON01・MBON03・MBON04 を標識する。これは SL86063 とまったく同じ組み合わせ**。
split-LexA が手に入らないときの代替になるが、U の中身が同じなので混同の判定は変わらない。

**R14C08-LexA は使いにくい**。MBON01・MBON03・MBON04 に加えて **PAM05・PAM06 も標識する**。
PAM05 → MBON04 は 1,513 シナプス、PAM06 → MBON04 は 1,749 シナプスあるので、
U の中に強い既存結合が入り込む。下流が MBON04 側になる向きでは使えない。

### PAM11 の代替: R71C03-LexA / R72D01-LexA

MBON31 → PAM11 の下流側は既に SL86059 / SL86053（特異性 2）があるので、
GMR-LexA を使う理由はない。R72D01 は MBON07 とも同じ文献に出てくる（α1 区画）ので、
**MBON07 と PAM11 の両方を標識する恐れがあり、むしろ避けるべき**。

### 論文本文での共起は 2 件、どちらも偽陽性

全文 365 本を文単位で走査して拾えたのは 2 件だけで、中身を読むと両方とも使えない。

| 拾った文 | 判定 |
|---|---|
| 「R58E02 と VT999036 の投射が γ4・γ5 で重なり、MBON21 の細胞に対応する」(PMC11611264) | **✕ 標識ではない**。R58E02 は PAM の系統で、MBON21 と同じ区画に**投射が重なる**と言っているだけ |
| 「R12G04-lexA, pJFRC216-13xLexAop2-… と PPL1-γ1pedc（マゼンタ）」(PMC4273436) | **✕ 図の凡例**。系統名と別の染色が同じ文に並んでいるだけ |

語の共起だけでは「その系統がその細胞を標識する」とは言えない例で、
**実質的な証拠は VFB の注釈のほうしかない**。

## 3. 落ちた 32 組は救えなかった

落ちた 32 組に出てくる細胞型 49 種類について、VFB の注釈は **0 件**、
論文本文の共起も上の 2 件（どちらも偽陽性）だけで、
GMR-LexA との結びつきが見つからなかった。

特に、別個体でも接触が再現し上流もその領域で出力している**上位 9 組**は次のとおりで、
どれも LexA の当てが無いまま残る。

| 組 | 接触 | 別個体（最小〜最大） | 主領域 | 下流のドライバー |
|---|---|---|---|---|
| MBON31 → SMP503 | 152 µm | 25〜55 µm | SMP(R) | SS32246〈GAL4〉 |
| VES012 → DNb05 | 145 µm | **35〜62 µm** | VES(R) | SS01051〈GAL4〉 |
| MBON30 → oviIN | 142 µm | **36〜68 µm** | CRE(R) | SS65422 / SS65423〈GAL4〉 |
| MBON21 → FB4X | 133 µm | **36〜63 µm** | CRE(R) | SS70711〈GAL4〉 |
| ExR2 → LAL098 | 103 µm | **47〜77 µm** | LAL(R) | SS49882〈GAL4〉 |
| AVLP478 → DNp06 | 95 µm | 28〜39 µm | AVLP(R) | SS02256 / SS01047〈GAL4〉 |
| PLP022 → DNp10 | 51 µm | 30〜51 µm | PLP(R) | SS01049 / SS01580〈GAL4〉 |
| PS161 → DNa04 | 62 µm | 10〜21 µm | SPS(R) | SS02384〈GAL4〉 |
| DL4_adPN → LHPV6c1 | 41 µm | 27〜37 µm | LH(R) | LH2220〈GAL4〉 |

**この 9 組は、別個体での接触だけ見れば全脳の本命 MBON04 → FB4N（10〜29 µm）を上回る。**
ExR2 → LAL098 の 47〜77 µm、MBON30 → oviIN の 36〜68 µm は、
キノコ体の第一候補 MBON04 ↔ PPL105（33〜38 µm）より良い。
**発現系が 1 つ足りないという理由だけで落ちているのは、もったいない。**

## 4. なぜ見つからないのか（結果の読み方）

**「GMR-LexA が無い」とは言い切れない。「分かっていない」が正しい。**

| 限界 | 中身 |
|---|---|
| 系統の網羅 | 確認できた GMR-LexA は 150 系統。実際には**約 1,400 系統**ある |
| 注釈の網羅 | 細胞型の注釈があるのは 286 の発現パターンだけ |
| 注釈の粒度 | FlyBase の解剖学名は粗い。「dopaminergic PPL1 neuron」止まりで PPL105 と PPL106 を区別しない |
| 型名の違い | LAL098・DNb05・SMP503 などは hemibrain のコネクトーム由来の型名で、FlyBase の解剖学には対応する語が無い |

つまり **「hemibrain の型名で呼ばれる細胞は、Gen1 の注釈体系に載っていない」**のが
見つからない主因であって、断片が標識していない証拠ではない。

## 5. 次にやるべきこと（有望な順）

### ① split-GAL4 → split-LexA の変換（最有力）

落ちた 9 組は、**下流に特異性の高い split-GAL4 が既にある**。
Shuai 2024 は split-GAL4 の DNA 結合ドメインを LexA に差し替える変換を
**34 系統試して 22 系統成功（65%）**させている。
変換すべき系統は下の 9 本（上の表の下流ドライバー）:

```
SS32246  SS01051  SS65422/SS65423  SS70711  SS49882
SS02256/SS01047  SS01049/SS01580  SS02384  LH2220
```

変換は「作れば解ける」種類の作業で、照合で解決する問題ではない。

### ② split-GAL4 のヘミドライバー断片を調べる → **実施済み。5 組が救えた**

`python hemidriver_lexa.py` / `results/hemidriver_lexa_report.md`。
以下は実施前の想定。

split-GAL4 は 2 つのエンハンサー断片（AD 側と DBD 側）の積で作られている。
例えば `SS49882 = Rxxxxx-p65ADZp ∩ VTxxxxxx-ZpGdbd` のような構成。
**その AD 側か DBD 側の断片に GMR-LexA 版があれば、その断片はその細胞を含む集合を標識する。**
特異性は split より落ちる（積ではなく片側だけになる）が、U を画像で確定できれば使える。
Janelia の split-GAL4 の遺伝子型表（Meissner 2024 eLife 98405 の補足表）から
ヘミドライバーを取り出し、`data/gmr_lexa_lines.txt` と突き合わせるのが次の一手。

### ③ FlyLight の画像で直接確認する

150 系統それぞれの発現画像を見て、候補の細胞が入っているかを判定する。
VFB の各系統ページから FlyLight / Gen1 MCFO / NeuronBridge へのリンクがある。

```
http://flweb.janelia.org/cgi-bin/view_flew_imagery.cgi?line=R15B01
https://gen1mcfo.janelia.org/cgi-bin/view_gen1mcfo_imagery.cgi?line=R15B01
```

NeuronBridge なら hemibrain の骨格から光学画像を検索できるので、
**「LAL098 に似た形の細胞を含む Gen1 系統」を逆引きできる**。これが
hemibrain 型名と Gen1 系統を結ぶ唯一の自動手段になりうる。

## 6. 使ったデータ

- GMR-LexA の系統名: VirtualFlyBrain solr
  `https://solr.virtualflybrain.org/solr/ontology/select?q=label:*lexA*GMR*`
  → FlyBase の allele `Ecol\lexA[GMRxxxxx]` 150 件
- 細胞型の注釈: VirtualFlyBrain neo4j `https://pdb.virtualflybrain.org/db/data/transaction/commit`
  （Gen1 MCFO の単一細胞 → `part_of` → 発現パターン → `expresses` → 系統）
- 論文本文: `data/literature_cache/`（Europe PMC のオープンアクセス全文）
- 候補ペア: `results/shortlist_final.csv`（306 組）
