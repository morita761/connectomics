# ヘミドライバー断片の LexA 版を使って、落ちた組を救う

2026-11-10。再現は `python hemidriver_lexa.py`。一覧は `results/hemidriver_lexa.csv`。

## 0. 結論

**9 組中 5 組が救える。**

系統まるごとの GMR-LexA では 1 組も救えなかったが（`results/gmr_lexa_report.md`）、
**split-GAL4 を構成する 2 つの断片のうち片方に LexA 版があるか**を見ると話が変わる。

| | 結果 |
|---|---|
| 調べた split-GAL4 | 69 系統（LexA が無くて落ちた 32 組のドライバー） |
| ヘミドライバーが分かったもの | **60 系統** |
| 片方の断片に LexA 版があるもの | **7 系統** |
| それで救える組 | **6 組**（うち接触が再現し上流も出力している上位の組が **5 組**） |

## 1. 考え方

split-GAL4 は 2 つのエンハンサー断片の**積**でできている。

```
SS49882 = VT045661-GAL4.DBD  ∩  VT039475-p65.AD
          ↑ この断片に LexA 版（Ecol\lexA[VT045661]）がある
```

積の片側だけを使うと、標識する細胞は**増える**（特異性は落ちる）。
しかし「LexA が無いので発現系が 2 つにならない」という壁は、それで外れる。
元の split-GAL4 がその細胞を標識している以上、**片側の断片も必ずその細胞を含む**。

## 2. 救える 5 組

| 判定 | 組 | 接触 | 別個体 | LexA にする側 | 使う断片 |
|---|---|---|---|---|---|
| ★ | **MBON30 → oviIN** | 142 µm | **36〜68 µm** | 下流 oviIN（SS65423） | **VT026347-LexA**（AD 側） |
| ★ | **MBON21 → FB4X** | 133 µm | **36〜63 µm** | 上流 MBON21（SS81353） | **GMR11C07-LexA**（DBD 側） |
| △ | **ExR2 → LAL098** | 103 µm | **47〜77 µm** | 下流 LAL098（SS49882） | **VT045661-LexA**（DBD 側） |
| ★ | **AVLP478 → DNp06** | 95 µm | 28〜39 µm | 下流 DNp06（SS02256 / SS01047） | **VT019018-LexA**（AD 側） |
| ★ | **DL4_adPN → LHPV6c1** | 41 µm | 27〜37 µm | 上流 DL4_adPN（SS10228 / SS32214） | **VT033006-LexA**（AD 側） |
| | LAL098 → DNb01 | 82 µm | 2〜4 µm | 上流 LAL098（SS49882） | VT045661-LexA |

最後の 1 組（LAL098 → DNb01）は別個体での接触が 2〜4 µm しかないので落ちたまま。

> **【追記】ExR2 ↔ LAL098 は EM で落ちた。** 膜の密着が、同じ ExR2 が 44 シナプスを作っている
> LAL184 との対照の 1/6〜1/10 しかない（直接隣接 1,348〜3,692 対 8,153、最接近の隙間 45 nm 対 0 nm）。
> 接触の長さが再現しても、膜が接しているとは限らないという例。
> 詳細と、向きを逆（LAL098 → ExR2）にしたほうが良い理由は `results/candidate_ExR2_LAL098.md`。
> **残り 4 組も同じ EM の測定にかけてから順位をつけるべき。**

**別個体での接触だけで見ると、上の 5 組のうち 3 組はこれまでの本命を上回る。**

| 候補 | 別個体での接触 |
|---|---|
| **ExR2 → LAL098** | **47〜77 µm** |
| **MBON30 → oviIN** | **36〜68 µm** |
| **MBON21 → FB4X** | **36〜63 µm** |
| MBON04 ↔ PPL105（キノコ体の第一候補） | 33〜38 µm |
| MBON04 → FB4N（全脳の第一候補） | 10〜29 µm |

### 断片の内訳

| split-GAL4 | 標的 | 構成 | LexA 版がある側 |
|---|---|---|---|
| SS65423 | oviIN | VT026035(DBD) ∩ **VT026347**(AD) | AD |
| SS81353 | MBON21 | **GMR11C07**(DBD) ∩ VT000877(AD) | DBD |
| SS49882 | LAL098 | **VT045661**(DBD) ∩ VT039475(AD) | DBD |
| SS02256 | DNp06 | VT017411(DBD) ∩ **VT019018**(AD) | AD |
| SS01047 | DNp06 | VT017647(DBD) ∩ **VT019018**(AD) | AD |
| SS10228 | DL4_adPN | GMR73F07(DBD) ∩ **VT033006**(AD) | AD |
| SS32214 | DL4_adPN | GMR73F02(DBD) ∩ **VT033006**(AD) | AD |

DNp06 は **2 つの独立な split-GAL4（SS02256 と SS01047）が同じ VT019018 を AD 側に持つ**。
別々に作られた系統が同じ断片に行き着いているので、この断片が DNp06 を含むことの確かさは高い。
DL4_adPN も同様に SS10228 と SS32214 が VT033006 を共有している。

## 3. 混同（U × D）は、分かる範囲ではすべて 0

断片だけを使うと標識する細胞の集合 U（または D）が広がるので、
**U × D のすべての組が 0** という条件を取り直す必要がある。
いま中身が分かっている分について neuPrint で数え直した。

| 組 | hemibrain | male-CNS | 備考 |
|---|---|---|---|
| MBON21 ↔ FB4X | 0 / 0 | 0 / 0 | 目標 |
| **PAM10 ↔ FB4X** | **0 / 0** | **0 / 0** | GMR11C07 は PAM10 も標識する（VFB の注釈） |
| MBON30 ↔ oviIN | 0 / 0 | 0 / 0 | 目標 |
| ExR2 ↔ LAL098 | 0 / 0 | 0 / 0 | 目標 |
| AVLP478 ↔ DNp06 | 0 / 0 | 0 / 0 | 目標 |
| **MBON11 ↔ DNp06** | **0 / 0** | **0 / 0** | SS01181（AVLP478 の GAL4）は MBON11 も標識する |
| DL4_adPN ↔ LHPV6c1 | 0 / 0 | 0 / 0 | 目標 |

**分かっている混同はすべて 0**。ただし下の限界を読むこと。

## 4. 限界 — ここから先は画像を見ないと決まらない

**断片が標識する細胞の全体（U または D の中身）は、まだ分かっていない。**

- VFB に細胞型の注釈があるのは **GMR11C07 だけ**（PAM10）。
  VT045661・VT033006・VT026347・VT019018 には注釈が無い。
- したがって上の「混同 0」は、**分かっている相手に限った話**である。
  断片が他に何十という細胞型を標識していれば、その中に既存結合があるかもしれない。

確かめるには FlyLight の発現画像を見て、標識細胞を数え上げるしかない。

```
http://flweb.janelia.org/cgi-bin/view_flew_imagery.cgi?line=VT045661
http://flweb.janelia.org/cgi-bin/view_flew_imagery.cgi?line=VT019018
http://flweb.janelia.org/cgi-bin/view_flew_imagery.cgi?line=VT026347
http://flweb.janelia.org/cgi-bin/view_flew_imagery.cgi?line=VT033006
http://flweb.janelia.org/cgi-bin/view_flew_imagery.cgi?line=R11C07
```

NeuronBridge なら hemibrain の骨格から光学画像を逆引きできるので、
「この断片の発現に、どのコネクトームの細胞が含まれるか」を半自動で出せる。
VFB にはコネクトームのニューロンと系統の発現パターンを結ぶ
`has_similar_morphology_to_part_of` という関係があり（NBLAST の結果）、
実際に DNb05 と SS01051 の発現パターンが結ばれている。**これを使えば U の推定を自動化できる。**

## 5. 実務上の注意

1. **断片の LexA 版は活性が違うことがある**。Pfeiffer 2010 は
   「LexA は GAL4 より転写活性化因子として弱い」と報告している。
   表記も `lexA::GAD` / `lexA::p65` / `lexA::VP16` と混在しており、強さが異なる。
   今回見つかった 5 つは `Ecol\lexA[VT045661]`、`Ecol\lexA[p65.VT033006]`、
   `Ecol\lexA[GAD.VT026347]`、`Ecol\lexA[GAD.VT019018]`、`Ecol\lexA[GMR11C07]` の形。
2. **在庫の確認が必要**。FlyBase に allele として登録があっても、
   配布系統があるとは限らない。BDSC / VDRC に照会すること。
3. **split-LexA への変換と比べる**。断片 1 つの LexA は特異性が落ちるので、
   Shuai 2024 の方法で split-GAL4 をそのまま split-LexA に変換するほうが
   特異性は保てる（成功率 22/34）。**断片の LexA は「今すぐ試せる」利点**がある。

## 6. 使ったデータ

- ヘミドライバーの構成: VirtualFlyBrain neo4j
  `MATCH (n:Individual:Split)-[:INSTANCEOF]->(c)` — split-GAL4 の試料が属する
  「P{A-GAL4.DBD} ∩ P{B-p65.AD} expression pattern」クラスから断片名を取る。
- LexA 版がある断片: VirtualFlyBrain solr の `Ecol\lexA[...]` allele 名 **190 件**
  （`GMRxxxxx` / `VTxxxxxx`、`GAD.` `p65.` `DBD.` などの接頭辞を吸収）。
- 候補ペア: `results/shortlist_final.csv`
- 混同の数え直し: neuPrint hemibrain:v1.2.1 / male-cns:v1.0
