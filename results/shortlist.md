# 強制シナプス形成に向けた最終候補（shortlist.py）

条件: 合計シナプス 0 個以下。「接触の再現」は別個体の各半球で 10 µm 以上の接触があること。
シナプス数は hemibrain の順/逆と、照合できた各個体の両方向の合計。

| 領域 | 上流 → 下流 | hemibrain 接触 | シナプス合計 | 照合した個体 | 別個体の接触 | 上流ドライバー | 下流ドライバー |
|---|---|---|---|---|---|---|---|
| 上前大脳 / クランプ | MBON31 → SMP503 | 152 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 55µm, male-cns R 25µm | SS52050〈GAL4〉[特異性2] / SS80896〈GAL4〉[特異性2] | SS32246〈GAL4〉[特異性1] / SS01239〈GAL4〉(+MBON11)[特異性0] |
| 側方付属葉 / 楔形領域 | VES012 → DNb05 | 145 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 35µm, male-cns R 62µm | SS50738〈GAL4〉[特異性2] | SS01051〈GAL4〉[特異性2] / SS01051〈GAL4〉[特異性2] |
| クランプ | MBON30 → oviIN | 142 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 68µm, male-cns R 36µm | SS76094〈GAL4〉[特異性2] | SS65422〈GAL4〉[特異性1] / SS65423〈GAL4〉[特異性1] |
| キノコ体の中だけ | PPL106 → PAM12 | 137 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 25µm, male-cns R 14µm | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] | MB443B〈GAL4〉[特異性2] / MB029B〈GAL4〉[特異性2] |
| クランプ | MBON21 → FB4X | 133 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 63µm, male-cns R 36µm | SS46348〈GAL4〉[特異性2] / SS81353〈GAL4〉[特異性2] | SS70711〈GAL4〉[特異性1] / SS70711〈GAL4〉 |
| 上前大脳 / LexA 12 型を軸に再探索 | MBON04 → FB4N | 109 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 29µm, male-cns R 10µm | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] | SS66004〈GAL4〉[特異性1] / SS95185〈GAL4〉 |
| 中心複合体 | ER3a_d → FB4L | 104 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 11µm, male-cns R 22µm | SS04710〈GAL4〉[特異性1] / SS04710〈GAL4〉 | SS48802〈GAL4〉[特異性2] / SS61133〈GAL4〉[特異性2] |
| 側方付属葉 | ExR2 → LAL098 | 103 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 47µm, male-cns R 77µm | SS58826〈GAL4〉[特異性2] / SS59679〈GAL4〉[特異性2] | SS49882〈GAL4〉(+VES019)[特異性2] / SS49882〈GAL4〉(+VES019)[特異性1] |
| 腹外側前大脳 | AVLP478 → DNp06 | 95 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 39µm, male-cns R 28µm | SS01181〈GAL4〉(+MBON11)[特異性0] | SS02256〈GAL4〉[特異性2] / SS01047〈GAL4〉[特異性1] |
| キノコ体 / LexA 12 型を軸に再探索 / キノコ体の中だけ | PPL105 → MBON04 | 83 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 24µm, male-cns R 14µm | MB059B〈GAL4〉[特異性1] / MB320B〈GAL4〉(+PPL101)[特異性1] | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] |
| 後斜面 | PS161 → DNa04 | 62 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 10µm, male-cns R 21µm | SS22505〈GAL4〉(+PS159)[特異性1] | SS02384〈GAL4〉(+DNa10)[特異性1] |
| 腹外側前大脳 | PLP022 → DNp10 | 51 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 30µm, male-cns R 51µm | SS73660〈GAL4〉(+PLP015,WEDPN2B?)[特異性2] | SS01049〈GAL4〉[特異性2] / SS01580〈GAL4〉[特異性2] |
| キノコ体 | DL4_adPN → LHPV6c1 | 41 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 27µm, male-cns R 37µm | SS10228〈GAL4〉(+VA1d_adPN)[特異性1] / SS32214〈GAL4〉(+VA1d_adPN)[特異性0] | LH2220〈GAL4〉[特異性1] |
| キノコ体の中だけ | MBON07 → PPL106 | 37 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 19µm, male-cns R 19µm | MB319B〈GAL4〉[特異性2] / MB319C〈GAL4〉[特異性2] | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] |
| キノコ体の中だけ | MBON31 → PAM11 | 25 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 23µm, male-cns R 18µm | SS52050〈GAL4〉[特異性2] / SS80896〈GAL4〉[特異性2] | SL86059〈LexA〉[特異性2] / SL86053〈LexA〉[特異性2] |
| キノコ体の中だけ | MBON02 → PPL106 | 25 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 19µm, male-cns R 25µm | MB517B〈GAL4〉[特異性2] / SS01337〈GAL4〉[特異性1] | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] |
| 上前大脳 / クランプ | MBON31 → SMP077 | 258 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 120µm, male-cns R 61µm | SS52050〈GAL4〉[特異性2] / SS80896〈GAL4〉[特異性2] | - |
| 上前大脳 / クランプ | MBON35 → ATL006 | 211 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 104µm, male-cns R 70µm | - | - |
| 後斜面 | LAL016 → PS233 | 174 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 66µm, male-cns R 81µm | - | - |
| 後斜面 | PS137 → DNa05 | 167 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 34µm, male-cns R 95µm | - | SS95335〈GAL4〉[特異性2] / SS02536〈GAL4〉[特異性2] |
| 上前大脳 | LHPV5l1 → SMP237 | 165 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 36µm, male-cns R 59µm | - | - |
| 後斜面 | LAL099 → PLP245 | 153 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 91µm, male-cns R 69µm | - | - |
| 上前大脳 | LHCENT3 → SLP209 | 152 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 82µm, male-cns R 96µm | - | - |
| 上前大脳 | LHCENT2 → SMP353 | 151 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 55µm, male-cns R 40µm | - | SS81072〈GAL4〉[特異性2] / SS81071〈GAL4〉[特異性2] |
| クランプ | MBON31 → LHCENT8 | 148 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 120µm, male-cns R 125µm | SS52050〈GAL4〉[特異性2] / SS80896〈GAL4〉[特異性2] | - |
| 上前大脳 / クランプ | MBON32 → SMP156 | 145 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 84µm, male-cns R 56µm | - | - |
| 上前大脳 | LHPV5e2 → SMP337 | 141 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 25µm, male-cns R 18µm | - | - |
| 後斜面 | PS156 → DNp19 | 140 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 86µm, male-cns R 54µm | - | SS81965〈GAL4〉[特異性2] / SS12653〈GAL4〉(+DNp15)[特異性1] |
| 上前大脳 / クランプ / LexA 12 型を軸に再探索 | MBON01 → CRE050 | 139 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 18µm, male-cns R 19µm | SS01308〈GAL4〉[特異性2] / SL86047〈LexA〉[特異性2] | - |
| 上前大脳 | LHPV6m1 → SMP527 | 133 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 75µm, male-cns R 95µm | - | SS32144〈GAL4〉[特異性2] |
| 上前大脳 | MBON33 → SMP199 | 132 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 70µm, male-cns R 69µm | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] | - |
| クランプ | MBON21 → LHPV5e3 | 128 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 46µm, male-cns R 39µm | SS46348〈GAL4〉[特異性2] / SS81353〈GAL4〉[特異性2] | - |
| 上前大脳 / クランプ | MBON33 → SMP156 | 126 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 113µm, male-cns R 55µm | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] | - |
| 上前大脳 | LHAV2o1 → SLP011 | 124 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 46µm, male-cns R 42µm | - | - |
| 上前大脳 | LHCENT9 → CRE011 | 122 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 30µm, male-cns R 39µm | - | SS45245〈GAL4〉[特異性2] / SS45222〈GAL4〉[特異性2] |
| 側方付属葉 | CRE041 → DNb01 | 121 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 68µm, male-cns R 41µm | - | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| 上前大脳 / クランプ | MBON27 → CRE077 | 121 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 60µm, male-cns R 44µm | - | - |
| 上前大脳 | LHAV3a1_c → SLP458 | 120 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 50µm, male-cns R 37µm | - | - |
| 側方付属葉 | LAL159 → DNa01 | 119 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 22µm, male-cns R 27µm | - | SS00731〈GAL4〉[特異性1] / SS00731〈GAL4〉[特異性0] |
| 上前大脳 / クランプ / LexA 12 型を軸に再探索 | MBON01 → CRE022 | 117 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 36µm, male-cns R 13µm | SS01308〈GAL4〉[特異性2] / SL86047〈LexA〉[特異性2] | - |
| 上前大脳 | MBON33 → SMP272 | 116 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 51µm, male-cns R 58µm | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] | - |
| クランプ | MBON21 → ATL017 | 115 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 35µm, male-cns R 29µm | SS46348〈GAL4〉[特異性2] / SS81353〈GAL4〉[特異性2] | - |
| 上前大脳 | MBON31 → CRE069 | 112 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 87µm, male-cns R 47µm | SS52050〈GAL4〉[特異性2] / SS80896〈GAL4〉[特異性2] | - |
| 側方付属葉 / 楔形領域 | LAL142 → DNa03 | 112 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 75µm, male-cns R 34µm | - | - |
| 上前大脳 | LHCENT3 → CRE045 | 111 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 81µm, male-cns R 77µm | - | - |
| 上前大脳 | MBON26 → SMP180 | 110 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 26µm, male-cns R 24µm | SS80958〈GAL4〉[特異性1] | - |
| 上前大脳 / クランプ | MBON35 → AL-MBDL1 | 107 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 34µm, male-cns R 34µm | - | - |
| 上前大脳 | SMP085 → DNp32 | 105 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 39µm, male-cns R 49µm | - | SS02394〈GAL4〉[特異性1] / SS02394〈GAL4〉[特異性0] |
| 側方付属葉 | LAL045 → DNb01 | 97 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 21µm, male-cns R 16µm | - | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| 外側角 | M_vPNml87 → LHAV2m1 | 96 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 41µm, male-cns R 48µm | SS00587〈GAL4〉[特異性1] / SS00586〈GAL4〉[特異性1] | - |
| 外側角 | M_vPNml86 → LHPD3a5 | 93 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 16µm, male-cns R 47µm | - | - |
| 側方付属葉 | ExR3 → LAL195 | 89 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 34µm, male-cns R 46µm | SS16969〈GAL4〉[特異性1] / SS82798〈GAL4〉 | - |
| 上前大脳 / LexA 12 型を軸に再探索 | MBON04 → SIP067 | 88 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 30µm, male-cns R 51µm | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] | - |
| クランプ | MBON31 → CRE072 | 88 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 64µm, male-cns R 23µm | SS52050〈GAL4〉[特異性2] / SS80896〈GAL4〉[特異性2] | - |
| 後斜面 | LAL099 → AOTU049 | 86 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 70µm, male-cns R 31µm | - | - |
| クランプ | MBON27 → ATL027 | 85 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 41µm, male-cns R 26µm | - | - |
| 側方付属葉 | LAL109 → DNa02 | 85 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 66µm, male-cns R 64µm | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 上前大脳 | MBON26 → SMP049 | 80 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 20µm, male-cns R 21µm | SS80958〈GAL4〉[特異性1] | - |
| クランプ | MBON29 → CRE022 | 80 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 22µm, male-cns R 30µm | SS67662〈GAL4〉[特異性2] / SS77450〈GAL4〉[特異性2] | - |
| 上前大脳 | MBON26 → SIP069 | 79 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 22µm, male-cns R 37µm | SS80958〈GAL4〉[特異性1] | - |
| 後斜面 | IB058 → DNp27 | 79 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 10µm, male-cns R 17µm | - | SS01876〈GAL4〉[特異性2] / MB202B〈GAL4〉[特異性2] |
| 後斜面 | LAL018 → PS140 | 78 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 32µm, male-cns R 17µm | - | SS72289〈GAL4〉[特異性2] |
| 上前大脳 | MBON32 → SMP054 | 78 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 26µm, male-cns R 33µm | - | - |
| クランプ | MBON29 → CRE004 | 76 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 17µm, male-cns R 21µm | SS67662〈GAL4〉[特異性2] / SS77450〈GAL4〉[特異性2] | - |
| 外側角 | DC2_adPN → LHAV5a8 | 73 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 29µm, male-cns R 36µm | - | - |
| 後斜面 / 腹外側前大脳 | PLP208 → DNp05 | 72 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 45µm, male-cns R 51µm | - | SS00865〈GAL4〉[特異性2] / SS00865〈GAL4〉[特異性1] |
| 腹外側前大脳 | AVLP432 → DNp42 | 71 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 45µm, male-cns R 68µm | - | SS57807〈GAL4〉[特異性2] / SS57049〈GAL4〉[特異性2] |
| 外側角 | WEDPN4 → LHPV9b1 | 70 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 17µm, male-cns R 12µm | - | - |
| クランプ | MBON33 → SMP504 | 70 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 57µm, male-cns R 24µm | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] | - |
| クランプ | MBON26 → SMP156 | 68 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 15µm, male-cns R 20µm | SS80958〈GAL4〉[特異性1] | - |
| 後斜面 | IB092 → DNp49 | 67 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 35µm, male-cns R 37µm | - | - |
| クランプ / LexA 12 型を軸に再探索 | MBON01 → SMP053 | 66 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 17µm, male-cns R 27µm | SS01308〈GAL4〉[特異性2] / SL86047〈LexA〉[特異性2] | - |
| クランプ | MBON32 → LHCENT5 | 66 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 48µm, male-cns R 37µm | - | - |
| クランプ | MBON27 → ATL017 | 64 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 13µm, male-cns R 29µm | - | - |
| 側方付属葉 | LAL155 → DNa01 | 64 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 14µm, male-cns R 11µm | - | SS00731〈GAL4〉[特異性1] / SS00731〈GAL4〉[特異性0] |
| 外側角 | M_l2PNl22 → LHAV3e1 | 63 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 31µm, male-cns R 41µm | SS32193〈GAL4〉[特異性2] / SS32191〈GAL4〉[特異性2] | - |
| 外側角 | M_imPNl92 → LHAV2k8 | 63 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 36µm, male-cns R 31µm | SS01097〈GAL4〉[特異性0] / SS01092〈GAL4〉(+SLP374)[特異性1] | - |
| 上前大脳 | MBON33 → CRE102 | 63 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 53µm, male-cns R 32µm | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] | - |
| 腹外側前大脳 | CL065 → DNp06 | 58 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 18µm, male-cns R 21µm | - | SS02256〈GAL4〉[特異性2] / SS01047〈GAL4〉[特異性1] |
| 外側角 | V_l2PN → LHAV3h1 | 58 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 47µm, male-cns R 48µm | SS01336〈GAL4〉[特異性1] / SS02076〈GAL4〉[特異性1] | - |
| 後斜面 | VES013 → DNp27 | 55 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 24µm, male-cns R 28µm | - | SS01876〈GAL4〉[特異性2] / MB202B〈GAL4〉[特異性2] |
| 側方付属葉 / 楔形領域 | VES030 → DNb05 | 55 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 13µm, male-cns R 35µm | - | SS01051〈GAL4〉[特異性2] / SS01051〈GAL4〉[特異性2] |
| 外側角 | M_vPNml55 → LHPD5d1 | 54 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 20µm, male-cns R 28µm | - | - |
| 後斜面 | PS060 → DNb03 | 54 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 38µm, male-cns R 34µm | - | SS02370〈GAL4〉[特異性2] / SS02370〈GAL4〉[特異性2] |
| 側方付属葉 | VES018 → DNa02 | 54 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 20µm, male-cns R 15µm | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 視覚 | LT79 → AVLP372 | 54 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 50µm, male-cns R 44µm | - | - |
| 外側角 | M_lvPNm39 → LHAD1c3 | 54 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 49µm, male-cns R 30µm | - | - |
| 側方付属葉 | LAL163 → DNa01 | 53 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 12µm, male-cns R 16µm | - | SS00731〈GAL4〉[特異性1] / SS00731〈GAL4〉[特異性0] |
| 外側角 | M_l2PNl22 → LHPV4m1 | 53 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 17µm, male-cns R 25µm | SS32193〈GAL4〉[特異性2] / SS32191〈GAL4〉[特異性2] | - |
| 外側角 | DM4_vPN → LHAV4a2 | 53 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 13µm, male-cns R 38µm | - | - |
| 側方付属葉 | LAL142 → DNb01 | 52 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 23µm, male-cns R 17µm | - | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| 外側角 | M_vPNml65 → LHCENT4 | 52 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 28µm, male-cns R 42µm | SS10947〈GAL4〉(+M_vPNml64,M_vPNml68)[特異性2] / SS01832〈GAL4〉(+M_vPNml63,M_vPNml64)[特異性1] | - |
| 外側角 | M_lvPNm39 → LHPD5c1 | 51 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 16µm, male-cns R 30µm | - | - |
| 後斜面 / 腹外側前大脳 | PVLP022 → DNp07 | 51 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 11µm, male-cns R 35µm | - | SS02276〈GAL4〉[特異性2] / SS02612〈GAL4〉[特異性2] |
| 後斜面 | PS065 → DNb02 | 51 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 22µm, male-cns R 17µm | - | SS01060〈GAL4〉[特異性2] / SS01060〈GAL4〉[特異性2] |
| 外側角 | V_l2PN → LHAD3e1_a | 50 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 42µm, male-cns R 39µm | SS01336〈GAL4〉[特異性1] / SS02076〈GAL4〉[特異性1] | - |
| 外側角 | DM4_vPN → LHPV4a8 | 50 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 41µm, male-cns R 40µm | - | - |
| 腹外側前大脳 | CL263 → DNp10 | 49 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 43µm, male-cns R 48µm | - | SS01049〈GAL4〉[特異性2] / SS01580〈GAL4〉[特異性2] |
| 外側角 | DM3_vPN → LHAV4j1 | 49 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 25µm, male-cns R 28µm | - | - |
| 腹外側前大脳 | PLP130 → DNp27 | 49 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 29µm, male-cns R 24µm | - | SS01876〈GAL4〉[特異性2] / MB202B〈GAL4〉[特異性2] |
| 外側角 | VP1m_l2PN → LHPV4j3 | 48 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 46µm, male-cns R 24µm | - | - |
| 外側角 | M_l2PNm15 → LHPD5d1 | 48 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 17µm, male-cns R 14µm | - | - |
| 側方付属葉 | SIP022 → DNb01 | 47 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 26µm, male-cns R 12µm | - | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| 外側角 | VA1d_vPN → LHAV2g2_b | 47 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 13µm, male-cns R 28µm | SS34963〈GAL4〉[特異性1] / SS01114〈GAL4〉[特異性1] | - |
| 後斜面 | VES001 → DNp10 | 46 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 24µm, male-cns R 16µm | - | SS01049〈GAL4〉[特異性2] / SS01580〈GAL4〉[特異性2] |
| 外側角 | M_l2PNl22 → LHAV2g1 | 46 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 24µm, male-cns R 37µm | SS32193〈GAL4〉[特異性2] / SS32191〈GAL4〉[特異性2] | - |
| 外側角 | M_vPNml79 → LHAV4g14 | 45 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 21µm, male-cns R 24µm | - | - |
| 側方付属葉 | CRE004 → DNa02 | 45 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 36µm, male-cns R 18µm | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 外側角 | DC1_adPN → LHPD4d1 | 45 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 18µm, male-cns R 25µm | - | - |
| 側方付属葉 | MBON35 → DNa02 | 45 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 20µm, male-cns R 32µm | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 外側角 | VA2_adPN → LHPD5b1 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 15µm, male-cns R 26µm | - | SS24617〈GAL4〉[特異性2] |
| 楔形領域 | LAL142 → DNa02 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 33µm, male-cns R 20µm | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 後斜面 | PS137 → DNb03 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 30µm, male-cns R 20µm | - | SS02370〈GAL4〉[特異性2] / SS02370〈GAL4〉[特異性2] |
| 視覚 | LT34 → AOTU035 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 38µm, male-cns R 18µm | SS80350〈GAL4〉[特異性1] | - |
| 外側角 | M_l2PNl21 → LHPV1c2 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 27µm, male-cns R 57µm | SS32192〈GAL4〉(+M_l2PNl22)[特異性1] | - |
| 後斜面 | PS106 → DNp49 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 42µm, male-cns R 17µm | - | - |
| 外側角 | M_vPNml86 → LHAD1c3 | 43 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 33µm, male-cns R 38µm | - | - |
| 外側角 | V_l2PN → LHPV4g2 | 43 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 28µm, male-cns R 36µm | SS01336〈GAL4〉[特異性1] / SS02076〈GAL4〉[特異性1] | - |
| 外側角 | VA1d_vPN → LHAD1f2 | 42 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 30µm, male-cns R 62µm | SS34963〈GAL4〉[特異性1] / SS01114〈GAL4〉[特異性1] | - |
| 視覚 | LT75 → AVLP464 | 41 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 23µm, male-cns R 21µm | - | - |
| キノコ体 | VP2_adPN → LHCENT5 | 35 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 24µm | - | - |
| キノコ体 | M_adPNm4 → LHAV3f1 | 31 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 11µm, male-cns R 19µm | - | - |
| キノコ体 | DA4l_adPN → LHPV6g1 | 29 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 11µm, male-cns R 14µm | - | - |
| 中心複合体 | FB6L → PFL2 | 238 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 14µm, male-cns R 0µm | SS95206〈GAL4〉 / SS95750〈GAL4〉 | SS82331〈GAL4〉[特異性2] / SS82329〈GAL4〉 |
| 側方付属葉 | LAL119 → DNb01 | 218 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 5µm, male-cns R 10µm | SS50337〈GAL4〉[特異性2] | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| 中心複合体 | FB6L → PFL3 | 196 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 10µm, male-cns R 3µm | SS95206〈GAL4〉 / SS95750〈GAL4〉 | SS82335〈GAL4〉 / SS60291〈GAL4〉 |
| キノコ体の中だけ | MBON04 → MBON15-like | 107 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] | MB257C〈GAL4〉(+MBON19)[特異性1] / MB540B〈GAL4〉(+MBON18)[特異性1] |
| 側方付属葉 | AOTU042 → DNp17 | 107 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | SS85999〈GAL4〉[特異性1] | SS02553〈GAL4〉[特異性0] |
| 中心複合体 | ExR2 → FB6M | 105 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 4µm, male-cns R 49µm | SS58826〈GAL4〉[特異性2] / SS59679〈GAL4〉[特異性2] | SS49931〈GAL4〉[特異性2] / SS72731〈GAL4〉[特異性2] |
| 側方付属葉 | LAL072 → DNa01 | 97 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 7µm, male-cns R 0µm | SS49386〈GAL4〉[特異性2] | SS00731〈GAL4〉[特異性1] / SS00731〈GAL4〉[特異性0] |
| キノコ体の中だけ | MBON15-like → MBON04 | 95 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | MB257C〈GAL4〉(+MBON19)[特異性1] / MB540B〈GAL4〉(+MBON18)[特異性1] | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] |
| キノコ体の中だけ | MBON21 → MBON03 | 94 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SS46348〈GAL4〉[特異性2] / SS81353〈GAL4〉[特異性2] | SS01127〈GAL4〉[特異性2] / SS01143〈GAL4〉[特異性2] |
| キノコ体の中だけ | MBON04 → PPL105 | 93 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] | MB059B〈GAL4〉[特異性1] / MB320B〈GAL4〉(+PPL101)[特異性1] |
| キノコ体の中だけ | PPL106 → PPL103 | 92 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 10µm, male-cns R 14µm | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] | MB302B〈GAL4〉[特異性2] / MB302B〈GAL4〉[特異性0] |
| 側方付属葉 | LAL098 → DNb01 | 82 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 4µm, male-cns R 2µm | SS49882〈GAL4〉(+VES019)[特異性2] / SS49882〈GAL4〉(+VES019)[特異性1] | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| キノコ体の中だけ | MBON03 → MBON21 | 78 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SS01127〈GAL4〉[特異性2] / SS01143〈GAL4〉[特異性2] | SS46348〈GAL4〉[特異性2] / SS81353〈GAL4〉[特異性2] |
| 後斜面 | PS124 → DNb01 | 68 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 3µm, male-cns R 15µm | SS67678〈GAL4〉[特異性2] / SS66732〈GAL4〉[特異性0] | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| キノコ体の中だけ | MBON21 → PPL106 | 60 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 13µm, male-cns R 4µm | SS46348〈GAL4〉[特異性2] / SS81353〈GAL4〉[特異性2] | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] |
| キノコ体の中だけ | PPL106 → MBON17 | 55 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] | SS01298〈GAL4〉(+MBON16)[特異性2] / SS01397〈GAL4〉(+MBON16)[特異性2] |
| キノコ体の中だけ | MBON04 → MBON20 | 50 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] | SS67741〈GAL4〉[特異性2] / SS95118〈GAL4〉[特異性2] |
| キノコ体の中だけ | MBON11 → MBON33 | 41 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 9µm, male-cns R 6µm | SL86058〈LexA〉[特異性2] / SS01359〈GAL4〉[特異性1] | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] |
| キノコ体の中だけ | MBON33 → MBON11 | 38 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 12µm, male-cns R 4µm | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] | SL86058〈LexA〉[特異性2] / SS01359〈GAL4〉[特異性1] |
| 側方付属葉 | ExR6 → DNb01 | 34 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 10µm, male-cns R 5µm | SS53617〈GAL4〉[特異性2] / SS53617〈GAL4〉 | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| キノコ体の中だけ | PPL103 → MBON07 | 30 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 3µm, male-cns R 4µm | MB302B〈GAL4〉[特異性2] / MB302B〈GAL4〉[特異性0] | MB319B〈GAL4〉[特異性2] / MB319C〈GAL4〉[特異性2] |
| キノコ体の中だけ | PPL101 → MBON01 | 28 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 3µm, male-cns R 0µm | SL86067〈LexA〉[特異性2] / MB503B〈GAL4〉[特異性1] | SS01308〈GAL4〉[特異性2] / SL86047〈LexA〉[特異性2] |
| キノコ体の中だけ | MBON20 → MBON04 | 28 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SS67741〈GAL4〉[特異性2] / SS95118〈GAL4〉[特異性2] | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] |
| キノコ体の中だけ | MBON01 → PPL104 | 26 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | 未確認 | SS01308〈GAL4〉[特異性2] / SL86047〈LexA〉[特異性2] | MB568B〈GAL4〉(+PPL105)[特異性0] / MB099C〈GAL4〉(+PPL103,PPL105,PPL106)[特異性2] |
| キノコ体の中だけ | MBON33 → PPL106 | 24 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm, male-cns R 12µm | SS75199〈GAL4〉[特異性1] / SS75200〈GAL4〉[特異性1] | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] |
| キノコ体の中だけ | PPL106 → MBON21 | 22 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm, male-cns R 1µm | SL86050〈LexA〉[特異性2] / SL86062〈LexA〉[特異性2] | SS46348〈GAL4〉[特異性2] / SS81353〈GAL4〉[特異性2] |
| 側方付属葉 | VES031 → DNp17 | 369 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| クランプ | MBON05 → CRE072 | 173 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm, male-cns R 0µm | SS01188〈GAL4〉[特異性2] / SS01297〈GAL4〉[特異性1] | - |
| 側方付属葉 | SAD036 → DNp17 | 138 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 上前大脳 | MBON06 → CRE024 | 133 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 9µm, male-cns R 21µm | SS01224〈GAL4〉[特異性1] / MB482C〈GAL4〉[特異性1] | - |
| クランプ | MBON26 → SMP031 | 124 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns R 6µm | SS80958〈GAL4〉[特異性1] | - |
| 触角葉 | DA4m_adPN → ALBN1 | 116 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | nan | - | - |
| 側方付属葉 / 後斜面 | LAL194 → DNa08 | 116 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm, male-cns R 18µm | - | SS02393〈GAL4〉[特異性2] / SS02393〈GAL4〉[特異性2] |
| 後斜面 | PS248 → DNb01 | 107 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 4µm, male-cns R 32µm | - | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| 側方付属葉 | ER6 → LAL152 | 104 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 8µm, male-cns R 14µm | SS58833〈GAL4〉[特異性2] / SS58833〈GAL4〉 | - |
| 側方付属葉 | VES003 → DNp17 | 99 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| キノコ体 | PAM13 → MBON32 | 91 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 13µm, male-cns R 8µm | MB023B〈GAL4〉[特異性1] / MB055B〈GAL4〉(+PAM14)[特異性2] | - |
| 側方付属葉 / 後斜面 | PS013 → DNp17 | 91 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm, male-cns R 4µm | - | SS02553〈GAL4〉[特異性0] |
| 側方付属葉 | FB1C → LAL048 | 85 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | SS63879〈GAL4〉[特異性2] / SS63927〈GAL4〉[特異性2] | - |
| 外側角 | M_vPNml86 → LHPV1c1 | 81 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 8µm, male-cns R 70µm | - | - |
| 外側角 | DA4m_adPN → LHPV4a11 | 75 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 6µm, male-cns R 10µm | - | - |
| クランプ / LexA 12 型を軸に再探索 | MBON04 → LAL175 | 74 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 6µm, male-cns R 11µm | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] | - |
| 側方付属葉 / 後斜面 | PS049 → DNp17 | 74 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 後斜面 | PS137 → DNa01 | 71 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 6µm, male-cns R 4µm | - | SS00731〈GAL4〉[特異性1] / SS00731〈GAL4〉[特異性0] |
| 側方付属葉 | VES034_b → DNb05 | 71 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 9µm, male-cns R 8µm | - | SS01051〈GAL4〉[特異性2] / SS01051〈GAL4〉[特異性2] |
| キノコ体 | PPL102 → MBON24 | 70 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 3µm, male-cns R 3µm | SS01252〈GAL4〉(+MBON03,MBON04)[特異性0] | - |
| 側方付属葉 / 後斜面 | VES056 → DNp17 | 69 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 1µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 外側角 | M_vPNml86 → LHAV2m1 | 66 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 4µm, male-cns R 19µm | - | - |
| 側方付属葉 | VES033 → DNp17 | 63 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| キノコ体 | KCab-p → MBON13 | 62 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 1µm, male-cns R 4µm | - | MB019B〈GAL4〉[特異性2] / MB037B〈GAL4〉[特異性2] |
| 側方付属葉 / 楔形領域 | VES025 → DNb05 | 62 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 9µm, male-cns R 13µm | - | SS01051〈GAL4〉[特異性2] / SS01051〈GAL4〉[特異性2] |
| 側方付属葉 | LAL110 → DNa01 | 60 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 7µm, male-cns R 3µm | - | SS00731〈GAL4〉[特異性1] / SS00731〈GAL4〉[特異性0] |
| キノコ体 / 外側角 | DA4m_adPN → LHPV9b1 | 57 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 14µm, male-cns R 7µm | - | - |
| 側方付属葉 | VES011 → DNp17 | 56 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 側方付属葉 | VES018 → DNp17 | 56 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 後斜面 | PS126 → DNa09 | 56 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | - |
| 外側角 | M_vPNml86 → LHPV6c1 | 53 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 7µm, male-cns R 30µm | - | LH2220〈GAL4〉[特異性1] |
| 視覚 | LT73 → LHAV2b3 | 51 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 10µm, male-cns R 29µm | SS80378〈GAL4〉[特異性1] | - |
| 側方付属葉 | PS185 → DNp17 | 51 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 側方付属葉 | LAL021 → DNp17 | 50 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 腹外側前大脳 | PVLP121 → DNp06 | 50 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 4µm, male-cns R 7µm | - | SS02256〈GAL4〉[特異性2] / SS01047〈GAL4〉[特異性1] |
| 側方付属葉 | VES050 → DNp17 | 48 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 後斜面 | PS127 → DNp17 | 48 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 側方付属葉 | VES047 → DNa02 | 48 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 25µm, male-cns R 4µm | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 後斜面 | PS126 → DNa04 | 47 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02384〈GAL4〉(+DNa10)[特異性1] |
| 外側角 | DA1_vPN → LHAV2k8 | 46 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 15µm, male-cns R 10µm | - | - |
| 外側角 | DA4m_adPN → LHPV2d1 | 45 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 14µm, male-cns R 6µm | - | LH1925〈GAL4〉[特異性2] / LH1899〈GAL4〉[特異性2] |
| 後斜面 | PS010 → DNp17 | 45 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 後斜面 | VS → PS262 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | SS04438〈GAL4〉[特異性2] | - |
| 側方付属葉 / 後斜面 | VES064 → DNp17 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 側方付属葉 | VES070 → DNa01 | 44 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | - | SS00731〈GAL4〉[特異性1] / SS00731〈GAL4〉[特異性0] |
| 外側角 | VA1v_vPN → LHAV2k12_b | 43 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 5µm, male-cns R 4µm | - | - |
| 外側角 | VL1_vPN → LHAD2d1 | 43 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 6µm, male-cns R 15µm | SS86025〈GAL4〉[特異性2] / SS34947〈GAL4〉[特異性2] | - |
| 側方付属葉 | PS214 → DNp17 | 42 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 11µm, male-cns R 0µm | - | SS02553〈GAL4〉[特異性0] |
| 楔形領域 | CRE074 → DNb05 | 37 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 6µm, male-cns R 17µm | - | SS01051〈GAL4〉[特異性2] / SS01051〈GAL4〉[特異性2] |
| キノコ体 | M_adPNm4 → LHPV10b1 | 35 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 1µm | - | - |
| キノコ体 | V_ilPN → LHPV6c1 | 32 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 14µm, male-cns R 9µm | - | LH2220〈GAL4〉[特異性1] |
| 後斜面 | VS → PS239 | 30 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 0µm, male-cns R 0µm | SS04438〈GAL4〉[特異性2] | - |
| 腹外側前大脳 | LT62 → PVLP094 | 27 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 15µm, male-cns R 3µm | - | - |
| 腹外側前大脳 | LT83 → DNp06 | 21 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 16µm, male-cns R 9µm | - | SS02256〈GAL4〉[特異性2] / SS01047〈GAL4〉[特異性1] |
| 後斜面 | HSS → PS234 | 16 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 7µm, male-cns R 4µm | SS22522〈GAL4〉(+HSE,HSN)[特異性2] | - |
| キノコ体 | VP1m_l2PN → KCg-s2 | 15 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm, male-cns R 1µm | - | - |
| 後斜面 | HSS → PS233 | 14 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 6µm, male-cns R 4µm | SS22522〈GAL4〉(+HSE,HSN)[特異性2] | - |
| キノコ体 | VP1d_il2PN → KCg-s2 | 11 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm, male-cns R 4µm | SS01104〈GAL4〉[特異性1] / SS01103〈GAL4〉[特異性1] | - |
| キノコ体 | VP2_adPN → KCg-s1 | 11 µm | 0 | flywire:783(合計 0/0), male-cns:v1.0(合計 0/0) | male-cns L 2µm | - | - |
| 上前大脳 / クランプ | PPL102 → CL138 | 218 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01252〈GAL4〉(+MBON03,MBON04)[特異性0] | SS88586〈GAL4〉[特異性2] |
| 側方付属葉 | FB5V → AOTU042 | 184 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS88825〈GAL4〉 / SS88827〈GAL4〉 | SS85999〈GAL4〉[特異性1] |
| 外側角 | VP1m+_lvPN → MB-C1 | 96 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | MB380C〈GAL4〉[特異性2] / LH2449〈GAL4〉[特異性2] |
| 後斜面 | CL033 → DNp31 | 68 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS47816〈GAL4〉[特異性2] | SS95362〈GAL4〉[特異性2] / SS95510〈GAL4〉[特異性2] |
| 側方付属葉 | CL033 → DNp17 | 63 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS47816〈GAL4〉[特異性2] | SS02553〈GAL4〉[特異性0] |
| 後斜面 | PS263 → DNp16 | 56 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS35837〈GAL4〉[特異性1] | SS00735〈GAL4〉[特異性1] / SS00735〈GAL4〉 |
| 後斜面 | PLP030 → DNp15 | 53 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS32832〈GAL4〉[特異性1] | SS01556〈GAL4〉[特異性1] / SS01556〈GAL4〉[特異性1] |
| 楔形領域 / 後斜面 / 腹外側前大脳 | PLP030 → DNp31 | 49 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS32832〈GAL4〉[特異性1] | SS95362〈GAL4〉[特異性2] / SS95510〈GAL4〉[特異性2] |
| 側方付属葉 | CL033 → DNp13 | 47 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS47816〈GAL4〉[特異性2] | SS61089〈GAL4〉[特異性2] / SS61090〈GAL4〉[特異性2] |
| 後斜面 | CL033 → DNb05 | 47 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS47816〈GAL4〉[特異性2] | SS01051〈GAL4〉[特異性2] / SS01051〈GAL4〉[特異性2] |
| 側方付属葉 | CL033 → DNa08 | 46 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS47816〈GAL4〉[特異性2] | SS02393〈GAL4〉[特異性2] / SS02393〈GAL4〉[特異性2] |
| 側方付属葉 / 後斜面 | PS196 → DNa02 | 44 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS48753〈GAL4〉[特異性2] | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 楔形領域 | OA-VUMa1 → DNp19 | 36 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS48831〈GAL4〉[特異性0] / SS48816〈GAL4〉(+OA-VUMa7)[特異性0] | SS81965〈GAL4〉[特異性2] / SS12653〈GAL4〉(+DNp15)[特異性1] |
| 触角葉 | VL1_ilPN → M_lv2PN9t49 | 224 µm | 0 | flywire:783(合計 0/0) | nan | MB546B〈GAL4〉[特異性2] / MB567B〈GAL4〉[特異性2] | - |
| 触角葉 | v2LN36 → SIFa | 191 µm | 0 | flywire:783(合計 0/0) | nan | - | MB015B〈GAL4〉[特異性2] |
| 後斜面 | LAL127 → PS039 | 166 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 上前大脳 | LHPV5e1 → SMP399 | 157 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 側方付属葉 / 後斜面 | PS271 → DNa02 | 155 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 後斜面 / 腹外側前大脳 | PLP195 → DNp07 | 151 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS02276〈GAL4〉[特異性2] / SS02612〈GAL4〉[特異性2] |
| 触角葉 | M_lvPNm41 → M_lv2PN9t49 | 147 µm | 0 | flywire:783(合計 0/0) | nan | - | - |
| 上前大脳 | LHPD5b1 → SMP585 | 145 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS24617〈GAL4〉[特異性2] | - |
| 触角葉 | VP5+VP3_l2PN → M_l2PNm14 | 142 µm | 0 | flywire:783(合計 0/0) | nan | SS32229〈GAL4〉[特異性2] / MB559C〈GAL4〉(+VP5+_l2PN,VP5+VP2_l2PN)[特異性0] | - |
| 後斜面 | LAL081 → VES026 | 140 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 | LAL099 → PS015 | 136 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 触角葉 | M_vPNml64 → M_lv2PN9t49 | 133 µm | 0 | flywire:783(合計 0/0) | nan | SS01823〈GAL4〉[特異性1] / SS10947〈GAL4〉(+M_vPNml65,M_vPNml68)[特異性2] | - |
| 触角葉 | VP5+VP2_l2PN → M_l2PNl23 | 132 µm | 0 | flywire:783(合計 0/0) | nan | MB559C〈GAL4〉(+VP5+_l2PN,VP5+VP3_l2PN)[特異性0] | - |
| 上前大脳 / クランプ | MBON32 → IB112 | 123 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS81529〈GAL4〉[特異性1] |
| 上前大脳 | LHAV1e1 → SMP029 | 120 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS35966〈GAL4〉[特異性2] |
| 後斜面 | PS275 → DNb01 | 116 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS02383〈GAL4〉[特異性2] / SS02542〈GAL4〉[特異性1] |
| 上前大脳 | LHCENT9 → pC1b | 113 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 | LAL141 → PS212 | 107 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 / クランプ | PS293 → DNp31 | 105 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS95362〈GAL4〉[特異性2] / SS95510〈GAL4〉[特異性2] |
| 後斜面 | LAL099 → PS064 | 102 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 側方付属葉 / 後斜面 | VES042 → DNp17 | 101 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS02553〈GAL4〉[特異性0] |
| 後斜面 | PS136 → DNb02 | 98 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS01060〈GAL4〉[特異性2] / SS01060〈GAL4〉[特異性2] |
| 腹外側前大脳 | LHAV4i2 → DNp32 | 97 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS02394〈GAL4〉[特異性1] / SS02394〈GAL4〉[特異性0] |
| 後斜面 | LAL099 → VES026 | 93 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 上前大脳 / LexA 12 型を軸に再探索 | MBON04 → PAL02 | 92 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS67728〈GAL4〉[特異性2] / SS67755〈GAL4〉[特異性2] | - |
| 側方付属葉 / 後斜面 | PS196 → DNa03 | 92 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS48753〈GAL4〉[特異性2] | - |
| 外側角 | VC3m_lvPN → LHAV2m1 | 89 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 外側角 | DM3_adPN → LHPD3a3_c | 85 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 | LAL046 → PS064 | 85 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 | LAL010 → PS163 | 82 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| キノコ体 / 外側角 | VP1m+_lvPN → LHCENT8 | 78 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | - |
| 外側角 | VP1m+_lvPN → LHPV4a11 | 78 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | - |
| 後斜面 | PS054 → DNp16 | 78 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS00735〈GAL4〉[特異性1] / SS00735〈GAL4〉 |
| 腹外側前大脳 | PLP070 → DNp27 | 78 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS01876〈GAL4〉[特異性2] / MB202B〈GAL4〉[特異性2] |
| 外側角 | VP1m+_lvPN → LHPV1c2 | 77 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | - |
| 外側角 | VP1d+VP4_l2PN1 → LHPV4a2 | 76 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | LH542〈GAL4〉(+LHPV4a1,LHPV4a3,LHPV4a4,LHPV4a5)[特異性2] / LH1475〈GAL4〉(+LHPV4a1,LHPV4a3,LHPV4a4,LHPV4a5)[特異性2] |
| 外側角 | VP1m+_lvPN → LHAV2m1 | 76 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | - |
| 後斜面 | 5-HTPMPV03 → DNp16/17 | 76 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| キノコ体 / 外側角 | VP1m+_lvPN → LHPV6g1 | 74 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | - |
| 後斜面 / 腹外側前大脳 | IB006 → DNb05 | 71 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS01051〈GAL4〉[特異性2] / SS01051〈GAL4〉[特異性2] |
| キノコ体 | VP1m+VP5_ilPN → LHPD5a1 | 68 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS32190〈GAL4〉[特異性0] / SS01100〈GAL4〉(+VP1d_il2PN)[特異性2] | - |
| 後斜面 | PS213 → DNp16 | 58 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS00735〈GAL4〉[特異性1] / SS00735〈GAL4〉 |
| 後斜面 | CL033 → DNp49 | 58 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS47816〈GAL4〉[特異性2] | - |
| 後斜面 | PS117 → DNp19 | 54 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS81965〈GAL4〉[特異性2] / SS12653〈GAL4〉(+DNp15)[特異性1] |
| 外側角 | M_vPNml64 → LHPV3b1_b | 52 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01823〈GAL4〉[特異性1] / SS10947〈GAL4〉(+M_vPNml65,M_vPNml68)[特異性2] | - |
| 外側角 | M_vPNml58 → LHAV4a5 | 52 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01202〈GAL4〉(+M_vPNml56)[特異性2] / SS01202〈GAL4〉(+M_vPNml56)[特異性2] | - |
| 後斜面 | PS239 → DNp16 | 51 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS00735〈GAL4〉[特異性1] / SS00735〈GAL4〉 |
| 腹外側前大脳 | PLP146 → DNp27 | 51 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS01876〈GAL4〉[特異性2] / MB202B〈GAL4〉[特異性2] |
| 外側角 | M_vPNml50 → LHAD1b2_a | 50 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 | PS018 → DNp07 | 50 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS02276〈GAL4〉[特異性2] / SS02612〈GAL4〉[特異性2] |
| 後斜面 | IB029 → DNp16/17 | 49 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 外側角 | VP1d+VP4_l2PN1 → LHAV3q1 | 49 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 | PS277 → DNp19 | 49 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS81965〈GAL4〉[特異性2] / SS12653〈GAL4〉(+DNp15)[特異性1] |
| 外側角 | M_vPNml85 → LHAV6a8 | 49 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01108〈GAL4〉(+VA1d_vPN,M_vPNml88,M_vPNml89)[特異性2] / SS01108〈GAL4〉(+VA1d_vPN,M_vPNml88,M_vPNml89)[特異性1] | - |
| 後斜面 | PLP242 → DNa09 | 48 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 外側角 | M_vPNml58 → LHCENT1 | 47 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01202〈GAL4〉(+M_vPNml56)[特異性2] / SS01202〈GAL4〉(+M_vPNml56)[特異性2] | - |
| 後斜面 | IB039 → DNp49 | 47 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 外側角 | VA1d_vPN → LHAD1b2_c | 47 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS34963〈GAL4〉[特異性1] / SS01114〈GAL4〉[特異性1] | - |
| 外側角 | M_vPNml85 → LHAD1a3 | 45 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01108〈GAL4〉(+VA1d_vPN,M_vPNml88,M_vPNml89)[特異性2] / SS01108〈GAL4〉(+VA1d_vPN,M_vPNml88,M_vPNml89)[特異性1] | - |
| 外側角 | VP1m+_lvPN → LHPV4b5 | 45 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | - |
| 後斜面 | IB006 → DNp26 | 45 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS01553〈GAL4〉[特異性2] / SS01553〈GAL4〉[特異性2] |
| 後斜面 | PLP195 → DNa02 | 45 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS00730〈GAL4〉[特異性1] / SS00730〈GAL4〉[特異性0] |
| 外側角 | V_l2PN → LHAV3e3 | 44 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01336〈GAL4〉[特異性1] / SS02076〈GAL4〉[特異性1] | - |
| 外側角 | DP1m_adPN → LHAV2k3 | 44 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 後斜面 | PS212 → DNp16/17 | 43 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 外側角 | M_l2PNl22 → LHAD1b2_a | 42 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS32193〈GAL4〉[特異性2] / SS32191〈GAL4〉[特異性2] | - |
| キノコ体 | VC3m_lvPN → LHPV10b1 | 41 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 視覚 | LT34 → LAL026 | 40 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS80350〈GAL4〉[特異性1] | - |
| 触角葉 | ORN_V → TRN_VP2 | 39 µm | 0 | flywire:783(合計 0/0) | nan | SS00525〈GAL4〉[特異性2] / SS00525〈GAL4〉[特異性1] | - |
| 触角葉 | ORN_VC3l → M_vPNml87 | 39 µm | 0 | flywire:783(合計 0/0) | nan | - | SS00587〈GAL4〉[特異性1] / SS00586〈GAL4〉[特異性1] |
| 触角葉 | ORN_VM5v → M_vPNml77 | 38 µm | 0 | flywire:783(合計 0/0) | nan | - | - |
| 楔形領域 | SAD016 → DNp11 | 36 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | SS49051〈GAL4〉[特異性2] / SS02891〈GAL4〉[特異性2] |
| キノコ体 | VP2+_adPN → LHPV8d1 | 32 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 触角葉 | ORN_VA5 → M_vPNml50 | 31 µm | 0 | flywire:783(合計 0/0) | nan | - | - |
| 触角葉 | ORN_VA6 → M_vPNml50 | 31 µm | 0 | flywire:783(合計 0/0) | nan | - | - |
| キノコ体 | VC3l_adPN → LHPV10b1 | 29 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS01119〈GAL4〉[特異性1] / SS34061〈GAL4〉[特異性1] | - |
| 触角葉 | ORN_VC3l → M_adPNm3 | 27 µm | 0 | flywire:783(合計 0/0) | nan | - | - |
| 腹外側前大脳 | LT36 → PS155 | 27 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 触角葉 | ORN_VM7d → M_vPNml50 | 24 µm | 0 | flywire:783(合計 0/0) | nan | - | - |
| 触角葉 | ORN_VA7m → VA2_adPN | 23 µm | 0 | flywire:783(合計 0/0) | nan | - | - |
| キノコ体 | VP2+_adPN → KCg-s1 | 23 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| 触角葉 | ORN_VM5v → M_imPNl92 | 20 µm | 0 | flywire:783(合計 0/0) | nan | - | SS01097〈GAL4〉[特異性0] / SS01092〈GAL4〉(+SLP374)[特異性1] |
| キノコ体 | VP1m+VP5_ilPN → KCg-s1 | 11 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS32190〈GAL4〉[特異性0] / SS01100〈GAL4〉(+VP1d_il2PN)[特異性2] | - |
| キノコ体 | VP1d+VP4_l2PN1 → KCg-s3 | 11 µm | 0 | flywire:783(合計 0/0) | 未確認 | - | - |
| キノコ体 | VP1m+_lvPN → KCg-d | 10 µm | 0 | flywire:783(合計 0/0) | 未確認 | SS33915〈GAL4〉[特異性2] | - |
