"""正解つき合成データ（アプローチ同士を定量比較するためのベンチマーク）。

生成モデル:
  1. 細胞タイプ: 同じタイプのニューロンは、樹状突起・軸索アーバーの中心を共有する（±3 µm のずれ）
  2. 各ニューロン = 細胞体 → 一次神経突起 → 樹状突起アーバー、そこから軸索幹 → 軸索アーバー
     （アーバーは 40 µm 立方の「ニューロピル」内のランダム分岐ウォーク）
  3. 軸索(A)–樹状突起(B) の接触ケーブル長 L_AB を計算し、Peters' rule に従って
         シナプス数 ~ Poisson(rate × L_AB × g_type(T_A,T_B) × g_AB)
     g_type はタイプの組ごとの「シナプスの作りやすさ」（ばらつき大）、g_AB は個体差（小）
  4. 正解（ペア単位）: L_AB が十分大きいペアの一部のシナプス数を 0〜weak_synapse_max に落とす
  5. 正解（パッチ単位, アプローチ F 用）: 普通に結合しているペアの中で、特定の接触パッチだけ
     シナプスを置かない（「同じペアの中でもシナプスを避ける部位」）

注意: 仕込み以外にも、偶然 Poisson で 0 になる近接ペアは存在する（現実と同じく、それは正解に含めない）。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import Thresholds
from .geometry import aggregate_contacts, cable_points, close_point_pairs, contact_patches, stack
from .model import AXON, DENDRITE, UNKNOWN, Dataset, Neuron

UM = 1000.0


def _arbor(rng, root, center, n_nodes, step, spread, branch_p):
    nodes, parents = [root], [-1]
    tips = [(0, rng.normal(size=3))]
    while len(nodes) < n_nodes:
        t = rng.integers(len(tips))
        idx, d = tips[t]
        p = nodes[idx]
        pull = (center - p) / spread
        nd = 0.75 * d / (np.linalg.norm(d) + 1e-9) + 0.55 * rng.normal(size=3) + 0.35 * pull
        nd /= np.linalg.norm(nd) + 1e-9
        nodes.append(p + step * nd)
        parents.append(idx)
        tips[t] = (len(nodes) - 1, nd)
        if rng.random() < branch_p:
            tips.append((len(nodes) - 1, rng.normal(size=3)))
    return np.array(nodes), np.array(parents)


def _line(a, b, step):
    k = max(1, int(np.ceil(np.linalg.norm(b - a) / step)))
    return a + (b - a) * (np.arange(1, k + 1) / k)[:, None]


def make_neuron(rng, nid, dc, ac, box_um=40.0, cable_um=250.0, step_um=1.0) -> Neuron:
    """dc, ac: 樹状突起 / 軸索アーバーの中心 [nm]"""
    box, step = box_um * UM, step_um * UM
    soma = np.r_[rng.uniform(0, box, 2), -20 * UM]
    n_arb = int(cable_um / step_um)

    nodes, parents, comp, rad = [soma[None]], [np.array([-1])], [np.array([UNKNOWN])], [np.array([1500.0])]
    # 一次神経突起
    pn = _line(soma, dc, step)
    nodes.append(pn); parents.append(np.arange(len(pn)))  # 各ノードの親 = 直前のノード（0 = 細胞体）
    comp.append(np.full(len(pn), UNKNOWN)); rad.append(np.full(len(pn), 400.0))
    off = 1 + len(pn)
    root_d = off - 1
    # 樹状突起
    dn, dp = _arbor(rng, dc, dc, n_arb, step, 8 * UM, 0.08)
    dn, dp = dn[1:], dp[1:]
    dp = np.where(dp == 0, root_d, dp - 1 + off)
    nodes.append(dn); parents.append(dp); comp.append(np.full(len(dn), DENDRITE))
    rad.append(rng.uniform(150, 350, len(dn)))
    off += len(dn)
    # 軸索幹 + 軸索
    tr = _line(dc, ac, step)
    nodes.append(tr); parents.append(np.r_[root_d, off + np.arange(len(tr) - 1)])
    comp.append(np.full(len(tr), AXON)); rad.append(np.full(len(tr), 250.0))
    off += len(tr)
    root_a = off - 1
    an, ap = _arbor(rng, ac, ac, n_arb, step, 8 * UM, 0.08)
    an, ap = an[1:], ap[1:]
    ap = np.where(ap == 0, root_a, ap - 1 + off)
    nodes.append(an); parents.append(ap); comp.append(np.full(len(an), AXON))
    rad.append(rng.uniform(100, 250, len(an)))

    return Neuron(id=nid, nodes=np.vstack(nodes), parents=np.concatenate(parents).astype(int),
                  radius=np.concatenate(rad), compartment=np.concatenate(comp).astype(np.int8))


def make_dataset(n_neurons=80, seed=0, th: Thresholds | None = None, rate_per_um=0.25,
                 neurons_per_type=5, type_sd=1.0, noise_sd=0.3, jitter_um=3.0, box_um=40.0,
                 planted_frac=0.15, planted_min_contact_um=10.0,
                 patch_link_um=2.0, patch_frac=0.5, patch_min_contact_um=5.0) -> Dataset:
    th = th or Thresholds()
    rng = np.random.default_rng(seed)
    ids = [10_000 + k for k in range(n_neurons)]

    # --- 細胞タイプとニューロン ---
    n_types = int(np.ceil(n_neurons / neurons_per_type))
    t_dc = rng.uniform(0, box_um * UM, (n_types, 3))
    t_ac = rng.uniform(0, box_um * UM, (n_types, 3))
    ntype = {i: k // neurons_per_type for k, i in enumerate(ids)}
    neurons = {}
    for i in ids:
        t = ntype[i]
        dc = t_dc[t] + rng.normal(0, jitter_um * UM, 3)
        ac = t_ac[t] + rng.normal(0, jitter_um * UM, 3)
        neurons[i] = make_neuron(rng, i, dc, ac, box_um=box_um)
        neurons[i].meta["type"] = f"T{t:02d}"
    g_type = rng.lognormal(0, type_sd, (n_types, n_types))

    # --- 正解の軸索/樹状突起で接触を計算 ---
    cps = {i: cable_points(n, th.cable_spacing_um * UM) for i, n in neurons.items()}
    cloud = stack(cps)
    pi, pj = close_point_pairs(cloud["points"], cloud["label"], th.skeleton_contact_um * UM,
                               roles=cloud["compartment"])
    con = aggregate_contacts(pi, pj, cloud, "pre", "post")
    con["L_um"] = (con.contact_a + con.contact_b) / 2 / UM
    gt = np.array([g_type[ntype[a], ntype[b]] for a, b in zip(con.pre, con.post)])
    con["g"] = gt * rng.lognormal(0, noise_sd, len(con))
    w = rng.poisson(rate_per_um * con.L_um * con.g)

    elig = np.flatnonzero(con.L_um >= planted_min_contact_um)
    planted = rng.choice(elig, size=int(round(planted_frac * len(elig))), replace=False)
    w[planted] = rng.integers(0, th.weak_synapse_max + 1, len(planted)) * (rng.random(len(planted)) < 0.4)
    con["weight"] = w
    con["planted"] = False
    con.loc[planted, "planted"] = True

    # --- 接触パッチと「シナプスを避けるパッチ」の仕込み ---
    patch = contact_patches(pi, pj, cloud, patch_link_um * UM)
    pairs = pd.DataFrame({"pre": cloud["label"][pi], "post": cloud["label"][pj], "i": pi, "j": pj,
                          "patch": patch})
    pw = cloud["weight"]
    plen = (pairs.drop_duplicates(["patch", "i"]).assign(w=lambda d: pw[d.i]).groupby("patch").w.sum()
            + pairs.drop_duplicates(["patch", "j"]).assign(w=lambda d: pw[d.j]).groupby("patch").w.sum()) / 2 / UM
    wmap = dict(zip(zip(con.pre, con.post), con.weight))
    avoided = set()
    patch_truth = []
    for (a, b), grp in pairs.groupby(["pre", "post"]):
        if wmap.get((a, b), 0) <= th.weak_synapse_max:
            continue
        ps = grp.patch.unique()
        big = [p for p in ps if plen[p] >= patch_min_contact_um]
        if len(ps) >= 2 and big and rng.random() < patch_frac:
            p = rng.choice(big)
            # 避けた後も、残りのパッチにシナプスを置ける余地があること
            if plen[ps].sum() - plen[p] >= patch_min_contact_um:
                avoided.add(p)
                m = cloud["points"][grp.i[grp.patch == p].values]
                patch_truth.append(dict(pre=a, post=b, patch_contact_um=plen[p],
                                        cx=m[:, 0].mean(), cy=m[:, 1].mean(), cz=m[:, 2].mean()))

    # --- シナプスの位置: 接触している点ペアからランダムに（避けるパッチは除く） ---
    pre_pts = {i: [] for i in ids}
    post_pts = {i: [] for i in ids}
    syn_rows = []
    for (a, b), grp in pairs.groupby(["pre", "post"]):
        k = wmap.get((a, b), 0)
        if k == 0:
            continue
        cand = grp[~grp.patch.isin(avoided)]
        sel = cand.sample(n=k, replace=len(cand) < k, random_state=int(rng.integers(1 << 31)))
        P = cloud["points"][sel.i.values] + rng.normal(0, 100, (k, 3))
        Q = cloud["points"][sel.j.values] + rng.normal(0, 100, (k, 3))
        pre_pts[a].append(P)
        post_pts[b].append(Q)
        mid = (P + Q) / 2
        syn_rows.append(pd.DataFrame({"pre": a, "post": b, "x": mid[:, 0], "y": mid[:, 1], "z": mid[:, 2]}))
    for i in ids:
        neurons[i].pre = np.vstack(pre_pts[i]) if pre_pts[i] else np.empty((0, 3))
        neurons[i].post = np.vstack(post_pts[i]) if post_pts[i] else np.empty((0, 3))

    edges = con.loc[con.weight > 0, ["pre", "post", "weight"]].reset_index(drop=True)
    truth = con.loc[con.planted, ["pre", "post", "L_um", "weight"]].reset_index(drop=True)
    return Dataset(neurons=neurons, edges=edges, name=f"synthetic_n{n_neurons}_s{seed}", truth=truth,
                   synapses=pd.concat(syn_rows, ignore_index=True),
                   meta=dict(rate_per_um=rate_per_um, contacts=con, patch_truth=pd.DataFrame(patch_truth)))
