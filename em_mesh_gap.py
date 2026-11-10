"""MBON04 と FB4N のメッシュ（実際の膜）同士の距離を測る。

骨格ベースの「接触 1 µm 以内」では、隙間にグリアが入れるかどうかは判断できない。
メッシュ頂点どうしの最短距離を見れば、膜が本当に向かい合っているかが分かる。
"""
import sys
sys.path.insert(0, r"C:\Users\morit\Programs-repo\connectomics")
import numpy as np
from common.loaders import neuprint_client, _retry
from neuprint import NeuronCriteria as NC
from navis.interfaces import neuprint as nvn
from scipy.spatial import cKDTree

VOX = 8.0  # nm / voxel
名 = {487925063: "MBON04_R", 613079053: "MBON04_L", 673746277: "FB4N_R"}
# 比較のため、既に結合している組も測る: MBON26_R -> FB4N_R, MBON04 -> FB4C
名.update({895441451: "MBON26_R"})

c = neuprint_client("hemibrain:v1.2.1")
meshes = {}
for b in 名:
    m = _retry(nvn.fetch_mesh_neuron, NC(bodyId=b), client=c, lod=2)
    m = m[0] if hasattr(m, "__len__") else m
    v = np.asarray(m.vertices, float) * VOX
    meshes[b] = v
    print(f"{名[b]:9s} 頂点 {len(v):7d}")

def gap(a, b, label):
    va, vb = meshes[a], meshes[b]
    t = cKDTree(vb)
    d, _ = t.query(va, k=1)
    close = d[d < 2000]  # 2 µm 以内の部分だけを「向かい合っている面」として見る
    if not len(close):
        print(f"\n{label}: 2 µm 以内の面が無い")
        return
    q = np.percentile(close, [1, 5, 25, 50])
    print(f"\n{label}")
    print(f"  最短 {d.min():7.0f} nm   50 nm 以内の頂点 {(d<50).sum():5d}   "
          f"100 nm 以内 {(d<100).sum():5d}   200 nm 以内 {(d<200).sum():5d}   500 nm 以内 {(d<500).sum():6d}")
    print(f"  向かい合う面の距離 1%={q[0]:.0f} 5%={q[1]:.0f} 25%={q[2]:.0f} 中央={q[3]:.0f} nm")

gap(613079053, 673746277, "MBON04_L ↔ FB4N_R（調べたい組・シナプス 0）")
gap(487925063, 673746277, "MBON04_R ↔ FB4N_R（調べたい組・シナプス 0）")
gap(895441451, 673746277, "MBON26_R ↔ FB4N_R（既に 24 シナプスある組＝対照）")
