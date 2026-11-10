"""実データの読み込み: neuPrint (hemibrain / MANC など) と FlyWire (fafbseg)。

どちらも common.model.Dataset に変換して返す（座標は nm に統一）。
認証トークンはコードに書かず、環境変数から読む:
  neuPrint : NEUPRINT_APPLICATION_CREDENTIALS （https://neuprint.janelia.org の Account → Auth Token）
  FlyWire  : fafbseg.flywire.set_chunkedgraph_secret(...) で一度保存しておく
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .model import Dataset, Neuron


# ----------------------------------------------------------------------
# neuPrint (hemibrain)
# ----------------------------------------------------------------------
def neuprint_token() -> str | None:
    """環境変数 NEUPRINT_APPLICATION_CREDENTIALS を読む。

    Windows で「ユーザー環境変数」に設定した直後は、すでに起動しているプロセスには反映されないので、
    レジストリ (HKEY_CURRENT_USER の Environment キー) からも読む。
    """
    import os
    tok = os.environ.get("NEUPRINT_APPLICATION_CREDENTIALS")
    if tok:
        return tok
    try:
        import winreg
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as k:
            return winreg.QueryValueEx(k, "NEUPRINT_APPLICATION_CREDENTIALS")[0]
    except (ImportError, OSError):
        return None


def neuprint_client(dataset: str = "hemibrain:v1.2.1", server: str = "neuprint.janelia.org"):
    from neuprint import Client
    tok = neuprint_token()
    if not tok:
        raise RuntimeError("neuPrint のトークンが見つかりません。環境変数 NEUPRINT_APPLICATION_CREDENTIALS を設定してください")
    return Client(server, dataset=dataset, token=tok)


def _retry(fn, *args, **kwargs):
    """neuPrint は数千個を問い合わせる間に接続を切ってくることがあるので、待って再試行する。"""
    import time
    last = None
    for wait in (0, 3, 10, 30, 60, 120):
        if wait:
            time.sleep(wait)
        try:
            return fn(*args, **kwargs)
        except Exception as e:  # 接続エラー・504 など
            last = e
    raise RuntimeError(f"neuPrint への問い合わせ {getattr(fn, '__name__', fn)} が失敗しました: {last}")


def load_neuprint(body_ids, dataset: str = "hemibrain:v1.2.1",
                  server: str = "neuprint.janelia.org", with_meshes: bool = False,
                  voxel_nm: float = 8.0, min_confidence: float = 0.5) -> Dataset:
    from neuprint import NeuronCriteria as NC, SynapseCriteria as SC
    from neuprint import fetch_adjacencies, fetch_neurons, fetch_synapse_connections, fetch_synapses

    client = neuprint_client(dataset, server)
    body_ids = [int(b) for b in body_ids]

    meta_df, _ = _retry(fetch_neurons, NC(bodyId=body_ids), client=client)
    types = dict(zip(meta_df.bodyId, meta_df.type))
    # instance は左右などを含む名前（例: "..._R" / "..._L"）。E を左右別に行うときに使う
    instances = dict(zip(meta_df.bodyId, meta_df.instance)) if "instance" in meta_df else {}
    extra = {c: dict(zip(meta_df.bodyId, meta_df[c])) for c in ("status", "cropped") if c in meta_df}

    syn = _retry(fetch_synapses, NC(bodyId=body_ids), SC(confidence=min_confidence), client=client)

    neurons = {}
    # スケルトンは 1 個ずつしか取れないので、並列に取得する。
    # 数百〜数千個を取る間に一時的な接続エラーが起きることがあるので、待って再試行する。
    from concurrent.futures import ThreadPoolExecutor

    def one(b):
        return _retry(client.fetch_skeleton, b, format="pandas", heal=True)

    with ThreadPoolExecutor(max_workers=8) as ex:
        skels = dict(zip(body_ids, ex.map(one, body_ids)))
    for bid in body_ids:
        sk = skels[bid]
        row_to_idx = {r: k for k, r in enumerate(sk.rowId)}
        parents = np.array([row_to_idx.get(l, -1) for l in sk.link], dtype=int)
        s = syn[syn.bodyId == bid]
        neurons[bid] = Neuron(
            id=bid,
            nodes=sk[["x", "y", "z"]].to_numpy(float) * voxel_nm,
            parents=parents,
            radius=sk.radius.to_numpy(float) * voxel_nm,
            pre=s.loc[s.type == "pre", ["x", "y", "z"]].to_numpy(float) * voxel_nm,
            post=s.loc[s.type == "post", ["x", "y", "z"]].to_numpy(float) * voxel_nm,
            meta={"type": types.get(bid), "instance": instances.get(bid),
                  **{c: m.get(bid) for c, m in extra.items()}},
        )

    if with_meshes:
        from navis.interfaces import neuprint as nvn
        meshes = nvn.fetch_mesh_neuron(NC(bodyId=body_ids), client=client, lod=None)
        for m in meshes:
            # navis は neuPrint のメッシュをボクセル単位に変換して返す → nm に戻す
            n = neurons[int(m.id)]
            n.mesh_vertices = np.asarray(m.vertices, float) * voxel_nm
            n.mesh_faces = np.asarray(m.faces, int)

    _, conn = _retry(fetch_adjacencies, sources=NC(bodyId=body_ids), targets=NC(bodyId=body_ids), client=client)
    edges = (conn.groupby(["bodyId_pre", "bodyId_post"]).weight.sum().reset_index()
             .rename(columns={"bodyId_pre": "pre", "bodyId_post": "post"}))

    # 個々のシナプス（アプローチ F 用）: T-bar と PSD の中点を位置とする
    sc = _retry(fetch_synapse_connections, NC(bodyId=body_ids), NC(bodyId=body_ids),
                SC(confidence=min_confidence), client=client)
    synapses = pd.DataFrame({
        "pre": sc.bodyId_pre, "post": sc.bodyId_post,
        "x": (sc.x_pre + sc.x_post) / 2 * voxel_nm,
        "y": (sc.y_pre + sc.y_post) / 2 * voxel_nm,
        "z": (sc.z_pre + sc.z_post) / 2 * voxel_nm,
    })
    return Dataset(neurons=neurons, edges=edges, synapses=synapses, name=f"neuprint_{dataset}".replace(":", "_"))


# ----------------------------------------------------------------------
# FlyWire (FAFB, 公開版 materialization 783)
# ----------------------------------------------------------------------
def _xyz(df: pd.DataFrame, prefix: str) -> np.ndarray:
    for cols in ([f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"], [f"{prefix}_pt_position"]):
        if all(c in df.columns for c in cols):
            if len(cols) == 1:
                return np.vstack(df[cols[0]].values).astype(float)
            return df[cols].to_numpy(float)
    raise KeyError(f"cannot find {prefix} coordinates in columns {list(df.columns)}")


def load_flywire(root_ids, materialization: int = 783, with_meshes: bool = False,
                 min_cleft_score: int = 50, types: dict | None = None) -> Dataset:
    """※ fafbseg の返り値の列名はバージョンで変わることがあるので、初回は中身を確認すること。

    types: {root_id: 細胞タイプ名}（アプローチ E 用）。FlyWire Codex からダウンロードできる
           細胞タイプの表（consolidated_cell_types.csv など）から作る。
    """
    from fafbseg import flywire

    root_ids = [int(r) for r in root_ids]
    sks = flywire.get_skeletons(root_ids, dataset=materialization)
    syn = flywire.get_synapses(root_ids, pre=True, post=True, materialization=materialization,
                               filtered=True, min_score=min_cleft_score)

    neurons = {}
    for sk in sks:
        rid = int(sk.id)
        nd = sk.nodes
        idx = {nid: k for k, nid in enumerate(nd.node_id)}
        parents = np.array([idx.get(p, -1) for p in nd.parent_id], dtype=int)
        radius = nd.radius.to_numpy(float) if "radius" in nd else None
        out_syn = syn[syn.pre == rid]
        in_syn = syn[syn.post == rid]
        neurons[rid] = Neuron(id=rid, nodes=nd[["x", "y", "z"]].to_numpy(float), parents=parents,
                              radius=radius, pre=_xyz(out_syn, "pre"), post=_xyz(in_syn, "post"),
                              meta={"type": (types or {}).get(rid)})

    if with_meshes:
        for rid in root_ids:
            m = flywire.get_mesh_neuron(rid)
            neurons[rid].mesh_vertices = np.asarray(m.vertices, float)
            neurons[rid].mesh_faces = np.asarray(m.faces, int)

    inside = syn[syn.pre.isin(root_ids) & syn.post.isin(root_ids)]
    edges = inside.groupby(["pre", "post"]).size().rename("weight").reset_index()
    mid = (_xyz(inside, "pre") + _xyz(inside, "post")) / 2
    synapses = pd.DataFrame({"pre": inside.pre.values, "post": inside.post.values,
                             "x": mid[:, 0], "y": mid[:, 1], "z": mid[:, 2]})
    return Dataset(neurons=neurons, edges=edges, synapses=synapses, name=f"flywire_{materialization}")
