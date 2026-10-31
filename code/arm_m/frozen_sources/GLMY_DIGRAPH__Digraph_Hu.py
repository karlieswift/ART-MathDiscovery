#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import argparse
from typing import Dict, List, Tuple, Iterable, Optional, Set, Any

# 尝试导入可选依赖
try:
    import numpy as np
except Exception:
    np = None

try:
    from bitarray import bitarray  # pip install bitarray
except Exception:
    bitarray = None


# ========= 基础工具 =========

def xor_sorted(a: List[int], b: List[int]) -> List[int]:
    """对两个升序整型列表做对称差(异或)，返回升序列表（纯 Python 线性复杂度）。"""
    i, j, la, lb = 0, 0, len(a), len(b)
    out = []
    while i < la and j < lb:
        ai, bj = a[i], b[j]
        if ai == bj:
            i += 1; j += 1
        elif ai < bj:
            out.append(ai); i += 1
        else:
            out.append(bj); j += 1
    if i < la: out.extend(a[i:])
    if j < lb: out.extend(b[j:])
    return out


# ========= 三种后端的统一接口 =========

class BoundaryBackend:
    """边界矩阵后端抽象类。每一列代表一个链的边界（Z2），支持 pivot 与 XOR。"""
    def __init__(self, n_rows: int):
        self.n_rows = n_rows

    def make_col(self, indices: List[int]) -> Any:
        raise NotImplementedError

    def pivot(self, col: Any) -> int:
        """返回该列的最低 1 的行索引（按升序，则为最后一个索引）；若为空列，返 -1。"""
        raise NotImplementedError

    def xor_inplace(self, col_a: Any, col_b: Any) -> Any:
        """令 col_a ^= col_b，并返回 col_a。"""
        raise NotImplementedError


class BackendList(BoundaryBackend):
    """每列保存为升序整型索引列表（默认后端）。"""
    def make_col(self, indices: List[int]) -> List[int]:
        return sorted(indices)

    def pivot(self, col: List[int]) -> int:
        return col[-1] if col else -1

    def xor_inplace(self, col_a: List[int], col_b: List[int]) -> List[int]:
        return xor_sorted(col_a, col_b)


class BackendNumpy(BoundaryBackend):
    """每列保存为 numpy.bool_ 向量。内存开销约 n_rows bits/列。
    注意：需要 numpy 可用。"""
    def __init__(self, n_rows: int):
        if np is None:
            raise RuntimeError("numpy 不可用，请改用 --backend list 或安装 numpy。")
        super().__init__(n_rows)

    def make_col(self, indices: List[int]):
        col = np.zeros(self.n_rows, dtype=bool)
        if indices:
            col[np.array(indices, dtype=int)] = True
        return col

    def pivot(self, col) -> int:
        nz = np.flatnonzero(col)
        return int(nz[-1]) if nz.size else -1

    def xor_inplace(self, col_a, col_b):
        col_a ^= col_b
        return col_a


class BackendBitarray(BoundaryBackend):
    """每列保存为 bitarray。需要 pip install bitarray。
    注意：pivot 查找使用 rfind 近似实现。"""
    def __init__(self, n_rows: int):
        if bitarray is None:
            raise RuntimeError("bitarray 不可用，请改用 --backend list 或安装 bitarray。")
        super().__init__(n_rows)

    def make_col(self, indices: List[int]):
        col = bitarray(self.n_rows)
        col.setall(False)
        for i in indices:
            col[i] = True
        return col

    def pivot(self, col) -> int:
        # bitarray 没有直接 rfind(1)，用反向扫描
        # 这种写法是 O(n_rows)，但 XOR 很快；对大规模建议用 list/numpy 后端。
        try:
            # 高版本 bitarray 支持 .rindex(True)
            return col.rindex(True)
        except Exception:
            # 兼容保底
            for idx in range(self.n_rows - 1, -1, -1):
                if col[idx]:
                    return idx
            return -1

    def xor_inplace(self, col_a, col_b):
        col_a ^= col_b
        return col_a


# ========= Digraph 主体 =========

class Digraph:
    False_value = 99999999

    def __init__(
        self,
        V: Iterable[Any],
        WG: Iterable[Tuple[Any, Any, float]],
        dim: int,
        backend: str = "list",
        simple_paths: bool = False,
        deduplicate: bool = True,
        max_repeat: Optional[int] = None,
    ):
        """
        V: 顶点可迭代（元素可哈希）
        WG: 带权有向边 [(u, v, w), ...]，w 可转 float
        dim: 生成到的最大路径长度（dim=4 => 生成长度1..4的路径）
        backend: 'list'（默认）、'numpy'、'bitarray'
        simple_paths: True 时拒绝路径内重复顶点（强力抑制爆炸）
        deduplicate: True 时每个长度只保留唯一路径
        max_repeat: 顶点在同一路径中允许出现的最大次数（None 不限制；=1 等价 simple_paths）
        """
        self.V = list(V)
        self.dim = int(dim)
        self.simple_paths = bool(simple_paths)
        self.deduplicate = bool(deduplicate)
        self.max_repeat = max_repeat if max_repeat is None else int(max_repeat)

        # 出边表：tail -> List[(to, weight)]
        v_out: Dict[Any, List[Tuple[Any, float]]] = {v: [] for v in self.V}
        for u, v, w in WG:
            v_out[u].append((v, float(w)))
        self._v_out = v_out

        # 选择后端
        self.backend_name = backend.lower()
        self.backend_cls = {"list": BackendList, "numpy": BackendNumpy, "bitarray": BackendBitarray}[self.backend_name]

        # 生成路径与权重（保持你的原始规则）
        self.P, self.PW = self._build_paths()

        # 生成“关联复形”缺失的面（False_value）
        self._add_missing_faces()

        # 排序并构建边界矩阵列（只为真实列构造）
        self._sort_and_build_boundary()

        # 输出容器
        self.diagram: Dict[str, List[List[float]]] = {}

    # ---------- 路径生成 ----------

    def _check_repeat_ok(self, path: Tuple[Any, ...], new_v: Any) -> bool:
        """检查追加 new_v 后是否满足 simple_paths / max_repeat 约束。"""
        if self.simple_paths:
            return new_v not in path
        if self.max_repeat is not None:
            # 最多出现 max_repeat 次
            # 追加后次数:
            cnt = 1
            for x in path:
                if x == new_v:
                    cnt += 1
                    if cnt > self.max_repeat:
                        return False
            return True
        return True

    def _build_paths(self) -> Tuple[List[Tuple[Any, ...]], List[float]]:
        """按你的逻辑生成长度 1..dim 的路径与权重（路径权重=max 边权）。"""
        P_by_len: List[List[Tuple[Any, ...]]] = [[] for _ in range(self.dim + 1)]
        W_by_len: List[List[float]] = [[] for _ in range(self.dim + 1)]

        # L=1: 顶点
        for v in self.V:
            P_by_len[1].append((v,))
            W_by_len[1].append(0.0)

        # L=2: 边
        for u, outs in self._v_out.items():
            for v, w in outs:
                P_by_len[2].append((u, v))
                W_by_len[2].append(float(w))

        # L>=3: 扩展
        for L in range(3, self.dim + 1):
            prev_paths = P_by_len[L - 1]
            prev_w = W_by_len[L - 1]
            if not prev_paths:
                break

            curr_set: Set[Tuple[Any, ...]] = set() if self.deduplicate else None
            curr_paths: List[Tuple[Any, ...]] = []
            curr_w: List[float] = []

            for path, pw in zip(prev_paths, prev_w):
                tail = path[-1]
                outs = self._v_out.get(tail, ())
                if not outs:
                    continue
                for to, w in outs:
                    if not self._check_repeat_ok(path, to):
                        continue
                    new_path = path + (to,)
                    new_w = pw if pw >= w else w
                    if curr_set is not None:
                        if new_path in curr_set:
                            continue
                        curr_set.add(new_path)
                    curr_paths.append(new_path)
                    curr_w.append(new_w)

            P_by_len[L] = curr_paths
            W_by_len[L] = curr_w

        # 汇总
        P: List[Tuple[Any, ...]] = []
        W: List[float] = []
        for L in range(1, self.dim + 1):
            P.extend(P_by_len[L])
            W.extend(W_by_len[L])
        return P, W

    # ---------- 面补齐 ----------

    def _add_missing_faces(self) -> None:
        """为 L>=3 的路径补全所有去掉一个顶点的面；新面权重 = False_value。"""
        path2idx = {p: i for i, p in enumerate(self.P)}
        appended: List[Tuple[Any, ...]] = []
        appended_w: List[float] = []
        for p in self.P:
            L = len(p)
            if L >= 3:
                for k in range(L):
                    face = p[:k] + p[k + 1:]
                    if face not in path2idx:
                        path2idx[face] = len(self.P) + len(appended)
                        appended.append(face)
                        appended_w.append(self.False_value)
        if appended:
            self.P.extend(appended)
            self.PW.extend(appended_w)

    # ---------- 排序与边界矩阵 ----------

    def _sort_and_build_boundary(self) -> None:
        # 排序键： (weight, len(path))
        order = sorted(range(len(self.P)), key=lambda i: (self.PW[i], len(self.P[i])))
        self.P = [self.P[i] for i in order]
        self.PW = [self.PW[i] for i in order]

        # 真实列数（非 False_value）
        N = 0
        while N < len(self.PW) and self.PW[N] != self.False_value:
            N += 1
        self._N_real = N

        # 行数 = 全部链元个数
        n_rows = len(self.P)

        # 反向索引
        path2pos = {p: i for i, p in enumerate(self.P)}

        # 选择后端
        self.backend: BoundaryBackend = self.backend_cls(n_rows)

        # 构造前 N 列的边界
        boundary_cols: List[Any] = []
        for j in range(N):
            p = self.P[j]
            if len(p) <= 1:
                boundary_cols.append(self.backend.make_col([]))
            else:
                faces = [path2pos[p[:k] + p[k + 1:]] for k in range(len(p))]
                faces.sort()
                boundary_cols.append(self.backend.make_col(faces))

        self.boundary_matrix = boundary_cols

    # ---------- 持久同调：Z2 列消元 ----------

    def get_persistence(self) -> Dict[str, List[List[float]]]:
        """
        与原逻辑等价：Z2 上的列消元（低星序），生成 self.diagram:
        {"0": [[birth, death], ...], "1": [...], ...}；death=-1 表示无穷条形
        """
        N = self._N_real
        P = self.P
        PW = self.PW
        boundary = self.boundary_matrix

        # L[row] = col 的“低索引 -> 列”匹配表；-1 表示当前未被占用
        Lmatch = [-1] * len(P)

        # 空边界列 -> 候选无穷条形出生
        infinite_birth_idx: Set[int] = set()

        for j in range(N):
            col = boundary[j]
            while True:
                low = self.backend.pivot(col)
                if low == -1:
                    break
                locked = Lmatch[low]
                if locked != -1:
                    col = self.backend.xor_inplace(col, boundary[locked])
                else:
                    break
            boundary[j] = col  # 更新（对 numpy/bitarray 属于引用）

            low = self.backend.pivot(col)
            if low == -1:
                infinite_birth_idx.add(j)
            else:
                Lmatch[low] = j

        diagram: Dict[str, List[List[float]]] = {}

        # 有限条形配对
        for j in range(N):
            col = boundary[j]
            low = self.backend.pivot(col)
            if low != -1:
                if low in infinite_birth_idx:
                    infinite_birth_idx.remove(low)
                b = PW[low]
                d = PW[j]
                dim = len(P[low]) - 1
                if d > b:  # 与原代码一致
                    diagram.setdefault(str(dim), []).append([b, d])

        # 无穷条形
        for idx in sorted(infinite_birth_idx):
            dim = len(P[idx]) - 1
            diagram.setdefault(str(dim), []).append([PW[idx], -1])

        self.diagram = diagram
        return diagram


# ========= 演示 / 命令行 =========

def demo_graph() -> Tuple[List[int], List[Tuple[int, int, float]]]:
    """
    一个小的示例图（含环），便于 sanity check。
    顶点: 1..5
    边: 1->2(0.2), 2->3(0.5), 3->1(0.7), 2->4(0.1), 4->5(0.4)
    """
    V = [1, 2, 3, 4, 5]
    WG = [
        (1, 2, 0.2),
        (2, 3, 0.5),
        (3, 1, 0.7),  # 构成 1-2-3-1 环
        (2, 4, 0.1),
        (4, 5, 0.4),
    ]
    return V, WG


def main():
    ap = argparse.ArgumentParser(description="Path-complex persistent homology on weighted digraphs (optimized).")
    ap.add_argument("--dim", type=int, default=4, help="Max path length to generate (default: 4)")
    ap.add_argument("--backend", type=str, default="list", choices=["list", "numpy", "bitarray"],
                    help="Boundary matrix backend (default: list)")
    ap.add_argument("--simple-paths", action="store_true",
                    help="Disallow repeated vertices within a path (strong pruning)")
    ap.add_argument("--no-dedupe", action="store_true",
                    help="Disable per-length path deduplication (not recommended)")
    ap.add_argument("--max-repeat", type=int, default=None,
                    help="Max times a vertex can appear within a path (1 == simple paths)")
    ap.add_argument("--use-demo", action="store_true", help="Run with an internal demo graph")
    ap.add_argument("--edge-file", type=str, default=None,
                    help="Optional CSV file with columns: u,v,weight (headers allowed)")
    args = ap.parse_args()

    if args.use_demo and args.edge_file:
        raise SystemExit("只能二选一：--use-demo 或 --edge-file。")

    if args.use_demo:
        V, WG = demo_graph()
    elif args.edge_file:
        import csv
        WG = []
        Vset = set()
        with open(args.edge_file, "r", newline="") as f:
            reader = csv.reader(f)
            # 既兼容有表头，也兼容无表头；跳过非数值行
            for row in reader:
                if not row or len(row) < 3:
                    continue
                try:
                    u = row[0].strip()
                    v = row[1].strip()
                    w = float(row[2])
                except Exception:
                    # 尝试跳过头行
                    if row[0].lower().startswith("u") or row[0].lower().startswith("src"):
                        continue
                    else:
                        # 忽略解析失败行
                        continue
                WG.append((u, v, w))
                Vset.add(u); Vset.add(v)
        V = sorted(Vset)
    else:
        V, WG = demo_graph()

    g = Digraph(
        V=V,
        WG=WG,
        dim=args.dim,
        backend=args.backend,
        simple_paths=args.simple_paths or (args.max_repeat == 1),
        deduplicate=not args.no_dedupe,
        max_repeat=args.max_repeat,
    )

    diagram = g.get_persistence()

    # 打印摘要
    print("=== Summary ===")
    print(f"#cells (total) = {len(g.P)}; #real columns = {g._N_real}; backend = {args.backend}")
    dims = sorted(int(k) for k in diagram.keys())
    if not dims:
        print("No bars found.")
    else:
        for d in dims:
            bars = diagram[str(d)]
            finite = [b for b in bars if b[1] != -1]
            infinite = [b for b in bars if b[1] == -1]
            print(f"H_{d}: finite={len(finite)}, infinite={len(infinite)}")
    print("\nSample bars (up to 5 per dim):")
    for k, bars in sorted(diagram.items(), key=lambda kv: int(kv[0])):
        show = bars[:5]
        print(f"  dim {k}: {show}")

    # 若需要，暴露 P/PW/部分边界
    # print("\nFirst 10 cells (path, weight):")
    # for i in range(min(10, len(g.P))):
    #     print(i, g.P[i], g.PW[i])


if __name__ == "__main__":
    main()

