#!/usr/bin/env python3
"""Reusable canonical GLMY5 extractor.

This module isolates the final five-channel GLMY/Hu path-topology definition
from any DNA-, protein-, or task-specific graph construction.

Frozen GLMY5 order:
    (b0, b1, b2, c1, c2)

For one fixed graph G=(V,E):
    depth 2: b0 = infinite H0 count; c1 = top-dimensional H1 count
    depth 3: b1 = infinite H1 count; c2 = top-dimensional H2 count
    depth 4: b2 = infinite H2 count

The c channels are top-depth path-cycle abundances at a truncation depth. The b
channels are one-level-stabilized GLMY homology counts. These are quantities of
the study-specific GLMY/Hu face-completed path-complex implemented by
Digraph_Hu.py. They must not be described as standard GF(2) regular path-
homology Betti numbers.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

GLMY5_ORDER = ("b0", "b1", "b2", "c1", "c2")
GLMY8_PARENT_ORDER = ("b0", "c1", "b1", "c2", "b2", "c3", "b3", "c4")


@dataclass(frozen=True)
class GLMY5Result:
    b0: int
    b1: int
    b2: int
    c1: int
    c2: int
    qc: Mapping[str, bool]
    diagnostics: Mapping[str, int]

    @property
    def vector(self) -> Tuple[int, int, int, int, int]:
        return (self.b0, self.b1, self.b2, self.c1, self.c2)

    def feature_dict(self) -> Dict[str, int]:
        return dict(zip(GLMY5_ORDER, self.vector))

    def to_dict(self) -> Dict[str, Any]:
        out = asdict(self)
        out["order"] = list(GLMY5_ORDER)
        out["vector"] = list(self.vector)
        out["method"] = "GLMY/Hu face-completed path-topology homology"
        out["guardrail"] = (
            "Not standard GF(2) regular path-homology Betti numbers."
        )
        return out


def load_digraph_class(path: Path):
    path = Path(path).expanduser().resolve()
    spec = importlib.util.spec_from_file_location("Digraph_Hu_glmy5_final", str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load Digraph_Hu module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Digraph


def _diagram_counts(DigraphClass, vertices, edges, depth: int):
    dg = DigraphClass(list(vertices), list(edges), int(depth))
    diagram = dg.get_persistence()
    infinite: Dict[int, int] = {}
    total: Dict[int, int] = {}
    finite: Dict[int, int] = {}
    for dim in range(int(depth)):
        bars = diagram.get(str(dim), [])
        total[dim] = len(bars)
        infinite[dim] = sum(1 for _, death in bars if death == -1)
        finite[dim] = total[dim] - infinite[dim]
    return infinite, total, finite


def compute_glmy5(
    DigraphClass,
    vertices: Sequence[Any],
    weighted_edges: Iterable[Tuple[Any, Any, float]],
    *,
    fail_on_qc: bool = True,
) -> GLMY5Result:
    """Compute final GLMY5 channels on one already-constructed graph.

    Parameters
    ----------
    DigraphClass:
        Digraph class loaded from the audited Digraph_Hu.py.
    vertices:
        All graph vertices, including isolated vertices.
    weighted_edges:
        Iterable of (source, target, weight). The graph construction itself is
        domain-specific and must be frozen separately.
    fail_on_qc:
        Raise if the audited depth invariants fail.
    """
    vertices = list(vertices)
    edges = list(weighted_edges)

    if len(vertices) == 0:
        qc = {
            "c1_depth2_equals_total_H1": True,
            "c2_depth3_equals_total_H2": True,
            "b0_stable_depth2_3_4": True,
            "b1_stable_depth3_4": True,
            "b1_le_c1": True,
            "b2_le_c2": True,
        }
        return GLMY5Result(0, 0, 0, 0, 0, qc, {})

    inf2, total2, finite2 = _diagram_counts(DigraphClass, vertices, edges, 2)
    inf3, total3, finite3 = _diagram_counts(DigraphClass, vertices, edges, 3)
    inf4, total4, finite4 = _diagram_counts(DigraphClass, vertices, edges, 4)

    b0 = int(inf2.get(0, 0))
    c1 = int(inf2.get(1, 0))
    b1 = int(inf3.get(1, 0))
    c2 = int(inf3.get(2, 0))
    b2 = int(inf4.get(2, 0))

    qc = {
        "c1_depth2_equals_total_H1": c1 == int(total2.get(1, 0)),
        "c2_depth3_equals_total_H2": c2 == int(total3.get(2, 0)),
        "b0_stable_depth2_3_4": (
            int(inf2.get(0, 0)) == int(inf3.get(0, 0)) == int(inf4.get(0, 0))
        ),
        "b1_stable_depth3_4": int(inf3.get(1, 0)) == int(inf4.get(1, 0)),
        "b1_le_c1": b1 <= c1,
        "b2_le_c2": b2 <= c2,
    }

    diagnostics = {
        "depth2_H0_inf": int(inf2.get(0, 0)),
        "depth2_H1_inf": int(inf2.get(1, 0)),
        "depth2_H1_total": int(total2.get(1, 0)),
        "depth3_H0_inf": int(inf3.get(0, 0)),
        "depth3_H1_inf": int(inf3.get(1, 0)),
        "depth3_H1_total": int(total3.get(1, 0)),
        "depth3_H1_finite": int(finite3.get(1, 0)),
        "depth3_H2_inf": int(inf3.get(2, 0)),
        "depth3_H2_total": int(total3.get(2, 0)),
        "depth4_H0_inf": int(inf4.get(0, 0)),
        "depth4_H1_inf": int(inf4.get(1, 0)),
        "depth4_H2_inf": int(inf4.get(2, 0)),
        "depth4_H2_total": int(total4.get(2, 0)),
        "depth4_H2_finite": int(finite4.get(2, 0)),
    }

    if fail_on_qc and not all(qc.values()):
        bad = [name for name, ok in qc.items() if not ok]
        raise RuntimeError(f"GLMY5 invariant QC failed: {bad}; diagnostics={diagnostics}")

    return GLMY5Result(b0, b1, b2, c1, c2, qc, diagnostics)


def compute_glmy5_from_count(
    DigraphClass,
    n_vertices: int,
    weighted_edges: Iterable[Tuple[int, int, float]],
    *,
    fail_on_qc: bool = True,
) -> GLMY5Result:
    """Convenience wrapper for graphs whose vertices are 0..n_vertices-1."""
    return compute_glmy5(
        DigraphClass,
        list(range(int(n_vertices))),
        weighted_edges,
        fail_on_qc=fail_on_qc,
    )


def _self_test(digraph_file: Path) -> None:
    DigraphClass = load_digraph_class(digraph_file)

    isolated = compute_glmy5(DigraphClass, [0, 1, 2], [])
    assert isolated.vector == (3, 0, 0, 0, 0), isolated.vector

    # Audited TCGA graph example used in the canonical GLMY8 review.
    sequence = "AATGATGCACATCATGGTGGC"
    group = ("A", "GC")
    threshold = 4

    positions_by_token = {}
    for token in group:
        positions_by_token[token] = [
            i + 1 for i in range(len(sequence) - len(token) + 1)
            if sequence[i:i + len(token)] == token
        ]

    length_by_position: Dict[int, int] = {}
    for token in group:
        for pos in positions_by_token[token]:
            length_by_position[pos] = max(len(token), length_by_position.get(pos, 0))
    positions = sorted(length_by_position)

    edges: List[Tuple[int, int, float]] = []
    for ui, u in enumerate(positions):
        for vi, v in enumerate(positions):
            if ui == vi:
                continue
            lu = length_by_position[u]
            lv = length_by_position[v]
            distance = abs(v - u)
            overlap_len = min(lu, lv)
            if distance > overlap_len - 1:
                distance = distance - overlap_len + 1
            if distance <= threshold:
                edges.append((ui, vi, float(distance)))

    demo = compute_glmy5_from_count(DigraphClass, len(positions), edges)
    expected = (2, 0, 1, 12, 19)
    assert demo.vector == expected, (demo.vector, expected)
    assert all(demo.qc.values()), demo.qc
    print("[SELF-TEST] PASS")
    print("[SELF-TEST] GLMY5 order:", GLMY5_ORDER)
    print("[SELF-TEST] audited demo:", demo.vector)


def main() -> int:
    parser = argparse.ArgumentParser(description="Reusable final GLMY5 extractor self-test")
    parser.add_argument("--digraph-file", required=True)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if not args.self_test:
        raise SystemExit("This module is an API. Use --self-test or scripts/glmy5_edge_cli.py.")
    _self_test(Path(args.digraph_file))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
