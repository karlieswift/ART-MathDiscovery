#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
79A2_build_full_glmy5_fresh_features.py

Submission-oriented fresh GLMY feature builder for TCGA_HRD_NV.

This script recomputes the five frozen predictive GLMY channels FROM SCRATCH
for every token-group graph. It does not reuse Revision74 b0/b1/b2 matrices and
it does not reuse Revision79A c1/c2 matrices.

Frozen five-channel definition
------------------------------
For the same graph G=(V,E) at a selected adjusted-distance threshold:

    depth 2: b0 = # infinite H0 bars
             c1 = # top-dimensional H1 bars (= len(H1), fail-closed checked)

    depth 3: b1 = # infinite H1 bars
             c2 = # top-dimensional H2 bars (= len(H2), fail-closed checked)

    depth 4: b2 = # infinite H2 bars

Canonical channel order used everywhere in this script:

    [b0, b1, b2, c1, c2]

Interpretation
--------------
- b0/b1/b2: stable GLMY homology channels.
- c1/c2: transient top-dimensional path-cycle abundance.

The script also fail-closes on the structural invariants established in the
Rev78/78B/78C audits:

    b0(depth2) == b0(depth3) == b0(depth4)
    b1(depth3) == b1(depth4)
    c1(depth2) == len(H1_depth2)
    c2(depth3) == len(H2_depth3)
    b1 <= c1
    b2 <= c2

Frozen TCGA representation
--------------------------
- exact-event-deduplicated hg38-centered 21-mer SNV contexts
- frozen Revision60v2 2,003-patient model cohort
- reverse-complement canonicalization to the six pyrimidine mutation classes
- 20 DNA 1-/2-mer tokens
- singles + unordered pairs = 210 token groups
- adjusted-distance symmetric/bidirected graph
- distances d=2,3,4,5 by default
- native real contexts plus deterministic outer-dinucleotide-preserving shuffle
- patient aggregation:

      6 RC6 classes
    x 4 blocks (WT, MUT, delta, absdelta)
    x 2 token-group sizes
    x 5 GLMY channels
    x 3 token summaries (mean,std,max)
    x 3 mutation summaries (mean,std,max)
    = 2,160 patient features per distance
    =   432 features per GLMY channel per distance

No HRD value is used in feature calculation. HRD/LOH/LST/TAI columns from the
frozen model matrix are copied to output tables only as downstream labels.

Code-submission intent
----------------------
This file intentionally contains the canonicalization, dinucleotide-control,
token-graph construction, aggregation, and five-channel GLMY extraction logic
in one place. The only project-specific algorithm file it imports is
Digraph_Hu.py.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import time
import zlib
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


CANCERS = ("BRCA", "OV", "LUAD", "PRAD")
RC6_CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CHANNELS = ("b0", "b1", "b2", "c1", "c2")
TOKEN_SUMMARIES = ("mean", "std", "max")
MUTATION_SUMMARIES = ("mean", "std", "max")
BLOCKS = ("wt", "mut", "delta", "absdelta")
CONTROLS = ("real", "outer_dinuc_shuffle")
LABEL_COLS = ("HRD", "HRD_LOH", "HRD_LST", "HRD_TAI")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
RC_MAP = {
    "A>C": "T>G", "A>G": "T>C", "A>T": "T>A",
    "G>A": "C>T", "G>C": "C>G", "G>T": "C>A",
    "C>A": "C>A", "C>G": "C>G", "C>T": "C>T",
    "T>A": "T>A", "T>C": "T>C", "T>G": "T>G",
}

# Worker globals.
_W_DIGRAPH = None
_W_TOKEN_GROUPS = None
_W_TOKEN_LENGTHS = None
_W_GROUP_SIZES = None
_W_DISTANCES = None
_W_MAX_DISTANCE = None
_W_CONTROL_SEED = None
_W_CONTROLS = None
_cached_glmy5 = None


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_csv(value):
    return [x.strip() for x in str(value).split(",") if x.strip()]


def resolve_path(root: Path, value) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def normalize_patient(value):
    return str(value)[:12]


def clean_seq(value):
    if pd.isna(value):
        return ""
    return str(value).upper().replace("U", "T")


def revcomp(seq):
    return str(seq).translate(COMP)[::-1].upper()


def resolve_duplicate_key(frame):
    chrom_candidates = ("chrom", "chr", "chromosome", "Chromosome")
    pos_candidates = ("pos", "pos1", "position", "start_position", "Start_Position")
    chrom = next((x for x in chrom_candidates if x in frame.columns), None)
    pos = next((x for x in pos_candidates if x in frame.columns), None)
    if (
        chrom is None
        or pos is None
        or "patient_barcode" not in frame.columns
        or "ref" not in frame.columns
        or "alt" not in frame.columns
    ):
        return []
    return ["patient_barcode", chrom, pos, "ref", "alt"]


def canonicalize_context(wt, mut, ref, alt):
    wt = clean_seq(wt)
    mut = clean_seq(mut)
    ref = str(ref).upper()
    alt = str(alt).upper()

    if (
        len(wt) != 21
        or len(mut) != 21
        or len(ref) != 1
        or len(alt) != 1
        or ref not in "ACGT"
        or alt not in "ACGT"
        or ref == alt
    ):
        raise ValueError("invalid 21-mer SNV context")

    if wt[10] != ref or mut[10] != alt:
        raise ValueError("center/ref/alt mismatch")

    original_class = f"{ref}>{alt}"
    if original_class not in RC_MAP:
        raise ValueError(f"unsupported substitution {original_class}")

    if ref in {"A", "G"}:
        wt = revcomp(wt)
        mut = revcomp(mut)

    rc6 = RC_MAP[original_class]

    if wt[:10] != mut[:10] or wt[11:] != mut[11:]:
        raise ValueError("WT/MUT differ outside the center")

    return wt, mut, rc6


def read_vocab(path):
    frame = pd.read_csv(path, sep="\t")
    if "Token" not in frame.columns:
        raise SystemExit(f"{path}: missing Token")
    tokens = sorted(
        {
            str(value).strip().upper()
            for value in frame["Token"].dropna()
            if str(value).strip()
        },
        key=lambda x: (len(x), x),
    )
    if len(tokens) != 20:
        raise SystemExit(f"Expected 20 DNA 1-to-2mer tokens, got {len(tokens)}")
    return tokens


def make_token_groups(tokens):
    groups = []
    for size in (1, 2):
        groups.extend(tuple(group) for group in combinations(tokens, size))
    if len(groups) != 210:
        raise RuntimeError(f"Expected 210 token groups, got {len(groups)}")
    return groups


def stable_rng(seed, event_id, control):
    value = zlib.crc32(f"{event_id}|{control}".encode("utf-8"))
    return np.random.default_rng((int(seed) + int(value)) % (2**32))


def euler_dinuc_shuffle(seq, rng, tries=100):
    """Alternative Euler trail preserving the exact directed dinucleotide multiset."""
    seq = str(seq)
    if len(seq) <= 2:
        return seq, False

    edge_counts = Counter(zip(seq[:-1], seq[1:]))
    start = seq[0]
    end = seq[-1]

    for _ in range(tries):
        adjacency = defaultdict(list)
        for (source, target), count in edge_counts.items():
            adjacency[source].extend([target] * count)
        for source in list(adjacency):
            rng.shuffle(adjacency[source])

        stack = [start]
        reverse_path = []
        while stack:
            vertex = stack[-1]
            if adjacency[vertex]:
                stack.append(adjacency[vertex].pop())
            else:
                reverse_path.append(stack.pop())

        trail = reverse_path[::-1]
        if len(trail) != len(seq):
            continue
        candidate = "".join(trail)
        if candidate[0] != start or candidate[-1] != end:
            continue
        if Counter(zip(candidate[:-1], candidate[1:])) != edge_counts:
            continue
        if candidate != seq:
            return candidate, True

    return seq, False


def make_control_pair(wt, mut, control, rng):
    if mut[9] != wt[9] or mut[11] != wt[11]:
        raise ValueError("WT/MUT center flanks differ")

    left = wt[:9]
    right = wt[12:]

    if control != "outer_dinuc_shuffle":
        raise KeyError(control)

    left_new, left_changed = euler_dinuc_shuffle(left, rng)
    right_new, right_changed = euler_dinuc_shuffle(right, rng)

    wt_control = left_new + wt[9:12] + right_new
    mut_control = left_new + mut[9:12] + right_new
    return wt_control, mut_control, left_changed, right_changed


def mono_counts(seq):
    return Counter(seq)


def dinuc_counts(seq):
    return Counter(seq[i:i + 2] for i in range(len(seq) - 1))


def sbs96_key(wt, mut):
    return (wt[9], wt[10], mut[10], wt[11])


def get_position_maps(seq, token_lengths):
    maps = {}
    for length in token_lengths:
        local = defaultdict(list)
        for index in range(len(seq) - length + 1):
            local[seq[index:index + length]].append(index + 1)
        maps[length] = local
    return maps


def group_positions(position_maps, group):
    length_by_position = {}
    for token in group:
        length = len(token)
        for position in position_maps[length].get(token, ()):
            length_by_position[position] = max(
                length,
                length_by_position.get(position, 0),
            )
    return sorted(length_by_position), length_by_position


def adjusted_distance(u, v, lu, lv):
    distance = abs(v - u)
    overlap_len = min(lu, lv)
    if distance > overlap_len - 1:
        distance = distance - overlap_len + 1
    return int(distance)


def load_digraph_class(path):
    spec = importlib.util.spec_from_file_location("Digraph_Hu_glmy5_fresh", str(path))
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.Digraph


def inf_count(diagram, dimension):
    bars = diagram.get(str(dimension), [])
    return sum(1 for _, death in bars if death == -1), len(bars)


def calculate_glmy5(Digraph, n_vertices, weighted_edges):
    """Return canonical (b0,b1,b2,c1,c2) for one graph."""
    n_vertices = int(n_vertices)
    if n_vertices == 0:
        return (0, 0, 0, 0, 0)

    vertices = list(range(n_vertices))
    edges = list(weighted_edges)

    dg2 = Digraph(vertices, edges, 2)
    dg2.get_persistence()
    b0_2, _ = inf_count(dg2.diagram, 0)
    c1, h1_2_len = inf_count(dg2.diagram, 1)

    dg3 = Digraph(vertices, edges, 3)
    dg3.get_persistence()
    b0_3, _ = inf_count(dg3.diagram, 0)
    b1, _ = inf_count(dg3.diagram, 1)
    c2, h2_3_len = inf_count(dg3.diagram, 2)

    dg4 = Digraph(vertices, edges, 4)
    dg4.get_persistence()
    b0_4, _ = inf_count(dg4.diagram, 0)
    b1_4, _ = inf_count(dg4.diagram, 1)
    b2, _ = inf_count(dg4.diagram, 2)

    if c1 != h1_2_len:
        raise RuntimeError(
            f"GLMY invariant failed: depth2 c1={c1} != len(H1)={h1_2_len}"
        )
    if c2 != h2_3_len:
        raise RuntimeError(
            f"GLMY invariant failed: depth3 c2={c2} != len(H2)={h2_3_len}"
        )
    if not (b0_2 == b0_3 == b0_4):
        raise RuntimeError(
            f"GLMY invariant failed: b0 depth2/3/4={b0_2}/{b0_3}/{b0_4}"
        )
    if b1 != b1_4:
        raise RuntimeError(
            f"GLMY invariant failed: b1 depth3/4={b1}/{b1_4}"
        )
    if b1 > c1:
        raise RuntimeError(f"GLMY invariant failed: b1={b1} > c1={c1}")
    if b2 > c2:
        raise RuntimeError(f"GLMY invariant failed: b2={b2} > c2={c2}")

    return (int(b0_2), int(b1), int(b2), int(c1), int(c2))


def init_worker(digraph_file, vocab_file, distances, controls, control_seed, cache_size):
    global _W_DIGRAPH
    global _W_TOKEN_GROUPS
    global _W_TOKEN_LENGTHS
    global _W_GROUP_SIZES
    global _W_DISTANCES
    global _W_MAX_DISTANCE
    global _W_CONTROL_SEED
    global _W_CONTROLS
    global _cached_glmy5

    _W_DIGRAPH = load_digraph_class(Path(digraph_file))
    tokens = read_vocab(Path(vocab_file))
    _W_TOKEN_GROUPS = make_token_groups(tokens)
    _W_TOKEN_LENGTHS = sorted({len(token) for token in tokens})
    _W_GROUP_SIZES = np.asarray([len(group) for group in _W_TOKEN_GROUPS], dtype=np.int8)
    _W_DISTANCES = tuple(int(x) for x in distances)
    _W_MAX_DISTANCE = max(_W_DISTANCES)
    _W_CONTROL_SEED = int(control_seed)
    _W_CONTROLS = tuple(controls)

    @lru_cache(maxsize=int(cache_size))
    def cached_glmy5(n_vertices, weighted_edges):
        return calculate_glmy5(_W_DIGRAPH, n_vertices, weighted_edges)

    _cached_glmy5 = cached_glmy5


def sequence_glmy5_tensor(sequence):
    """Shape: n_distances x 210 token groups x 5 channels."""
    position_maps = get_position_maps(sequence, _W_TOKEN_LENGTHS)
    tensor = np.zeros(
        (len(_W_DISTANCES), len(_W_TOKEN_GROUPS), len(CHANNELS)),
        dtype=np.float32,
    )

    for group_index, group in enumerate(_W_TOKEN_GROUPS):
        positions, length_by_position = group_positions(position_maps, group)
        raw_nodes = sorted(positions)
        n_vertices = len(raw_nodes)
        if n_vertices == 0:
            continue

        max_edges = []
        for u_index, u in enumerate(raw_nodes):
            lu = length_by_position.get(u, 1)
            for v_index, v in enumerate(raw_nodes):
                if u_index == v_index:
                    continue
                lv = length_by_position.get(v, 1)
                d = adjusted_distance(u, v, lu, lv)
                if d <= _W_MAX_DISTANCE:
                    max_edges.append((u_index, v_index, int(d)))

        for distance_index, threshold in enumerate(_W_DISTANCES):
            weighted_edges = tuple(
                (u, v, edge_distance)
                for u, v, edge_distance in max_edges
                if edge_distance <= threshold
            )
            tensor[distance_index, group_index, :] = _cached_glmy5(
                n_vertices, weighted_edges
            )

    return tensor


def summarize_tokens(tensor, group_sizes):
    pieces = []
    names = []
    for group_size in (1, 2):
        mask = group_sizes == group_size
        sub = tensor[:, mask, :]
        if sub.shape[1] == 0:
            raise RuntimeError(f"No token groups for size {group_size}")

        for descriptor_index, descriptor in enumerate(CHANNELS):
            values = sub[:, :, descriptor_index]
            stats = {
                "mean": np.mean(values, axis=1),
                "std": np.std(values, axis=1, ddof=0),
                "max": np.max(values, axis=1),
            }
            for statistic in TOKEN_SUMMARIES:
                pieces.append(stats[statistic][:, None])
                names.append(f"size{group_size}_{descriptor}_tok{statistic}")

    return np.concatenate(pieces, axis=1).astype(np.float32, copy=False), names


def build_event_features(wt_tensor, mut_tensor, group_sizes):
    block_values = {
        "wt": wt_tensor,
        "mut": mut_tensor,
        "delta": mut_tensor - wt_tensor,
        "absdelta": np.abs(mut_tensor - wt_tensor),
    }
    pieces = []
    names = []
    for block in BLOCKS:
        values, local_names = summarize_tokens(block_values[block], group_sizes)
        pieces.append(values)
        names.extend(f"{block}_{name}" for name in local_names)
    return np.concatenate(pieces, axis=1).astype(np.float32, copy=False), names


class GroupAccumulator:
    """Exact Revision64 aggregation: RC6 x event-feature x mean/std/max."""

    def __init__(self, patients, event_feature_names):
        self.patients = list(patients)
        self.event_feature_names = list(event_feature_names)
        self.patient_index = {patient: i for i, patient in enumerate(self.patients)}
        self.class_index = {mutation_class: i for i, mutation_class in enumerate(RC6_CLASSES)}
        self.n_groups = len(self.patients) * len(RC6_CLASSES)
        self.n_features = len(self.event_feature_names)
        shape = (self.n_groups, self.n_features)
        self.sum = np.zeros(shape, dtype=np.float64)
        self.sumsq = np.zeros(shape, dtype=np.float64)
        self.maximum = np.full(shape, -np.inf, dtype=np.float32)
        self.count = np.zeros(self.n_groups, dtype=np.int64)

    def group_indices(self, patient_barcodes, mutation_classes):
        indices = []
        for patient, mutation_class in zip(patient_barcodes, mutation_classes):
            p = self.patient_index.get(normalize_patient(patient))
            c = self.class_index.get(str(mutation_class))
            indices.append(-1 if p is None or c is None else p * len(RC6_CLASSES) + c)
        return np.asarray(indices, dtype=np.int64)

    def update(self, group_indices, values):
        valid = group_indices >= 0
        if not np.any(valid):
            return 0
        group_indices = group_indices[valid]
        values = values[valid].astype(np.float32, copy=False)
        np.add.at(self.sum, group_indices, values)
        np.add.at(self.sumsq, group_indices, values.astype(np.float64) ** 2)
        np.maximum.at(self.maximum, group_indices, values)
        np.add.at(self.count, group_indices, 1)
        return int(len(group_indices))

    def finalize(self, prefix):
        counts = self.count.astype(np.float64)
        mean = np.divide(
            self.sum,
            counts[:, None],
            out=np.zeros_like(self.sum),
            where=counts[:, None] > 0,
        )
        variance = np.divide(
            self.sumsq,
            counts[:, None],
            out=np.zeros_like(self.sumsq),
            where=counts[:, None] > 0,
        ) - mean ** 2
        variance = np.maximum(variance, 0.0)
        std = np.sqrt(variance)
        maximum = self.maximum.astype(np.float64)
        maximum[~np.isfinite(maximum)] = 0.0
        mutation_stats = {"mean": mean, "std": std, "max": maximum}

        n_patients = len(self.patients)
        n_classes = len(RC6_CLASSES)
        output_columns = []
        output_arrays = []

        for class_index, mutation_class in enumerate(RC6_CLASSES):
            class_name = mutation_class.replace(">", "to")
            row_indices = np.arange(n_patients) * n_classes + class_index
            for event_feature_index, event_feature in enumerate(self.event_feature_names):
                for mutation_stat in MUTATION_SUMMARIES:
                    output_columns.append(
                        f"{prefix}_{class_name}_{event_feature}_mut{mutation_stat}"
                    )
                    output_arrays.append(
                        mutation_stats[mutation_stat][row_indices, event_feature_index][:, None]
                    )

        matrix = np.concatenate(output_arrays, axis=1).astype(np.float32, copy=False)
        return matrix, output_columns


def sequence_control_qc(wt, mut, wt_shuffle, mut_shuffle, left_changed, right_changed):
    result = {
        "sequence_changed": int(wt_shuffle != wt),
        "left_segment_changed": int(left_changed),
        "right_segment_changed": int(right_changed),
        "center3_preserved": int(wt_shuffle[9:12] == wt[9:12] and mut_shuffle[9:12] == mut[9:12]),
        "center5_preserved": int(wt_shuffle[8:13] == wt[8:13] and mut_shuffle[8:13] == mut[8:13]),
        "sbs1536_preserved": int(wt_shuffle[8:13] == wt[8:13] and mut_shuffle[8:13] == mut[8:13]),
        "sbs96_preserved": int(sbs96_key(wt_shuffle, mut_shuffle) == sbs96_key(wt, mut)),
        "wt_mono_preserved": int(mono_counts(wt_shuffle) == mono_counts(wt)),
        "mut_mono_preserved": int(mono_counts(mut_shuffle) == mono_counts(mut)),
        "wt_dinuc_preserved": int(dinuc_counts(wt_shuffle) == dinuc_counts(wt)),
        "mut_dinuc_preserved": int(dinuc_counts(mut_shuffle) == dinuc_counts(mut)),
        "wt_mut_outside_center_match": int(
            wt_shuffle[:10] == mut_shuffle[:10] and wt_shuffle[11:] == mut_shuffle[11:]
        ),
    }
    required = (
        "center3_preserved", "center5_preserved", "sbs96_preserved",
        "sbs1536_preserved", "wt_mono_preserved", "mut_mono_preserved",
        "wt_dinuc_preserved", "mut_dinuc_preserved", "wt_mut_outside_center_match",
    )
    if not all(result[key] == 1 for key in required):
        raise RuntimeError(f"Dinucleotide-control invariant failed: {result}")
    return result


def process_event(task):
    event_id, patient_barcode, wt_context, mut_context, ref, alt = task
    wt, mut, rc6 = canonicalize_context(wt_context, mut_context, ref, alt)

    event_feature_count = len(BLOCKS) * 2 * len(CHANNELS) * len(TOKEN_SUMMARIES)
    feature_cube = np.zeros(
        (len(_W_CONTROLS), len(_W_DISTANCES), event_feature_count),
        dtype=np.float32,
    )

    qc = None

    if "real" in _W_CONTROLS:
        real_wt = sequence_glmy5_tensor(wt)
        real_mut = sequence_glmy5_tensor(mut)
        ci = _W_CONTROLS.index("real")
        for di in range(len(_W_DISTANCES)):
            ef, _ = build_event_features(
                real_wt[di:di + 1], real_mut[di:di + 1], _W_GROUP_SIZES
            )
            feature_cube[ci, di, :] = ef[0]

    if "outer_dinuc_shuffle" in _W_CONTROLS:
        rng = stable_rng(_W_CONTROL_SEED, event_id, "outer_dinuc_shuffle")
        wt_s, mut_s, left_changed, right_changed = make_control_pair(
            wt, mut, "outer_dinuc_shuffle", rng
        )
        qc = sequence_control_qc(wt, mut, wt_s, mut_s, left_changed, right_changed)
        swt = sequence_glmy5_tensor(wt_s)
        smut = sequence_glmy5_tensor(mut_s)
        ci = _W_CONTROLS.index("outer_dinuc_shuffle")
        for di in range(len(_W_DISTANCES)):
            ef, _ = build_event_features(swt[di:di + 1], smut[di:di + 1], _W_GROUP_SIZES)
            feature_cube[ci, di, :] = ef[0]

    return patient_barcode, rc6, feature_cube, qc


def flush_results(buffer, controls, distances, accumulators, qc_sums):
    if not buffer:
        return 0

    patients = [row[0] for row in buffer]
    classes = [row[1] for row in buffer]
    reference = accumulators[(controls[0], distances[0])]
    group_indices = reference.group_indices(patients, classes)
    cube = np.stack([row[2] for row in buffer], axis=0)

    for ci, control in enumerate(controls):
        for di, distance in enumerate(distances):
            accumulators[(control, distance)].update(group_indices, cube[:, ci, di, :])

    if "outer_dinuc_shuffle" in controls:
        entry = qc_sums["outer_dinuc_shuffle"]
        for row in buffer:
            qc = row[3]
            if qc is None:
                raise RuntimeError("Missing shuffle QC record")
            entry["n"] += 1
            for key, value in qc.items():
                entry[key] = entry.get(key, 0.0) + float(value)

    count = len(buffer)
    buffer.clear()
    return count


def descriptor_from_feature_name(feature_name):
    for descriptor in CHANNELS:
        if f"_{descriptor}_" in feature_name:
            return descriptor
    return ""


def self_test(digraph_file):
    Digraph = load_digraph_class(digraph_file)
    sequence = "AATGATGCACATCATGGTGGC"
    group = ("A", "GC")
    threshold = 4
    token_lengths = [1, 2]
    position_maps = get_position_maps(sequence, token_lengths)
    positions, length_by_position = group_positions(position_maps, group)
    nodes = sorted(positions)
    edges = []
    for ui, u in enumerate(nodes):
        for vi, v in enumerate(nodes):
            if ui == vi:
                continue
            d = adjusted_distance(
                u, v, length_by_position.get(u, 1), length_by_position.get(v, 1)
            )
            if d <= threshold:
                edges.append((ui, vi, d))
    observed = calculate_glmy5(Digraph, len(nodes), tuple(edges))
    expected = (2, 0, 1, 12, 19)
    print(f"[SELF-TEST] sequence={sequence} group=A__GC d={threshold}")
    print(f"[SELF-TEST] channels={CHANNELS}")
    print(f"[SELF-TEST] observed={observed} expected={expected}")
    if observed != expected:
        raise SystemExit("GLMY5 self-test failed; Digraph_Hu.py may have changed")
    print("[SELF-TEST] PASS")


def main():
    parser = argparse.ArgumentParser(
        description="Fresh submission-grade five-channel GLMY patient feature build."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--baseline-dir", default="results/revision_baselines_60_dedup_v2")
    parser.add_argument("--digraph-file", required=True)
    parser.add_argument("--vocab-file", default="data/features/vocab_dna_1to2mer.tsv")
    parser.add_argument("--cancers", default=",".join(CANCERS))
    parser.add_argument("--distances", default="2,3,4,5")
    parser.add_argument("--controls", default="real,outer_dinuc_shuffle")
    parser.add_argument("--control-seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--pool-chunksize", type=int, default=4)
    parser.add_argument("--aggregate-batch-size", type=int, default=128)
    parser.add_argument("--glmy-cache-size", type=int, default=50000)
    parser.add_argument(
        "--max-events-per-cancer",
        type=int,
        default=0,
        help="Pilot only. 0 = all exact-deduplicated model-cohort events.",
    )
    parser.add_argument("--progress-every", type=int, default=5000)
    parser.add_argument("--outdir", default="results/revision79A2_glmy5_fresh_full")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    started = time.time()
    root = Path(args.project_root).expanduser().resolve()
    baseline_dir = resolve_path(root, args.baseline_dir)
    digraph_file = resolve_path(root, args.digraph_file)
    vocab_file = resolve_path(root, args.vocab_file)
    outdir = resolve_path(root, args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.self_test:
        self_test(digraph_file)
        return 0

    cancers = tuple(x.upper() for x in parse_csv(args.cancers))
    unknown_cancers = sorted(set(cancers) - set(CANCERS))
    if unknown_cancers:
        raise SystemExit(f"Unsupported cancers: {unknown_cancers}")

    distances = tuple(sorted({int(x) for x in parse_csv(args.distances)}))
    if not distances or min(distances) < 1:
        raise SystemExit("Distances must be positive integers")

    controls = tuple(parse_csv(args.controls))
    unsupported_controls = sorted(set(controls) - set(CONTROLS))
    if unsupported_controls:
        raise SystemExit(f"Unsupported controls: {unsupported_controls}")
    if not controls:
        raise SystemExit("No controls requested")

    tokens = read_vocab(vocab_file)
    token_groups = make_token_groups(tokens)
    group_sizes = np.asarray([len(group) for group in token_groups], dtype=np.int8)

    dummy = np.zeros((1, len(token_groups), len(CHANNELS)), dtype=np.float32)
    _, event_feature_names = build_event_features(dummy, dummy, group_sizes)
    expected_event_features = len(BLOCKS) * 2 * len(CHANNELS) * len(TOKEN_SUMMARIES)
    expected_patient_features = (
        len(RC6_CLASSES) * expected_event_features * len(MUTATION_SUMMARIES)
    )
    if len(event_feature_names) != 120 or expected_event_features != 120:
        raise SystemExit(
            f"Expected 120 event features for five channels, got {len(event_feature_names)}"
        )
    if expected_patient_features != 2160:
        raise SystemExit(f"Expected 2160 patient features, got {expected_patient_features}")

    workers = int(args.workers) if int(args.workers) > 0 else max(1, os.cpu_count() or 1)

    build_qc_rows = []
    dedup_qc_rows = []
    sequence_qc_rows = []
    manifest_rows = []

    print("=" * 92)
    print("Revision79A2 fresh five-channel GLMY build")
    print(f"channels={CHANNELS}")
    print(f"event_features={expected_event_features}; patient_features/d={expected_patient_features}")
    print(f"cancers={cancers}; controls={controls}; distances={distances}")
    print(f"Digraph_Hu SHA256={sha256_file(digraph_file)}")
    print("=" * 92)

    for cancer in cancers:
        cancer_started = time.time()
        print(f"\n{'=' * 36} {cancer} {'=' * 36}")

        context_path = root / "data" / "features" / f"TCGA_{cancer}_snv_context_w10.tsv"
        model_path = baseline_dir / "per_cancer" / f"{cancer}_mutation_spectrum_model_matrix.tsv"
        if not context_path.exists():
            raise SystemExit(f"Missing context: {context_path}")
        if not model_path.exists():
            raise SystemExit(f"Missing model matrix: {model_path}")

        context = pd.read_csv(context_path, sep="\t", low_memory=False)
        model = pd.read_csv(model_path, sep="\t", low_memory=False)
        context["patient_barcode"] = context["patient_barcode"].astype(str).str.slice(0, 12)
        model["patient_barcode"] = model["patient_barcode"].astype(str).str.slice(0, 12)

        if model["patient_barcode"].duplicated().any():
            raise SystemExit(f"{cancer}: duplicate model patient rows")

        model_patients = model["patient_barcode"].tolist()
        model_patient_set = set(model_patients)
        duplicate_key = resolve_duplicate_key(context)
        if not duplicate_key:
            raise SystemExit(f"{cancer}: could not resolve exact-event duplicate key")

        dedup_keep = ~context.duplicated(duplicate_key, keep="first")
        model_keep = context["patient_barcode"].isin(model_patient_set)
        final_keep = dedup_keep & model_keep
        selected_indices = np.flatnonzero(final_keep.to_numpy())
        full_model_events = int(len(selected_indices))

        pilot_mode = int(args.max_events_per_cancer) > 0
        if pilot_mode and len(selected_indices) > int(args.max_events_per_cancer):
            rng = np.random.default_rng(int(args.control_seed) + 10007 * CANCERS.index(cancer))
            selected_indices = np.sort(
                rng.choice(selected_indices, size=int(args.max_events_per_cancer), replace=False)
            )

        dedup_qc_rows.append({
            "cancer_type": cancer,
            "context_rows_raw": int(len(context)),
            "context_rows_deduplicated": int(dedup_keep.sum()),
            "exact_duplicates_removed": int(len(context) - dedup_keep.sum()),
            "full_model_patient_events": full_model_events,
            "events_selected_this_run": int(len(selected_indices)),
            "model_patients": int(len(model_patients)),
            "pilot_mode": int(pilot_mode),
            "duplicate_key": ",".join(duplicate_key),
        })

        accumulators = {
            (control, distance): GroupAccumulator(model_patients, event_feature_names)
            for control in controls
            for distance in distances
        }
        qc_sums = {"outer_dinuc_shuffle": {"n": 0}}

        tasks = (
            (
                f"{cancer}_{int(index)}",
                normalize_patient(context.iloc[index]["patient_barcode"]),
                context.iloc[index]["wt_context"],
                context.iloc[index]["mut_context"],
                context.iloc[index]["ref"],
                context.iloc[index]["alt"],
            )
            for index in selected_indices
        )

        processed = 0
        buffer = []
        ctx = mp.get_context("fork")
        print(
            f"{cancer}: events={len(selected_indices):,}, patients={len(model_patients)}, "
            f"workers={workers}, controls={controls}, distances={distances}"
        )

        with ctx.Pool(
            processes=workers,
            initializer=init_worker,
            initargs=(
                str(digraph_file), str(vocab_file), distances, controls,
                int(args.control_seed), int(args.glmy_cache_size),
            ),
            maxtasksperchild=1000,
        ) as pool:
            for result in pool.imap(process_event, tasks, chunksize=int(args.pool_chunksize)):
                buffer.append(result)
                if len(buffer) >= int(args.aggregate_batch_size):
                    processed += flush_results(buffer, controls, distances, accumulators, qc_sums)
                    if (
                        int(args.progress_every) > 0
                        and processed % int(args.progress_every) < int(args.aggregate_batch_size)
                    ):
                        print(f"[{cancer}] processed {processed:,}/{len(selected_indices):,}")
            processed += flush_results(buffer, controls, distances, accumulators, qc_sums)

        if processed != len(selected_indices):
            raise SystemExit(f"{cancer}: processed {processed}, expected {len(selected_indices)}")

        if "outer_dinuc_shuffle" in controls:
            entry = qc_sums["outer_dinuc_shuffle"]
            if entry["n"] != processed:
                raise SystemExit(f"{cancer}: sequence QC count {entry['n']} != {processed}")
            qc_row = {
                "cancer_type": cancer,
                "control": "outer_dinuc_shuffle",
                "n_events": int(entry["n"]),
            }
            for key, total in entry.items():
                if key != "n":
                    qc_row[key] = float(total) / float(entry["n"])
            sequence_qc_rows.append(qc_row)
            required = (
                "center3_preserved", "center5_preserved", "sbs96_preserved",
                "sbs1536_preserved", "wt_mono_preserved", "mut_mono_preserved",
                "wt_dinuc_preserved", "mut_dinuc_preserved", "wt_mut_outside_center_match",
            )
            for column in required:
                if not np.isclose(qc_row[column], 1.0, atol=0.0, rtol=0.0):
                    raise SystemExit(f"{cancer}: shuffle preservation failed: {column}={qc_row[column]}")

        for control in controls:
            for distance in distances:
                matrix, columns = accumulators[(control, distance)].finalize(
                    f"glmy5_d{distance}"
                )
                if len(columns) != 2160:
                    raise SystemExit(
                        f"{cancer} {control} d{distance}: expected 2160 features, got {len(columns)}"
                    )
                if not np.isfinite(matrix).all():
                    raise SystemExit(f"{cancer} {control} d{distance}: non-finite matrix")

                base_cols = ["patient_barcode", "cancer_type"] + [
                    x for x in LABEL_COLS if x in model.columns
                ]
                output = pd.concat(
                    [
                        model[base_cols].reset_index(drop=True),
                        pd.DataFrame(matrix, columns=columns),
                    ],
                    axis=1,
                )

                if control == "real":
                    filename = f"{cancer}_glmy5_fresh_d{distance}.tsv.gz"
                else:
                    filename = f"{cancer}_glmy5_fresh_{control}_d{distance}.tsv.gz"
                output_dir = outdir / "per_control" / control / f"d{distance}"
                output_dir.mkdir(parents=True, exist_ok=True)
                output_path = output_dir / filename
                output.to_csv(output_path, sep="\t", index=False, compression="gzip")

                build_qc_rows.append({
                    "cancer_type": cancer,
                    "control": control,
                    "distance": int(distance),
                    "events_processed": int(processed),
                    "model_patients": int(len(model_patients)),
                    "event_features": int(len(event_feature_names)),
                    "patient_features": int(len(columns)),
                    "output_rows": int(len(output)),
                    "pilot_mode": int(pilot_mode),
                    "all_glmy_invariants_fail_closed": 1,
                    "output_path": str(output_path),
                })

                for feature in columns:
                    manifest_rows.append({
                        "control": control,
                        "distance": int(distance),
                        "feature": feature,
                        "descriptor": descriptor_from_feature_name(feature),
                    })

                print(
                    f"[OK] {cancer} {control} d{distance}: "
                    f"patients={len(output)} features={len(columns)}"
                )

        print(f"{cancer}: completed in {time.time() - cancer_started:.1f} s")

    build_qc = pd.DataFrame(build_qc_rows)
    dedup_qc = pd.DataFrame(dedup_qc_rows)
    sequence_qc = pd.DataFrame(sequence_qc_rows)
    manifest_frame = (
        pd.DataFrame(manifest_rows)
        .drop_duplicates()
        .sort_values(["control", "distance", "feature"])
        .reset_index(drop=True)
    )

    build_qc.to_csv(outdir / "build_qc.tsv", sep="\t", index=False)
    dedup_qc.to_csv(outdir / "dedup_qc.tsv", sep="\t", index=False)
    sequence_qc.to_csv(outdir / "sequence_control_qc.tsv", sep="\t", index=False)

    sequence_summary_rows = []
    if not sequence_qc.empty:
        for control, group in sequence_qc.groupby("control", sort=False):
            weights = group["n_events"].to_numpy(dtype=float)
            row = {"control": control, "n_events": int(weights.sum())}
            for column in sequence_qc.columns:
                if column in ("cancer_type", "control", "n_events"):
                    continue
                row[column] = float(np.average(group[column], weights=weights))
            sequence_summary_rows.append(row)
    pd.DataFrame(sequence_summary_rows).to_csv(
        outdir / "sequence_control_qc_summary.tsv", sep="\t", index=False
    )

    manifest_frame.to_csv(outdir / "feature_manifest.tsv", sep="\t", index=False)
    dimension_summary = (
        manifest_frame.groupby(["control", "distance", "descriptor"], sort=False)
        .size()
        .reset_index(name="n_features")
    )
    dimension_summary.to_csv(outdir / "feature_dimension_summary.tsv", sep="\t", index=False)

    definition = pd.DataFrame([
        {"channel": "b0", "type": "stable_homology", "definition": "infinite H0 at path depth 2", "features_per_distance": 432},
        {"channel": "b1", "type": "stable_homology", "definition": "infinite H1 at path depth 3", "features_per_distance": 432},
        {"channel": "b2", "type": "stable_homology", "definition": "infinite H2 at path depth 4", "features_per_distance": 432},
        {"channel": "c1", "type": "transient_path_cycle", "definition": "top H1 at path depth 2", "features_per_distance": 432},
        {"channel": "c2", "type": "transient_path_cycle", "definition": "top H2 at path depth 3", "features_per_distance": 432},
    ])
    definition.to_csv(outdir / "glmy5_feature_definition.tsv", sep="\t", index=False)

    runtime = float(time.time() - started)
    script_path = Path(__file__).resolve()
    run_manifest = {
        "revision": "79A2",
        "method": "fresh five-channel GLMY stable-homology + transient-cycle representation",
        "canonical_channel_order": list(CHANNELS),
        "definitions": {
            "b0": "infinite H0 at depth2",
            "b1": "infinite H1 at depth3",
            "b2": "infinite H2 at depth4",
            "c1": "top-dimensional H1 at depth2",
            "c2": "top-dimensional H2 at depth3",
        },
        "fresh_recomputation": True,
        "reuses_revision74_stable_matrices": False,
        "reuses_revision79A_transient_matrices": False,
        "project_root": str(root),
        "baseline_dir": str(baseline_dir),
        "digraph_file": str(digraph_file),
        "vocab_file": str(vocab_file),
        "script_sha256": sha256_file(script_path),
        "digraph_sha256": sha256_file(digraph_file),
        "vocab_sha256": sha256_file(vocab_file),
        "cancers": list(cancers),
        "distances": list(distances),
        "controls": list(controls),
        "control_seed": int(args.control_seed),
        "model_patient_total_expected": 2003,
        "token_count": len(tokens),
        "token_group_count": len(token_groups),
        "event_features": 120,
        "patient_features_per_distance": 2160,
        "features_per_channel_per_distance": 432,
        "aggregation": "RC6 x WT/MUT/delta/absdelta x token-size x channel x token mean/std/max x mutation mean/std/max",
        "invariants": [
            "c1_depth2_equals_lenH1",
            "c2_depth3_equals_lenH2",
            "b0_stable_depth2_3_4",
            "b1_stable_depth3_4",
            "b1_le_c1",
            "b2_le_c2",
        ],
        "workers": workers,
        "max_events_per_cancer": int(args.max_events_per_cancer),
        "pilot_mode": bool(int(args.max_events_per_cancer) > 0),
        "runtime_seconds": runtime,
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    readme = f"""# Revision 79A2 fresh GLMY5 build

This directory is a full fresh recomputation of the five frozen predictive GLMY
channels. No historical GLMY patient matrix is reused.

Canonical channels: `{CHANNELS}`

- stable homology: `b0,b1,b2`
- transient path-cycle descriptors: `c1,c2`

Features per channel/distance: **432**

Features per distance: **2,160**

Distances: `{distances}`

Controls: `{controls}`

Cancers: `{cancers}`

Pilot mode: `{bool(int(args.max_events_per_cancer) > 0)}`

The program fail-closes if any audited GLMY depth invariant is violated. The
fresh matrices produced here should be used as the authoritative source for
Rev79 combination selection and for future code-submission packaging.
"""
    (outdir / "README_revision79A2.md").write_text(readme, encoding="utf-8")

    print("\n" + "=" * 92)
    print("Revision 79A2 completed")
    print(f"Output: {outdir}")
    print(f"Runtime seconds: {runtime:.3f}")
    print("Authoritative fresh channels: b0,b1,b2,c1,c2")
    print("=" * 92)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
