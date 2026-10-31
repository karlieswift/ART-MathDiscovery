#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute graph irregularity features for each BPE/k-mer token in each sequence.

Graph construction follows hu_glmy_b01.py / hu_glmy_c012.py style:
  1. For each token/k-mer, find all occurrence start positions in the sequence.
  2. Treat occurrences as graph vertices.
  3. Connect two occurrence vertices when their adjusted distance is <= --max-distance.

The graph is treated as a simple undirected graph.

Feature names:
  V1 = var(G)
  V2 = s(G)
  V3 = irrt(G) = sum_{i<=j} |d_i - d_j|
  V3SQ = irrt2(G) = sum_{i<=j} (d_i - d_j)^2
  V3SIGN = sorted signed irrt(G) = sum_{i<=j} (d_(i) - d_(j)), degrees sorted descending
  V4 = PI1(G)
  V5 = irrent(G) = ln(n) - H(d_i / 2m)
  V6 = IA(G), adjacency shifted-eigenvalue dot degree sequence
  V7 = IL(G), Laplacian shifted-eigenvalue dot degree sequence
  V8 = IQ(G), signless-Laplacian shifted-eigenvalue dot degree sequence
  BETA_a<alpha> = generalized first Betti number beta_1^alpha(G)
  LGBETA = beta_1(L(G)), first Betti number of the line graph of G
  L2GBETA = beta_1(L2(G)), first Betti number of the second iterated line graph of G

Examples:
  # Original behavior
  python hu_glmy_irregular_with_draw_v2_features.py \
    --input H4/dev.csv \
    --output H4/dev_irregular.csv \
    --vocab-file vocab_dna_1to5mer.tsv \
    --max-distance 10 \
    --features V1,V2,V3,V4

  # Recommended compact H4 test: old best + entropy + selected Betti alphas
  python hu_glmy_irregular_with_draw_v2_features.py \
    --input H4/dev.csv \
    --output H4/dev_irregular_plus_beta.csv \
    --vocab-file vocab_dna_1to5mer.tsv \
    --max-distance 10 \
    --features V1,V3,V4,V5,BETA \
    --betti-alphas 1.5,2,2.71,3,3.14

  # Add spectral features too. This is much slower.
  python hu_glmy_irregular_with_draw_v2_features.py \
    --input H4/dev.csv \
    --output H4/dev_irregular_plus_spectral.csv \
    --vocab-file vocab_dna_1to5mer.tsv \
    --max-distance 10 \
    --features V1,V3,V4,V5,V6,V7,V8,BETA \
    --spectral-max-vertices 200
"""

import argparse
import csv
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations

import pandas as pd
from tqdm import tqdm

csv.field_size_limit(sys.maxsize)


# ---------------- CLI ----------------
p = argparse.ArgumentParser(
    description="Compute degree, entropy, spectral, and generalized-Betti graph irregularity features."
)
p.add_argument("--input", required=True, help="Input CSV path; must contain sequence and label columns")
p.add_argument("--output", required=True, help="Output CSV path")
p.add_argument("--max-distance", type=int, default=2, required=True, help="Max graph edge distance")
p.add_argument("--vocab-file", default="bpe_vocab_1500.txt", help="TSV vocab file with Token column")
p.add_argument(
    "--combo-k", "--token-combo-k",
    dest="combo_k",
    type=int,
    default=1,
    help=(
        "Token combination size. 1 = original single-token graphs; "
        "2 = unordered token-pair graphs; 3 = unordered token-triple graphs."
    ),
)
p.add_argument(
    "--combo-include-smaller",
    action="store_true",
    help="If set, include all combination sizes 1..combo-k. Example combo-k=2 gives singles+pairs.",
)
p.add_argument("--workers", type=int, default=32, help="Number of parallel workers; default: 32")
p.add_argument(
    "--irrt-mode",
    choices=["abs", "square", "sign", "both", "all"],
    default="abs",
    help=(
        "How to output V3 when --features includes V3. "
        "abs = original irrt(G)=sum |d_i-d_j|; "
        "square = output V3SQ=sum (d_i-d_j)^2 instead of V3; "
        "sign = output V3SIGN=sum sorted(d_i)-sorted(d_j) for i<=j; "
        "both = output both V3 and V3SQ, preserving old behavior; "
        "all = output V3, V3SQ, and V3SIGN. Default: abs."
    ),
)

p.add_argument(
    "--features",
    default="V1,V2,V3,V4",
    help=(
        "Comma/plus/space separated features to output. "
        "Allowed: V1,V2,V3,V3SQ,V3SIGN,V4,V5,V6,V7,V8,BETA,LGBETA,L2GBETA,all. "
        "Aliases: var,s,irrt,irrt2,irrtsq,irrtsign,PI1,irrent,IA,IL,IQ,beta,linebeta,line2beta. "
        "Compact forms like V1V3V4 are supported. "
        "Default: V1,V2,V3,V4."
    ),
)
p.add_argument(
    "--betti-alphas",
    default="1.5,2,2.71,3,3.14",
    help="Comma-separated alpha values for BETA features when --features includes BETA. Default: 1.5,2,2.71,3,3.14",
)
p.add_argument(
    "--spectral-max-vertices",
    type=int,
    default=200,
    help=(
        "Max graph vertices for exact eigenvalue features V6/V7/V8. "
        "If a token graph is larger, spectral values are written as 0. "
        "Use 0 or negative for no limit. Default: 200."
    ),
)

# Drawing controls. Disabled by default to keep full-batch runs fast and clean.
p.add_argument("--draw-irregular", action="store_true", help="Draw token graphs when G is irregular.")
p.add_argument("--draw-dir", default="pic_irregular", help="Directory for graph PNG/CSV files")
p.add_argument("--draw-max-per-seq", type=int, default=0, help="Max irregular graphs to draw per sequence. 0 means no limit.")
p.add_argument("--draw-min-vertices", type=int, default=2, help="Only draw graphs with at least this many vertices.")
p.add_argument("--draw-max-vertices", type=int, default=80, help="Skip drawing graphs above this many vertices.")
p.add_argument("--irregular-eps", type=float, default=1e-12, help="Tolerance for deciding irregularity from var(G).")
args = p.parse_args()

input_csv = args.input
output_csv = args.output
MAX_DISTANCE = args.max_distance


# ---------------- Feature selection ----------------
BASE_FEATURES = ["V1", "V2", "V3", "V3SQ", "V3SIGN", "V4", "V5", "V6", "V7", "V8"]
FEATURE_ALIASES = {
    "V1": "V1", "VAR": "V1", "VARG": "V1", "VAR(G)": "V1",
    "V2": "V2", "S": "V2", "SG": "V2", "S(G)": "V2",
    "V3": "V3", "IRRT": "V3", "IRRTG": "V3", "IRRT(G)": "V3",
    "V3SQ": "V3SQ", "V3SQUARE": "V3SQ", "IRRT2": "V3SQ", "IRRTSQ": "V3SQ", "IRRTSQUARE": "V3SQ", "SQUAREIRRT": "V3SQ",
    "V3SIGN": "V3SIGN", "V3SIGNED": "V3SIGN", "IRRTSIGN": "V3SIGN", "IRRTSIGNED": "V3SIGN", "SIGNEDIRRT": "V3SIGN",
    "V4": "V4", "PI1": "V4", "PI1G": "V4", "PI1(G)": "V4",
    "V5": "V5", "IRRENT": "V5", "IRRENTG": "V5", "IRRENT(G)": "V5", "ENTROPY": "V5",
    "V6": "V6", "IA": "V6", "IAG": "V6", "IA(G)": "V6",
    "V7": "V7", "IL": "V7", "ILG": "V7", "IL(G)": "V7",
    "V8": "V8", "IQ": "V8", "IQG": "V8", "IQ(G)": "V8",
    "BETA": "BETA", "BETTI": "BETA", "BETTI1": "BETA", "BETAG": "BETA",
    "LGBETA": "LGBETA", "LINEBETA": "LGBETA", "BETALINE": "LGBETA",
    "BETA_LINE": "LGBETA", "BETA_LG": "LGBETA", "BETA1LG": "LGBETA",
    "L2GBETA": "L2GBETA", "LGBETA2": "L2GBETA", "LINE2BETA": "L2GBETA",
    "SECONDLINEBETA": "L2GBETA", "BETALINE2": "L2GBETA", "BETA2LINE": "L2GBETA",
    "BETA_L2G": "L2GBETA", "BETA1L2G": "L2GBETA", "BETA1L2": "L2GBETA",
}


def parse_alpha_list(x):
    alphas = []
    seen = set()
    for part in re.split(r"[,;+|\s]+", str(x).strip()):
        if not part:
            continue
        a = float(part)
        if a < 1:
            raise ValueError(f"Betti alpha must be >= 1, got {a}")
        key = f"{a:.12g}"
        if key not in seen:
            alphas.append(a)
            seen.add(key)
    if not alphas:
        raise ValueError("--betti-alphas produced no valid alpha values")
    return alphas


def compact_feature_string(raw):
    return str(raw).upper().replace("_", "").replace("-", "").replace(" ", "")


def parse_selected_features(x):
    raw = str(x).strip()
    if raw == "" or raw.lower() == "all":
        return BASE_FEATURES + ["BETA", "LGBETA"]

    compact = compact_feature_string(raw)
    compact_map = {
        "V1V2V3V4": ["V1", "V2", "V3", "V4"],
        "V1234": ["V1", "V2", "V3", "V4"],
        "V1V3V4": ["V1", "V3", "V4"],
        "V134": ["V1", "V3", "V4"],
        "V1V3V4V5": ["V1", "V3", "V4", "V5"],
        "V1345": ["V1", "V3", "V4", "V5"],
        "V1V3V4V5BETA": ["V1", "V3", "V4", "V5", "BETA"],
        "V1345BETA": ["V1", "V3", "V4", "V5", "BETA"],
        "V1V2V3V4V5V6V7V8BETALGBETA": ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "BETA", "LGBETA"],
        "V12345678BETALGBETA": ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "BETA", "LGBETA"],
    }
    if compact in compact_map:
        return compact_map[compact]

    parts = re.split(r"[,+;|\s]+", raw)
    selected = []
    seen = set()
    for part in parts:
        if not part:
            continue
        key = part.strip().upper().replace("_", "").replace("-", "")
        feat = FEATURE_ALIASES.get(key)
        if feat is None:
            raise ValueError(
                f"Unknown feature '{part}'. Allowed: V1,V2,V3,V3SQ,V3SIGN,V4,V5,V6,V7,V8,BETA,LGBETA,L2GBETA "
                "or aliases var,s,irrt,irrt2,irrtsq,irrtsign,PI1,irrent,IA,IL,IQ,beta,linebeta,line2beta"
            )
        if feat not in seen:
            selected.append(feat)
            seen.add(feat)
    if not selected:
        raise ValueError("No valid features selected.")
    return selected


SELECTED_FEATURES = parse_selected_features(args.features)
BETTI_ALPHAS = parse_alpha_list(args.betti_alphas)
NEED_SPECTRAL = any(f in SELECTED_FEATURES for f in ("V6", "V7", "V8"))
NEED_BETTI = "BETA" in SELECTED_FEATURES
NEED_LGBETA = "LGBETA" in SELECTED_FEATURES
NEED_L2GBETA = "L2GBETA" in SELECTED_FEATURES


def alpha_name(a):
    # Make stable, column-safe names: 1.5 -> 1p5, 2 -> 2, 3.14 -> 3p14
    s = f"{a:.12g}"
    return s.replace(".", "p").replace("-", "m")


# ---------------- Load vocab ----------------
def load_token_list(vocab_file):
    """Load Token column from TSV vocabulary; keep tokens with length >= 1."""
    try:
        vocab_df = pd.read_csv(vocab_file, sep="\t", usecols=["Token"])
    except Exception as e:
        sys.exit(f"Error reading vocab file: {e}")

    tokens = sorted(tok for tok in vocab_df["Token"].astype(str) if len(tok) >= 1)
    if not tokens:
        sys.exit("Error: no valid tokens found in vocab file.")
    return tokens


KMER_LIST = load_token_list(args.vocab_file)


def make_token_groups(tokens, combo_k, include_smaller=False):
    if combo_k < 1:
        raise ValueError("--combo-k must be >= 1")
    if combo_k > len(tokens):
        raise ValueError(f"--combo-k={combo_k} exceeds vocab size {len(tokens)}")
    sizes = range(1, combo_k + 1) if include_smaller else [combo_k]
    groups = []
    for r in sizes:
        groups.extend(tuple(x) for x in combinations(tokens, r))
    return groups


def token_group_name(group):
    return group[0] if len(group) == 1 else "__".join(group)


TOKEN_GROUPS = make_token_groups(KMER_LIST, args.combo_k, args.combo_include_smaller)
TOKEN_LENGTHS = sorted(set(len(tok) for group in TOKEN_GROUPS for tok in group))


# ---------------- Graph helpers ----------------
def get_positions(seq, k):
    """Return dict: token -> 1-based occurrence start positions."""
    pos = {}
    if k <= 0 or len(seq) < k:
        return pos
    for i in range(len(seq) - k + 1):
        km = seq[i : i + k]
        pos.setdefault(km, []).append(i + 1)
    return pos


def get_group_positions(pos_maps, group):
    """Return sorted unique start positions and max token length observed at each position."""
    length_by_pos = {}
    for tok in group:
        L = len(tok)
        for p0 in pos_maps.get(L, {}).get(tok, []):
            old = length_by_pos.get(p0, 0)
            if L > old:
                length_by_pos[p0] = L
    return sorted(length_by_pos), length_by_pos


def adjusted_distance(u, v, lu, lv):
    """Generalizes the original same-length rule using min token length; identical when lu==lv."""
    d = abs(v - u)
    overlap_len = min(lu, lv)
    if d > overlap_len - 1:
        d = d - overlap_len + 1
    return d


def build_simple_graph(positions, length_by_pos=None, k_len=None):
    """
    Build an undirected simple graph from occurrence start positions.

    For combo-k=1, all vertices have the same token length and this reproduces
    the old distance rule exactly. For mixed-length token combinations, each
    vertex stores the max matched token length at that start position and the
    pairwise adjusted distance uses min(length_u, length_v).
    """
    raw_nodes = sorted(set(positions))
    n = len(raw_nodes)
    edges = []

    if length_by_pos is None:
        if k_len is None:
            k_len = 1
        length_by_pos = {p: k_len for p in raw_nodes}

    for i in range(n):
        u = raw_nodes[i]
        lu = length_by_pos.get(u, k_len or 1)
        for j in range(i + 1, n):
            v = raw_nodes[j]
            lv = length_by_pos.get(v, k_len or 1)
            d = adjusted_distance(u, v, lu, lv)
            if d <= MAX_DISTANCE:
                edges.append((i, j, d))

    return raw_nodes, edges

def component_count(n, edges):
    if n <= 0:
        return 0
    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    for u, v, _w in edges:
        union(u, v)
    return len({find(i) for i in range(n)})


def degree_features(n, edges):
    """Compute V1-V5 and degrees."""
    if n <= 0:
        return {"V1": 0.0, "V2": 0.0, "V3": 0.0, "V3SQ": 0.0, "V3SIGN": 0.0, "V4": 0.0, "V5": 0.0}, []

    deg = [0] * n
    for u, v, _weight in edges:
        deg[u] += 1
        deg[v] += 1

    m = len(edges)
    avg_deg = (2.0 * m) / n

    var_g = sum((d - avg_deg) ** 2 for d in deg) / n
    s_g = sum(abs(d - avg_deg) for d in deg)

    irrt_g = 0.0
    irrt_sq_g = 0.0
    for i in range(n):
        di = deg[i]
        for j in range(i, n):
            diff = di - deg[j]
            irrt_g += abs(diff)
            irrt_sq_g += diff * diff

    # Signed no-absolute variant. To keep this a graph invariant, sort degrees
    # in non-increasing order first; then d_i - d_j >= 0 for i <= j.
    deg_sorted = sorted(deg, reverse=True)
    irrt_sign_g = 0.0
    for i in range(n):
        di = deg_sorted[i]
        for j in range(i, n):
            irrt_sign_g += di - deg_sorted[j]

    denom = 2.0 * m + n
    if denom <= 0:
        pi1_g = 0.0
    else:
        log_term = n * math.log(n / denom) + sum(math.log(d + 1.0) for d in deg)
        try:
            pi1_g = 1.0 - math.exp(log_term)
        except OverflowError:
            pi1_g = float("-inf")

    # irrent(G) = ln(n) - entropy(p), p_i = d_i/(2m). Skip p_i=0.
    # For m=0, degree distribution is undefined; use 0 as a neutral numeric value.
    if n <= 1 or m <= 0:
        irrent_g = 0.0
    else:
        total_deg = 2.0 * m
        entropy = 0.0
        for d in deg:
            if d > 0:
                p = d / total_deg
                entropy -= p * math.log(p)
        irrent_g = math.log(n) + entropy

    return {"V1": var_g, "V2": s_g, "V3": irrt_g, "V3SQ": irrt_sq_g, "V3SIGN": irrt_sign_g, "V4": pi1_g, "V5": irrent_g}, deg


def betti_alpha_features(n, edges, degrees, alphas):
    """Compute generalized beta_1^alpha(G) = c + 1/2 sum d_i^alpha - sum d_i^(alpha-1)."""
    c = component_count(n, edges)
    out = {}
    for a in alphas:
        if abs(a - 1.0) < 1e-12:
            # d^0 is interpreted as 1 for every vertex to recover beta_1 = m - n + c.
            sum_prev = float(n)
        else:
            sum_prev = sum((float(d) ** (a - 1.0)) for d in degrees)
        val = c + 0.5 * sum((float(d) ** a) for d in degrees) - sum_prev
        out[f"BETA_a{alpha_name(a)}"] = val
    return out


def line_graph_betti_feature(n, edges, degrees):
    """Compute beta_1(L(G)) = c - n0 + 1/2 * sum_i d_i * (d_i - 2).

    Here c is the number of connected components of G and n0 is the number
    of isolated vertices of G. Isolated vertices of G create no vertices in L(G).
    """
    if n <= 0:
        return {"LGBETA": 0.0}
    c = component_count(n, edges)
    n0 = sum(1 for d in degrees if d == 0)
    val = (c - n0) + 0.5 * sum(float(d) * (float(d) - 2.0) for d in degrees)
    return {"LGBETA": val}


def second_line_graph_betti_feature(n, edges, degrees):
    """Compute beta_1(L2(G)).

    Formula from Irregularity_v4:
      beta_1(L2(G)) = c - n0 - m0
        + 1/2 * sum_v (d_v^3 - 6 d_v^2 + 4 d_v)
        + sum_{uv in E(G)} d_u d_v

    where c is the number of connected components of G, n0 is the number of
    isolated vertices, and m0 is the number of isolated K2 edge components.
    In a simple graph, an isolated K2 component is exactly an edge whose two
    endpoints both have degree 1.
    """
    if n <= 0:
        return {"L2GBETA": 0.0}
    c = component_count(n, edges)
    n0 = sum(1 for d in degrees if d == 0)
    m0 = sum(1 for u, v, _w in edges if degrees[u] == 1 and degrees[v] == 1)
    vertex_term = 0.5 * sum(
        (float(d) ** 3) - 6.0 * (float(d) ** 2) + 4.0 * float(d)
        for d in degrees
    )
    edge_term = sum(float(degrees[u]) * float(degrees[v]) for u, v, _w in edges)
    val = (c - n0 - m0) + vertex_term + edge_term
    return {"L2GBETA": val}


def spectral_features(n, edges, degrees):
    """Compute V6=IA, V7=IL, V8=IQ by exact dense eigendecomposition."""
    if n <= 0:
        return {"V6": 0.0, "V7": 0.0, "V8": 0.0}
    if args.spectral_max_vertices > 0 and n > args.spectral_max_vertices:
        return {"V6": 0.0, "V7": 0.0, "V8": 0.0}

    try:
        import numpy as np
    except Exception as e:
        raise RuntimeError(f"numpy is required for spectral features: {e}")

    A = np.zeros((n, n), dtype=float)
    for u, v, _weight in edges:
        A[u, v] = 1.0
        A[v, u] = 1.0
    d = np.array(degrees, dtype=float)
    deg_sorted = np.sort(d)[::-1]
    D = np.diag(d)
    L = D - A
    Q = D + A

    def shifted_dot(M):
        eig = np.linalg.eigvalsh(M)[::-1]
        shifted = eig - (np.trace(M) / n)
        return float(np.dot(shifted, deg_sorted))

    return {"V6": shifted_dot(A), "V7": shifted_dot(L), "V8": shifted_dot(Q)}


def is_irregular_from_degrees(degrees):
    if not degrees:
        return False
    return any(d != degrees[0] for d in degrees)


def safe_filename_token(x, max_len=80):
    s = str(x)
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    s = s.strip("._")
    if not s:
        s = "NA"
    return s[:max_len]


def draw_graph_for_kmer(seq_id, seq, label, kmer, raw_nodes, edges, k_len, features, degrees):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import networkx as nx
    except Exception as e:
        print(f"[WARN] Drawing skipped because networkx/matplotlib import failed: {e}", file=sys.stderr)
        return None

    os.makedirs(args.draw_dir, exist_ok=True)

    G = nx.Graph()
    for idx, original_pos in enumerate(raw_nodes):
        G.add_node(original_pos, degree=degrees[idx])
    for u_idx, v_idx, weight in edges:
        G.add_edge(raw_nodes[u_idx], raw_nodes[v_idx], weight=weight)

    plt.figure(figsize=(12, 8))
    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=520, node_color="lightblue", alpha=0.9, edgecolors="black")
    nx.draw_networkx_edges(G, pos, edge_color="green", width=1.5, alpha=0.8)
    nx.draw_networkx_labels(G, pos, labels={node: str(node) for node in G.nodes()}, font_size=10)
    edge_labels = {(u, v): f"d={data['weight']}" for u, v, data in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

    title_seq = safe_filename_token(seq, max_len=30)
    title_parts = [f"{k}={v:.6g}" for k, v in features.items() if k in ("V1", "V2", "V3", "V4", "V5")]
    plt.title(
        f"Irregular graph | row={seq_id} label={label} token='{kmer}' L={k_len}\n"
        f"n={len(raw_nodes)} m={len(edges)} degrees={degrees}\n" + " ".join(title_parts),
        fontsize=12,
        pad=18,
    )
    plt.axis("off")
    plt.tight_layout()

    file_base = os.path.join(
        args.draw_dir,
        f"row{seq_id}_label{safe_filename_token(label, 30)}_seq{title_seq}_kmer{safe_filename_token(kmer, 50)}_L{k_len}",
    )
    png_file = f"{file_base}.png"
    edge_file = f"{file_base}_edges.csv"
    node_file = f"{file_base}_nodes.csv"

    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.close()

    with open(edge_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["source", "target", "weight"])
        for u_idx, v_idx, weight in edges:
            writer.writerow([raw_nodes[u_idx], raw_nodes[v_idx], weight])

    with open(node_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "position", "degree"])
        for idx, original_pos in enumerate(raw_nodes):
            writer.writerow([original_pos, original_pos, degrees[idx]])

    print(f"[DRAW] Irregular graph saved: {png_file}")
    return file_base


# ---------------- Per-sequence computation ----------------
def output_feature_names():
    names = []
    for feat in SELECTED_FEATURES:
        if feat == "BETA":
            names.extend([f"BETA_a{alpha_name(a)}" for a in BETTI_ALPHAS])
        elif feat == "V3":
            if args.irrt_mode == "abs":
                names.append("V3")
            elif args.irrt_mode == "square":
                names.append("V3SQ")
            elif args.irrt_mode == "sign":
                names.append("V3SIGN")
            elif args.irrt_mode == "both":
                names.extend(["V3", "V3SQ"])
            elif args.irrt_mode == "all":
                names.extend(["V3", "V3SQ", "V3SIGN"])
        else:
            names.append(feat)
    return names


OUTPUT_FEATURES = output_feature_names()


def compute_features_for_graph(n, edges):
    feats, deg = degree_features(n, edges)
    if NEED_BETTI:
        feats.update(betti_alpha_features(n, edges, deg, BETTI_ALPHAS))
    if NEED_LGBETA:
        feats.update(line_graph_betti_feature(n, edges, deg))
    if NEED_L2GBETA:
        feats.update(second_line_graph_betti_feature(n, edges, deg))
    if NEED_SPECTRAL:
        feats.update(spectral_features(n, edges, deg))
    return feats, deg


def compute_irregular_for_sequence(seq, label, seq_id):
    seq = str(seq)
    pos_maps = {L: get_positions(seq, L) for L in TOKEN_LENGTHS}
    out = [seq, label]
    drawn_for_this_seq = 0

    for group in TOKEN_GROUPS:
        group_label = token_group_name(group)
        verts, length_by_pos = get_group_positions(pos_maps, group)
        if len(verts) < 1:
            out.extend([0.0] * len(OUTPUT_FEATURES))
            continue

        try:
            raw_nodes, edges = build_simple_graph(verts, length_by_pos=length_by_pos)
            n = len(raw_nodes)
            feats, degrees = compute_features_for_graph(n, edges)

            should_draw = (
                args.draw_irregular
                and n >= args.draw_min_vertices
                and n <= args.draw_max_vertices
                and is_irregular_from_degrees(degrees)
                and abs(feats.get("V1", 0.0)) > args.irregular_eps
                and (args.draw_max_per_seq <= 0 or drawn_for_this_seq < args.draw_max_per_seq)
            )
            if should_draw:
                draw_graph_for_kmer(seq_id, seq, label, group_label, raw_nodes, edges, max(len(x) for x in group), feats, degrees)
                drawn_for_this_seq += 1

        except Exception as e:
            print(f"[ERROR] compute features failed for row={seq_id}, group={group_label}: {e}", file=sys.stderr)
            feats = {name: 0.0 for name in OUTPUT_FEATURES}

        out.extend([feats.get(name, 0.0) for name in OUTPUT_FEATURES])

    return out


def compute_irregular_for_sequence_wrapper(row):
    seq, label, seq_id = row
    return compute_irregular_for_sequence(seq, label, seq_id)


# ---------------- Main ----------------
def read_rows(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "sequence" not in reader.fieldnames or "label" not in reader.fieldnames:
            raise ValueError("Input CSV must contain columns: sequence,label")
        for row_idx, r in enumerate(reader):
            rows.append((r["sequence"], r["label"], row_idx))
    return rows


def main():
    header = ["sequence", "label"]
    for group in TOKEN_GROUPS:
        group_label = token_group_name(group)
        for feat in OUTPUT_FEATURES:
            header.append(f"{group_label}_{feat}")

    print(f"[INFO] Vocab tokens: {len(KMER_LIST)}")
    print(f"[INFO] combo-k: {args.combo_k}; include-smaller={args.combo_include_smaller}; token graphs={len(TOKEN_GROUPS)}")
    print(f"[INFO] Selected features: {','.join(SELECTED_FEATURES)}")
    print(f"[INFO] irrt-mode: {args.irrt_mode}")
    if NEED_BETTI:
        print(f"[INFO] Betti alphas: {','.join(str(a) for a in BETTI_ALPHAS)}")
    if NEED_SPECTRAL:
        print(f"[INFO] Spectral features enabled; spectral-max-vertices={args.spectral_max_vertices}")

    rows = read_rows(input_csv)

    out_dir = os.path.dirname(os.path.abspath(output_csv))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    if args.draw_irregular:
        os.makedirs(args.draw_dir, exist_ok=True)

    if args.draw_irregular and args.workers != 1:
        print(
            "[WARN] --draw-irregular with multiprocessing may create many files in parallel. "
            "For debugging/inspection, consider --workers 1.",
            file=sys.stderr,
        )

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = list(
            tqdm(
                executor.map(compute_irregular_for_sequence_wrapper, rows),
                total=len(rows),
                desc="Processing sequences",
            )
        )

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(results)

    print(f"[OK] {output_csv} generated.")
    if args.draw_irregular:
        print(f"[OK] Irregular graph drawings/data were saved under: {args.draw_dir}")


if __name__ == "__main__":
    main()
