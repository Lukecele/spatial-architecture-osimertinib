#!/usr/bin/env python3
"""Apply only mathematically identifiable diagnostics to the bundled cache.

This module deliberately does not infer payoffs, causality, or temporal
dynamics.  It evaluates spatial graph geometry and records which proposed
mathematical tools are identifiable from a single spatial snapshot.
"""

import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csgraph
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "clinical_pan_cancer_cache.npz"
OUTPUT = ROOT / "results" / "MATHEMATICAL_APPLICABILITY_AUDIT.json"

MARKERS = {
    "breast": (("EPCAM", "KRT8", "KRT18"), ("VIM", "CDH2", "AXL")),
    "melanoma": (("MITF", "MLANA", "PMEL"), ("VIM", "AXL")),
    "lung": (("EGFR", "KRT5", "TP63"), ("VIM", "CDH2", "AXL")),
}


def knn_graph(coords, k=6):
    """Return an unweighted symmetric k-nearest-neighbour graph."""
    distances, indices = cKDTree(coords).query(coords, k=k + 1)
    rows = np.repeat(np.arange(len(coords)), k)
    cols = indices[:, 1:].reshape(-1)
    data = np.ones(len(rows), dtype=float)
    directed = coo_matrix((data, (rows, cols)), shape=(len(coords), len(coords)))
    graph = directed.maximum(directed.T).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    return graph


def phenotype_signal(cache, cohort):
    positive, resistant = MARKERS[cohort]

    def module(genes):
        values = [
            cache[f"{cohort}_gene_{gene}"]
            for gene in genes
            if f"{cohort}_gene_{gene}" in cache
        ]
        if not values:
            raise KeyError(f"No marker found for {cohort}: {genes}")
        return np.mean(values, axis=0)

    scores = np.vstack((module(positive), module(resistant)))
    scores -= scores.max(axis=0, keepdims=True)
    frequencies = np.exp(scores)
    frequencies /= frequencies.sum(axis=0, keepdims=True)
    return frequencies[0]


def graph_diagnostics(coords, signal, k=6, n_eigenvalues=24):
    graph = knn_graph(coords, k=k)
    laplacian = csgraph.laplacian(graph, normed=True)
    n_eigenvalues = min(n_eigenvalues, len(coords) - 2)
    eigenvalues = eigsh(
        laplacian,
        k=n_eigenvalues + 1,
        which="SM",
        return_eigenvectors=False,
        tol=1e-7,
    )
    eigenvalues = np.sort(np.maximum(eigenvalues, 0.0))
    positive = eigenvalues[eigenvalues > 1e-8]

    centered = signal - signal.mean()
    degree = np.asarray(graph.sum(axis=1)).ravel()
    dirichlet = float(centered @ (graph @ centered) / max(graph.nnz, 1))
    variance = float(np.var(centered))
    moran_like = float(
        len(signal) * (centered @ (graph @ centered))
        / max(graph.nnz * (centered @ centered), 1e-15)
    )

    # This is a truncated low-frequency spectral zeta, not the full zeta
    # function of the graph. It is a geometry descriptor only.
    spectral_zeta = {
        str(s): float(np.sum(positive ** (-s))) for s in (1, 2)
    }
    normalized = positive / max(float(positive.sum()), 1e-15)
    spectral_entropy = float(-np.sum(normalized * np.log(normalized)))

    return {
        "n_nodes": int(len(signal)),
        "n_edges": int(graph.nnz // 2),
        "n_connected_components": int(csgraph.connected_components(graph, directed=False)[0]),
        "k": int(k),
        "mean_degree": float(degree.mean()),
        "low_frequency_eigenvalues": eigenvalues.tolist(),
        "truncated_spectral_zeta": spectral_zeta,
        "truncated_spectral_entropy": spectral_entropy,
        "phenotype_dirichlet_energy": dirichlet,
        "phenotype_variance": variance,
        "phenotype_moran_like": moran_like,
    }


def applicability_registry():
    return {
        "identifiable_now": {
            "pair_correlation_and_Ripley_K": "spatial clustering and interaction scale",
            "graph_Laplacian": "low-frequency tissue geometry and smoothness",
            "truncated_spectral_zeta": "descriptive summary of the graph spectrum",
            "conditional_information": "only with careful permutation and bias control",
            "hierarchical_bootstrap": "patient-level uncertainty when patients exist",
        },
        "requires_new_data": {
            "Z_transform": "requires repeated measurements or a discrete time series",
            "Koopman_or_SDE": "requires observed transitions over time",
            "Pontryagin_control": "requires intervention response and a calibrated state model",
            "genetic_selection": "requires mutation, CNV, lineage, or matched single-cell data",
        },
        "not_a_valid_current_claim": {
            "Riemann_zeta": "no biological observable or null hypothesis is supplied by the snapshot",
            "quantum_or_relativistic_models": "no measured phenomenon would distinguish them here",
            "Millennium_prize_analogies": "mathematical prestige does not identify biological parameters",
        },
    }


def main():
    with np.load(CACHE) as cache:
        cohorts = {}
        for cohort in ("breast", "melanoma", "lung"):
            signal = phenotype_signal(cache, cohort)
            cohorts[cohort] = graph_diagnostics(cache[f"{cohort}_coords"], signal)

    output = {
        "purpose": "Identifiable mathematical diagnostics; not causal or temporal validation",
        "methods": applicability_registry(),
        "cohorts": cohorts,
        "interpretation": [
            "The truncated spectral zeta describes graph geometry, not a payoff matrix.",
            "No result here supports a Riemann-zeta biological mechanism.",
            "A spectral or spatial association must still beat composition and spatial baselines in held-out patients.",
        ],
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
