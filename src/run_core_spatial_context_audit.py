#!/usr/bin/env python3
"""Run the minimal non-game spatial-context audit on the bundled cache."""

import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "clinical_pan_cancer_cache.npz"
OUTPUT = ROOT / "results" / "CORE_SPATIAL_CONTEXT_AUDIT.json"
RADIUS_UM = 150.0

SPECS = {
    "breast": {
        "states": (("EPCAM", "KRT8", "CDH1"), ("VIM", "CDH2", "AXL")),
        "outcome": ("MKI67", "PCNA", "TOP2A", "CDK1"),
        "environment": (("HIF1A", "VEGFA", "CA9"), ("PECAM1", "VWF")),
    },
    "melanoma": {
        "states": (("MITF", "MLANA", "PMEL"), ("VIM", "AXL")),
        "outcome": ("MKI67", "PCNA", "TOP2A", "CDK1"),
        "environment": (("HIF1A", "VEGFA", "CA9"), ("PECAM1", "VWF")),
    },
    "lung": {
        "states": (("EGFR", "KRT5", "TP63"), ("VIM", "CDH2", "AXL")),
        "outcome": ("MKI67", "PCNA", "TOP2A", "CDK1"),
        "environment": (("HIF1A", "VEGFA", "CA9"), ("PECAM1", "VWF")),
    },
}


def module(cache, cohort, genes):
    values = [cache[f"{cohort}_gene_{gene}"] for gene in genes]
    return np.mean(values, axis=0)


def prepare(cache, cohort):
    spec = SPECS[cohort]
    state_scores = np.vstack(
        [module(cache, cohort, genes) for genes in spec["states"]]
    )
    state_scores -= state_scores.max(axis=0, keepdims=True)
    frequencies = np.exp(state_scores)
    frequencies /= frequencies.sum(axis=0, keepdims=True)

    outcome = module(cache, cohort, spec["outcome"])
    environment = np.column_stack(
        [module(cache, cohort, genes) for genes in spec["environment"]]
    )
    outcome = outcome - LinearRegression().fit(environment, outcome).predict(environment)
    return frequencies.T, outcome


def neighborhood(coords, frequencies):
    coords_um = coords * (100.0 / np.median(cKDTree(coords).query(coords, k=2)[0][:, 1]))
    indices = cKDTree(coords_um).query_ball_point(coords_um, RADIUS_UM)
    local = np.empty_like(frequencies)
    for index, neighbors in enumerate(indices):
        neighbors = [neighbor for neighbor in neighbors if neighbor != index]
        local[index] = frequencies[neighbors].mean(axis=0) if neighbors else frequencies[index]
    return local


def evaluate(design, outcome, coords):
    in_sample = LinearRegression(fit_intercept=False).fit(design, outcome)
    folds = KMeans(n_clusters=5, random_state=42, n_init=10).fit_predict(coords)
    prediction = np.empty(len(outcome))
    for fold in range(5):
        train, test = folds != fold, folds == fold
        model = LinearRegression(fit_intercept=False).fit(design[train], outcome[train])
        prediction[test] = model.predict(design[test])
    return {
        "in_sample_r2": float(in_sample.score(design, outcome)),
        "spatial_block_oof_r2": float(r2_score(outcome, prediction)),
    }


def oof_score(design, outcome, folds):
    prediction = np.empty(len(outcome))
    for fold in range(5):
        train, test = folds != fold, folds == fold
        model = LinearRegression(fit_intercept=False).fit(design[train], outcome[train])
        prediction[test] = model.predict(design[test])
    return float(r2_score(outcome, prediction))


def oof_increment(focal, context, outcome, folds):
    focal_score = oof_score(focal, outcome, folds)
    context_score = oof_score(context, outcome, folds)
    return context_score - focal_score


def permutation_null(focal, local, outcome, folds, permutations=100):
    observed = oof_increment(
        focal,
        np.column_stack([focal, local]),
        outcome,
        folds,
    )
    rng = np.random.default_rng(42)
    null = np.asarray([
        oof_increment(
            focal,
            np.column_stack([focal, local[rng.permutation(len(local))]]),
            outcome,
            folds,
        )
        for _ in range(permutations)
    ])
    return {
        "observed_oof_increment": float(observed),
        "null_mean_oof_increment": float(null.mean()),
        "null_max_oof_increment": float(null.max()),
        "empirical_p_value": float(
            (1 + np.sum(null >= observed)) / (permutations + 1)
        ),
        "permutations": permutations,
    }


def main():
    results = {}
    with np.load(CACHE) as cache:
        for cohort in SPECS:
            frequencies, outcome = prepare(cache, cohort)
            local = neighborhood(cache[f"{cohort}_coords"], frequencies)
            focal = frequencies
            context = np.column_stack([frequencies, local])
            folds = KMeans(n_clusters=5, random_state=42, n_init=10).fit_predict(
                cache[f"{cohort}_coords"]
            )
            results[cohort] = {
                "spots": int(len(outcome)),
                "radius_um": RADIUS_UM,
                "focal_composition": evaluate(
                    focal, outcome, cache[f"{cohort}_coords"]
                ),
                "focal_plus_neighbourhood": evaluate(
                    context, outcome, cache[f"{cohort}_coords"]
                ),
                "neighbourhood_permutation_null": permutation_null(
                    focal,
                    local,
                    outcome,
                    folds,
                ),
                "interpretation": (
                    "Non-game spatial-context audit; not a payoff or causal model."
                ),
            }
    output = {
        "purpose": "Minimal validated-core audit",
        "radius_um": RADIUS_UM,
        "cohorts": results,
        "limitations": [
            "The cache contains restricted marker panels, not whole-transcriptome data.",
            "Spatial blocks are not independent patients.",
            "The outcome is a molecular proliferation proxy.",
            "No result establishes causality or biological payoff parameters.",
        ],
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
