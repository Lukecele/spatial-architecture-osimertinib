#!/usr/bin/env python3
"""Compare composition-only and game baselines on the paired public pilot.

This is an exploratory audit. It does not pool patients or estimate external
validity. The outcome is the same environmental-residual fitness proxy used by
the existing assay.
"""

import json
from pathlib import Path

import h5py
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
from sklearn.metrics import r2_score


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "extension_pilot"
GENE_SPEC = {
    "T_S": ["EGFR", "KRT5", "TP63", "EPCAM"],
    "T_R": ["VIM", "AXL", "CDH2", "MMP2"],
    "S_CAF": ["ACTA2", "COL1A1", "FAP"],
    "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
    "VASCULAR": ["PECAM1", "VWF"],
}


def load_sample(sample: str, gene_spec=GENE_SPEC):
    sample_dir = PILOT / f"pilot_{sample}"
    matrix_path = sample_dir / f"{sample}-filtered_feature_bc_matrix.h5"
    with h5py.File(matrix_path, "r") as handle:
        matrix = handle["matrix"]
        names = [value.decode() for value in matrix["features"]["name"][:]]
        name_to_row = {name: index for index, name in enumerate(names)}
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        barcodes = [value.decode() for value in matrix["barcodes"][:]]

        all_genes = sorted(set(sum(gene_spec.values(), [])))
        counts = {gene: np.zeros(len(barcodes), dtype=float) for gene in all_genes}
        total_umi = np.zeros(len(barcodes), dtype=float)
        for column in range(len(barcodes)):
            start, stop = indptr[column], indptr[column + 1]
            rows = indices[start:stop]
            values = data[start:stop]
            total_umi[column] = values.sum()
            for gene in all_genes:
                row = name_to_row.get(gene)
                if row is not None:
                    hits = np.flatnonzero(rows == row)
                    if len(hits):
                        counts[gene][column] = values[hits[0]]

    coordinates = {}
    for line in (sample_dir / "tissue_positions_list.csv").read_text().splitlines():
        fields = line.split(",")
        if fields:
            coordinates[fields[0]] = [int(fields[4]), int(fields[5])]
    coords = np.asarray([coordinates[barcode] for barcode in barcodes], dtype=float)
    normalized = {
        gene: np.log1p((values / np.maximum(total_umi, 1.0)) * 1e4)
        for gene, values in counts.items()
    }

    def module(genes):
        return np.mean([normalized[gene] for gene in genes], axis=0)

    scores = {name: module(genes) for name, genes in gene_spec.items()}
    expression = np.vstack([scores["T_S"], scores["T_R"], scores["S_CAF"]])
    expression -= expression.max(axis=0, keepdims=True)
    frequencies = np.exp(expression)
    frequencies /= frequencies.sum(axis=0, keepdims=True)
    environment = np.column_stack([scores["HYPOXIA"], scores["VASCULAR"]])
    residual = scores["FITNESS"] - LinearRegression().fit(
        environment, scores["FITNESS"]
    ).predict(environment)
    return coords, frequencies.T, residual


def neighborhood_frequencies(coords, frequencies):
    from scipy.spatial import Delaunay

    neighbors = [set() for _ in range(len(coords))]
    for simplex in Delaunay(coords).simplices:
        for left in simplex:
            neighbors[left].update(index for index in simplex if index != left)
    local = np.empty_like(frequencies)
    for index, adjacent in enumerate(neighbors):
        local[index] = frequencies[list(adjacent)].mean(axis=0) if adjacent else frequencies[index]
    return local


def fit_r2(design, outcome):
    return float(LinearRegression(fit_intercept=False).fit(design, outcome).score(design, outcome))


def block_cv(design, outcome, coords, n_blocks=5):
    folds = KMeans(n_clusters=n_blocks, random_state=42, n_init=10).fit_predict(coords)
    predictions = np.full(len(outcome), np.nan)
    for fold in range(n_blocks):
        train = folds != fold
        test = folds == fold
        model = LinearRegression(fit_intercept=False).fit(design[train], outcome[train])
        predictions[test] = model.predict(design[test])
    return float(r2_score(outcome, predictions)), int(n_blocks)


def shuffled_neighbor_null(frequencies, local, outcome, n_permutations=100):
    observed = fit_r2(
        np.column_stack(
            [frequencies[:, row, None] * local for row in range(frequencies.shape[1])]
        ),
        outcome,
    )
    rng = np.random.default_rng(42)
    null = []
    for _ in range(n_permutations):
        shuffled = local[rng.permutation(len(local))]
        design = np.column_stack(
            [
                frequencies[:, row, None] * shuffled
                for row in range(frequencies.shape[1])
            ]
        )
        null.append(fit_r2(design, outcome))
    return {
        "observed_r2": observed,
        "null_mean_r2": float(np.mean(null)),
        "null_max_r2": float(np.max(null)),
        "empirical_p_value": float((1 + np.sum(np.asarray(null) >= observed)) / (n_permutations + 1)),
        "permutations": n_permutations,
    }


def main():
    output = {}
    for sample in ("P10_B1", "P10_T1"):
        coords, frequencies, outcome = load_sample(sample)
        local = neighborhood_frequencies(coords, frequencies)
        composition = np.column_stack([frequencies, local])
        game = np.column_stack(
            [frequencies[:, row, None] * local for row in range(frequencies.shape[1])]
        )
        output[sample] = {
            "spots": int(len(outcome)),
            "mean_only_r2": fit_r2(np.ones((len(outcome), 1)), outcome),
            "composition_and_neighborhood_r2": fit_r2(composition, outcome),
            "game_r2": fit_r2(game, outcome),
            "composition_and_neighborhood_block_cv_r2": block_cv(
                composition, outcome, coords
            )[0],
            "game_block_cv_r2": block_cv(game, outcome, coords)[0],
            "game_neighbor_shuffle_null": shuffled_neighbor_null(
                frequencies, local, outcome
            ),
            "interpretation": "exploratory in-sample comparison; not external validation",
        }
    destination = PILOT / "pilot_baseline_comparison.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(destination)
    for sample, values in output.items():
        print(sample, values)


if __name__ == "__main__":
    main()
