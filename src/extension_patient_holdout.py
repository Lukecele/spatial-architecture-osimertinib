#!/usr/bin/env python3
"""Patient-held-out audit for the public NSCLC extension pilot.

This intentionally uses the same marker specification and outcome proxy as
the existing assay, while fitting only on other patients for each test
patient. It is a validation audit, not a replacement estimator.
"""

import json
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import Delaunay
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from extension_pilot_baselines import GENE_SPEC


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "extension_pilot"
SAMPLES = (
    "P10_B1", "P10_T1", "P11_B1", "P11_T1",
    "P15_B1", "P15_T1", "P16_B1", "P16_T1",
    "P17_B1", "P17_T1", "P19_B1", "P19_T1",
    "P24_B1", "P24_T1", "P25_B1", "P25_T1",
)


def load_sample(sample):
    directory = PILOT / f"pilot_{sample}"
    with h5py.File(directory / f"{sample}-filtered_feature_bc_matrix.h5", "r") as handle:
        matrix = handle["matrix"]
        names = [x.decode() for x in matrix["features"]["name"][:]]
        rows = {name: index for index, name in enumerate(names)}
        data, indices, indptr = matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]
        n = len(matrix["barcodes"])
        genes = sorted(set(sum(GENE_SPEC.values(), [])))
        counts = {gene: np.zeros(n) for gene in genes}
        total = np.zeros(n)
        for column in range(n):
            start, stop = indptr[column], indptr[column + 1]
            observed_rows, values = indices[start:stop], data[start:stop]
            total[column] = values.sum()
            for gene in genes:
                row = rows.get(gene)
                if row is not None:
                    hit = np.flatnonzero(observed_rows == row)
                    if len(hit):
                        counts[gene][column] = values[hit[0]]
        barcodes = [x.decode() for x in matrix["barcodes"][:]]

    positions = {}
    for line in (directory / "tissue_positions_list.csv").read_text().splitlines():
        fields = line.split(",")
        positions[fields[0]] = [int(fields[4]), int(fields[5])]
    coords = np.asarray([positions[x] for x in barcodes], dtype=float)
    normalized = {
        gene: np.log1p(values / np.maximum(total, 1.0) * 1e4)
        for gene, values in counts.items()
    }

    def module(group):
        return np.mean([normalized[gene] for gene in GENE_SPEC[group]], axis=0)

    scores = {group: module(group) for group in GENE_SPEC}
    frequencies = np.vstack([scores["T_S"], scores["T_R"], scores["S_CAF"]])
    frequencies -= frequencies.max(axis=0, keepdims=True)
    frequencies = np.exp(frequencies)
    frequencies /= frequencies.sum(axis=0, keepdims=True)
    environment = np.column_stack([scores["HYPOXIA"], scores["VASCULAR"]])
    residual = scores["FITNESS"] - LinearRegression().fit(
        environment, scores["FITNESS"]
    ).predict(environment)
    return coords, frequencies.T, residual


def neighborhood(coords, frequencies):
    graph = [set() for _ in coords]
    for simplex in Delaunay(coords).simplices:
        for node in simplex:
            graph[node].update(x for x in simplex if x != node)
    return np.asarray([
        frequencies[list(neighbors)].mean(axis=0) if neighbors else frequencies[index]
        for index, neighbors in enumerate(graph)
    ])


def design(coords, frequencies):
    local = neighborhood(coords, frequencies)
    game = np.column_stack([
        frequencies[:, row, None] * local for row in range(frequencies.shape[1])
    ])
    composition = np.column_stack([frequencies, local])
    return composition, game


def patient(sample):
    return sample.split("_")[0]


def main():
    data = {sample: load_sample(sample) for sample in SAMPLES}
    rows = []
    for held_out in sorted({patient(sample) for sample in SAMPLES}):
        train_samples = [sample for sample in SAMPLES if patient(sample) != held_out]
        test_samples = [sample for sample in SAMPLES if patient(sample) == held_out]
        train_composition = np.vstack(
            [design(*data[sample][:2])[0] for sample in train_samples]
        )
        train_game = np.vstack(
            [design(*data[sample][:2])[1] for sample in train_samples]
        )
        train_y = np.concatenate([data[sample][2] for sample in train_samples])
        composition_model = LinearRegression(fit_intercept=False).fit(
            train_composition, train_y
        )
        game_model = LinearRegression(fit_intercept=False).fit(train_game, train_y)
        for sample in test_samples:
            coords, frequencies, outcome = data[sample]
            test_composition, test_game = design(coords, frequencies)
            composition_prediction = composition_model.predict(test_composition)
            game_prediction = game_model.predict(test_game)
            rows.append({
                "held_out_patient": held_out,
                "sample": sample,
                "spots": int(len(outcome)),
                "composition_r2": float(r2_score(outcome, composition_prediction)),
                "game_r2": float(r2_score(outcome, game_prediction)),
                "game_minus_composition_r2": float(
                    r2_score(outcome, game_prediction)
                    - r2_score(outcome, composition_prediction)
                ),
                "game_rmse": float(
                    np.sqrt(np.mean((outcome - game_prediction) ** 2))
                ),
            })
    output = {
        "samples": list(SAMPLES),
        "patients": sorted({patient(sample) for sample in SAMPLES}),
        "rows": rows,
        "interpretation": "patient-held-out exploratory validation; not clinical validation",
    }
    path = PILOT / "patient_holdout_results.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(path)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
