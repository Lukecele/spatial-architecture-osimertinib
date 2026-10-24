#!/usr/bin/env python3
"""Patient-held-out test of the pre-specified 150--300 um interaction window."""

import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from extension_patient_holdout import SAMPLES, load_sample


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "extension_pilot"
RADII_UM = (150.0, 200.0, 250.0, 300.0)


def patient(sample):
    return sample.split("_")[0]


def scale_coordinates(coords):
    distances, _ = cKDTree(coords).query(coords, k=2)
    median_nearest_neighbor = float(np.median(distances[:, 1]))
    return coords * (100.0 / median_nearest_neighbor), median_nearest_neighbor


def local_frequencies(coords_um, frequencies, radius_um):
    tree = cKDTree(coords_um)
    neighborhoods = tree.query_ball_point(coords_um, radius_um)
    local = np.empty_like(frequencies)
    for index, neighbors in enumerate(neighborhoods):
        neighbors = [neighbor for neighbor in neighbors if neighbor != index]
        local[index] = (
            frequencies[neighbors].mean(axis=0)
            if neighbors
            else frequencies[index]
        )
    return local


def designs(coords, frequencies, radius_um):
    coords_um, median_nearest_neighbor = scale_coordinates(coords)
    local = local_frequencies(coords_um, frequencies, radius_um)
    composition = np.column_stack([frequencies, local])
    game = np.column_stack(
        [frequencies[:, row, None] * local for row in range(frequencies.shape[1])]
    )
    return composition, game, median_nearest_neighbor


def bootstrap_mean(values, seed=42, replicates=10000):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    draws = rng.choice(values, size=(replicates, len(values)), replace=True)
    means = draws.mean(axis=1)
    return [
        float(np.quantile(means, 0.025)),
        float(np.quantile(means, 0.975)),
    ]


def main():
    data = {sample: load_sample(sample) for sample in SAMPLES}
    patients = sorted({patient(sample) for sample in SAMPLES})
    rows = []
    summaries = {}

    for radius_um in RADII_UM:
        prepared = {}
        for sample, (coords, frequencies, outcome) in data.items():
            composition, game, median_nn = designs(coords, frequencies, radius_um)
            prepared[sample] = (composition, game, outcome, median_nn)

        for held_out in patients:
            train_samples = [
                sample for sample in SAMPLES if patient(sample) != held_out
            ]
            test_samples = [
                sample for sample in SAMPLES if patient(sample) == held_out
            ]
            train_composition = np.vstack(
                [prepared[sample][0] for sample in train_samples]
            )
            train_game = np.vstack([prepared[sample][1] for sample in train_samples])
            train_y = np.concatenate([prepared[sample][2] for sample in train_samples])
            composition_model = LinearRegression(fit_intercept=False).fit(
                train_composition, train_y
            )
            game_model = LinearRegression(fit_intercept=False).fit(
                train_game, train_y
            )

            patient_rows = []
            for sample in test_samples:
                composition, game, outcome, median_nn = prepared[sample]
                composition_r2 = r2_score(
                    outcome, composition_model.predict(composition)
                )
                game_r2 = r2_score(outcome, game_model.predict(game))
                row = {
                    "radius_um": radius_um,
                    "held_out_patient": held_out,
                    "sample": sample,
                    "spots": int(len(outcome)),
                    "median_nearest_neighbor_px": median_nn,
                    "composition_r2": float(composition_r2),
                    "game_r2": float(game_r2),
                    "game_minus_composition_r2": float(game_r2 - composition_r2),
                }
                rows.append(row)
                patient_rows.append(row["game_minus_composition_r2"])

            patient_delta = float(np.mean(patient_rows))
            summaries.setdefault(str(int(radius_um)), []).append(patient_delta)

    summary = {}
    for radius, values in summaries.items():
        values = np.asarray(values)
        summary[radius] = {
            "n_patients": len(values),
            "mean_patient_delta_r2": float(values.mean()),
            "median_patient_delta_r2": float(np.median(values)),
            "patients_with_positive_delta": int(np.sum(values > 0)),
            "bootstrap_95ci_mean_delta_r2": bootstrap_mean(values),
        }

    output = {
        "radii_um": list(RADII_UM),
        "patients": patients,
        "samples": list(SAMPLES),
        "rows": rows,
        "summary_by_radius_um": summary,
        "pre_specification": (
            "150--300 um window selected before this patient-held-out analysis "
            "from the prior same-section multiscale scan"
        ),
        "interpretation": (
            "Patient-held-out exploratory validation; no radius is selected "
            "after inspecting these results."
        ),
    }
    destination = PILOT / "multiscale_patient_holdout_results.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(destination)
    for radius, values in summary.items():
        print(radius, values)


if __name__ == "__main__":
    main()
