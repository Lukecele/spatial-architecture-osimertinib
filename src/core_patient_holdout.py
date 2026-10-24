#!/usr/bin/env python3
"""Patient-held-out validation of the minimal non-game spatial-context model."""

import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from extension_patient_holdout import SAMPLES, load_sample


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "extension_pilot" / "core_patient_holdout.json"
RADIUS_UM = 150.0


def patient(sample):
    return sample.split("_")[0]


def local_context(coords, frequencies):
    tree = cKDTree(coords * (100.0 / np.median(tree_distances(coords))))
    neighborhoods = tree.query_ball_point(tree.data, RADIUS_UM)
    local = np.empty_like(frequencies)
    for index, neighbors in enumerate(neighborhoods):
        neighbors = [neighbor for neighbor in neighbors if neighbor != index]
        local[index] = frequencies[neighbors].mean(axis=0) if neighbors else frequencies[index]
    return local


def tree_distances(coords):
    return cKDTree(coords).query(coords, k=2)[0][:, 1]


def bootstrap(values, seed=42, replicates=10000):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    means = rng.choice(values, size=(replicates, len(values)), replace=True).mean(axis=1)
    return [
        float(np.quantile(means, 0.025)),
        float(np.quantile(means, 0.975)),
    ]


def main():
    data = {sample: load_sample(sample) for sample in SAMPLES}
    patients = sorted({patient(sample) for sample in SAMPLES})
    rows = []

    for held_out in patients:
        train_samples = [sample for sample in SAMPLES if patient(sample) != held_out]
        test_samples = [sample for sample in SAMPLES if patient(sample) == held_out]
        train_focal = []
        train_context = []
        train_outcome = []

        prepared = {}
        for sample, (coords, frequencies, outcome) in data.items():
            local = local_context(coords, frequencies)
            focal = frequencies
            context = np.column_stack([frequencies, local])
            prepared[sample] = (focal, context, outcome)
            if sample in train_samples:
                train_focal.append(focal)
                train_context.append(context)
                train_outcome.append(outcome)

        focal_model = LinearRegression(fit_intercept=False).fit(
            np.vstack(train_focal), np.concatenate(train_outcome)
        )
        context_model = LinearRegression(fit_intercept=False).fit(
            np.vstack(train_context), np.concatenate(train_outcome)
        )

        patient_deltas = []
        for sample in test_samples:
            focal, context, outcome = prepared[sample]
            focal_r2 = r2_score(outcome, focal_model.predict(focal))
            context_r2 = r2_score(outcome, context_model.predict(context))
            delta = context_r2 - focal_r2
            patient_deltas.append(delta)
            rows.append({
                "held_out_patient": held_out,
                "sample": sample,
                "spots": int(len(outcome)),
                "focal_r2": float(focal_r2),
                "focal_plus_neighbourhood_r2": float(context_r2),
                "context_minus_focal_r2": float(delta),
            })

    patient_summary = []
    for held_out in patients:
        values = [
            row["context_minus_focal_r2"]
            for row in rows
            if row["held_out_patient"] == held_out
        ]
        patient_summary.append({
            "patient": held_out,
            "mean_delta_r2": float(np.mean(values)),
            "samples": len(values),
        })
    deltas = np.asarray([row["mean_delta_r2"] for row in patient_summary])

    output = {
        "radius_um": RADIUS_UM,
        "patients": patients,
        "rows": rows,
        "patient_summary": patient_summary,
        "mean_patient_delta_r2": float(deltas.mean()),
        "median_patient_delta_r2": float(np.median(deltas)),
        "patients_with_positive_delta": int(np.sum(deltas > 0)),
        "bootstrap_95ci_mean_delta_r2": bootstrap(deltas),
        "interpretation": (
            "Patient-held-out validation of non-game spatial context; "
            "molecular proxy only, not causal or clinical validation."
        ),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(OUTPUT)
    print(json.dumps({
        "mean_patient_delta_r2": output["mean_patient_delta_r2"],
        "patients_with_positive_delta": output["patients_with_positive_delta"],
        "bootstrap_95ci_mean_delta_r2": output["bootstrap_95ci_mean_delta_r2"],
    }, indent=2))


if __name__ == "__main__":
    main()
