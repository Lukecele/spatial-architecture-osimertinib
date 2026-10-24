#!/usr/bin/env python3
"""Dependency-free hierarchical audit for the NSCLC extension.

This is a fixed-effect/clustered audit, not a full random-effects model. It
keeps spots nested within samples and patients by evaluating patient-level
held-out predictions and reports patient-cluster summaries.
"""

import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge

from extension_patient_holdout import SAMPLES, load_sample, patient
from extension_pilot_baselines import neighborhood_frequencies


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "extension_pilot"


def design(coords, frequencies, median_umi, median_features, game):
    local = neighborhood_frequencies(coords, frequencies)
    composition = np.column_stack([frequencies, local])
    interaction = np.column_stack([
        frequencies[:, row, None] * local
        for row in range(frequencies.shape[1])
    ])
    selected = interaction if game else composition
    log_qc = np.column_stack([
        np.ones(len(coords)),
        np.full(len(coords), np.log1p(median_umi)),
        np.full(len(coords), np.log1p(median_features)),
    ])
    return np.column_stack([selected, log_qc])


def qc_lookup():
    import csv

    with (PILOT / "pilot_qc.tsv").open() as handle:
        return {
            row["sample_id"]: (
                float(row["median_umi"]),
                float(row["median_detected_features"]),
            )
            for row in csv.DictReader(handle, delimiter="\t")
        }


def main():
    qc = qc_lookup()
    data = {sample: load_sample(sample) for sample in SAMPLES}
    patients = sorted({patient(sample) for sample in SAMPLES})
    rows = []
    for held_out in patients:
        train = [sample for sample in SAMPLES if patient(sample) != held_out]
        test = [sample for sample in SAMPLES if patient(sample) == held_out]
        models = {}
        for game in (False, True):
            x_train = np.vstack([
                design(*data[sample][:2], *qc[sample], game)
                for sample in train
            ])
            y_train = np.concatenate([data[sample][2] for sample in train])
            models[game] = Ridge(alpha=1e-3, fit_intercept=False).fit(
                x_train, y_train
            )
        for sample in test:
            coords, frequencies, outcome = data[sample]
            predictions = {
                game: models[game].predict(
                    design(*data[sample][:2], *qc[sample], game)
                )
                for game in (False, True)
            }
            variance = np.sum((outcome - outcome.mean()) ** 2)
            r2 = {
                game: float(
                    1 - np.sum((outcome - predictions[game]) ** 2) / variance
                )
                for game in (False, True)
            }
            rows.append({
                "patient": held_out,
                "sample": sample,
                "baseline_r2": r2[False],
                "game_r2": r2[True],
                "game_minus_baseline_r2": r2[True] - r2[False],
            })

    patient_summary = []
    for held_out in patients:
        values = [
            row["game_minus_baseline_r2"]
            for row in rows if row["patient"] == held_out
        ]
        patient_summary.append({
            "patient": held_out,
            "mean_delta_r2": float(np.mean(values)),
        })
    deltas = np.asarray([row["mean_delta_r2"] for row in patient_summary])
    rng = np.random.default_rng(42)
    bootstrap = np.asarray([
        rng.choice(deltas, len(deltas), replace=True).mean()
        for _ in range(10000)
    ])
    output = {
        "rows": rows,
        "patient_summary": patient_summary,
        "mean_patient_delta_r2": float(deltas.mean()),
        "bootstrap_95ci_mean_delta_r2": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "covariates": ["log_median_umi", "log_median_detected_features"],
        "method": (
            "patient-held-out ridge audit with QC covariates; not a full "
            "random-effects mixed model"
        ),
    }
    destination = PILOT / "hierarchical_qc_audit.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(destination)
    print("mean delta:", output["mean_patient_delta_r2"])
    print("95% CI:", output["bootstrap_95ci_mean_delta_r2"])


if __name__ == "__main__":
    main()
