#!/usr/bin/env python3
"""Summarize patient-held-out performance with patient-level bootstrap."""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "extension_pilot" / "patient_holdout_results.json"
OUTPUT = ROOT / "data" / "extension_pilot" / "patient_holdout_summary.json"


def main():
    report = json.loads(INPUT.read_text())
    by_patient = {}
    for row in report["rows"]:
        by_patient.setdefault(row["held_out_patient"], []).append(row)
    patient_rows = []
    for patient, rows in sorted(by_patient.items()):
        deltas = np.asarray([row["game_minus_composition_r2"] for row in rows])
        patient_rows.append({
            "patient": patient,
            "n_samples": len(rows),
            "mean_delta_r2": float(deltas.mean()),
            "mean_game_r2": float(np.mean([row["game_r2"] for row in rows])),
            "mean_composition_r2": float(
                np.mean([row["composition_r2"] for row in rows])
            ),
        })

    deltas = np.asarray([row["mean_delta_r2"] for row in patient_rows])
    rng = np.random.default_rng(42)
    boot = np.asarray([
        rng.choice(deltas, size=len(deltas), replace=True).mean()
        for _ in range(10000)
    ])
    summary = {
        "patients": patient_rows,
        "n_patients": int(len(deltas)),
        "mean_patient_delta_r2": float(deltas.mean()),
        "median_patient_delta_r2": float(np.median(deltas)),
        "patients_with_positive_delta": int(np.sum(deltas > 0)),
        "bootstrap_95ci_mean_delta_r2": [
            float(np.quantile(boot, 0.025)),
            float(np.quantile(boot, 0.975)),
        ],
        "bootstrap_replicates": 10000,
        "seed": 42,
        "interpretation": (
            "patient-level exploratory summary; not clinical validation and "
            "not evidence of causal payoff parameters"
        ),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n")
    print(OUTPUT)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
