#!/usr/bin/env python3
"""Summarize the existing multiscale scan without adding a new claim.

The scan evaluates several radii on the same tissue. This script applies
Benjamini-Hochberg correction to the in-sample exploratory p-values and reports
the held-out spatial-block effect separately. It does not treat the best
radius as independently validated.
"""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "MULTISCALE_SPATIAL_SCAN_REPORT.json"
OUTPUT = ROOT / "results" / "MULTISCALE_AUDIT_SUMMARY.json"


def benjamini_hochberg(p_values):
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values)
    running = 1.0
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * len(p_values) / rank)
        adjusted[index] = running
    return adjusted


def summarize_cohort(cohort, data):
    rows = data["multiscale_scan"]
    adjusted = benjamini_hochberg([row["p_value"] for row in rows])
    for row, q_value in zip(rows, adjusted):
        row["q_value_bh"] = float(q_value)

    best = max(rows, key=lambda row: row["delta_oof_pct"])
    positive_oof = [row for row in rows if row["delta_oof_pct"] > 0]
    return {
        "cohort": cohort,
        "n_scales": len(rows),
        "correlation_length_xi_um": data["correlation_length_xi_um"],
        "best_oof_radius_um": best["radius_um"],
        "best_delta_oof_pct": best["delta_oof_pct"],
        "contact_scale_delta_oof_pct": rows[0]["delta_oof_pct"],
        "positive_oof_scales": len(positive_oof),
        "all_in_sample_q_below_0_05": bool(all(row["q_value_bh"] < 0.05 for row in rows)),
        "scale_results": [
            {
                "radius_um": row["radius_um"],
                "delta_oof_pct": row["delta_oof_pct"],
                "delta_r2_pct": row["delta_r2_pct"],
                "p_value": row["p_value"],
                "q_value_bh": row["q_value_bh"],
            }
            for row in rows
        ],
        "interpretation": (
            "Exploratory same-section multiscale scan; the selected radius is "
            "not independently validated and does not establish causality."
        ),
    }


def main():
    source = json.loads(INPUT.read_text())
    output = {
        "purpose": "Multiple-scale spatial audit with separate held-out reporting",
        "cohorts": {
            cohort: summarize_cohort(cohort, data)
            for cohort, data in source.items()
        },
        "limitations": [
            "Radii were evaluated on the same tissue, so the best radius is selected, not independent.",
            "Benjamini-Hochberg q-values apply only to exploratory in-sample p-values.",
            "Spatial-block OOF performance, not in-sample p-values, is the relevant generalization diagnostic.",
            "The scan does not provide patient-level replication or causal identification.",
        ],
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
