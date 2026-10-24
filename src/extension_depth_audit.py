#!/usr/bin/env python3
"""Audit whether held-out game performance tracks sequencing depth."""

import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "extension_pilot"


def main():
    holdout = json.loads((PILOT / "patient_holdout_results.json").read_text())
    qc = {}
    with (PILOT / "pilot_qc.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            qc[row["sample_id"]] = row

    rows = []
    for result in holdout["rows"]:
        sample = result["sample"]
        if sample not in qc:
            raise KeyError(f"Missing QC row for {sample}")
        row = {
            **result,
            "median_umi": float(qc[sample]["median_umi"]),
            "median_detected_features": float(
                qc[sample]["median_detected_features"]
            ),
        }
        rows.append(row)

    def corr(key):
        x = np.asarray([row[key] for row in rows], dtype=float)
        y = np.asarray([row["game_minus_composition_r2"] for row in rows])
        return float(np.corrcoef(x, y)[0, 1])

    output = {
        "rows": rows,
        "correlation_depth_delta_r2": corr("median_umi"),
        "correlation_detected_features_delta_r2": corr(
            "median_detected_features"
        ),
        "interpretation": (
            "diagnostic association only; sequencing depth is not randomized "
            "and this does not establish or remove biological confounding"
        ),
    }
    destination = PILOT / "depth_confounding_audit.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(destination)
    print(
        "corr median UMI / delta R2 =",
        output["correlation_depth_delta_r2"],
    )
    print(
        "corr detected features / delta R2 =",
        output["correlation_detected_features_delta_r2"],
    )


if __name__ == "__main__":
    main()
