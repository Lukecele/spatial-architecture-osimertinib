#!/usr/bin/env python3
"""Pre-specified multiscale held-out audit for the lung cohort."""

import json
from pathlib import Path

import numpy as np

import validate_gse189487 as lung


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "public_validation" / "GSE189487_multiscale.json"
RADII = [55.0, 100.0, 150.0, 200.0, 300.0, 500.0]


def main():
    samples = lung.samples()
    loaded = {
        sample: lung.load_sample(sample)
        for sample in samples
    }
    results = {}
    for radius in RADII:
        lung.RADIUS_UM = radius
        data = {}
        for sample, (coords, composition, outcome, environment) in loaded.items():
            data[sample] = (
                sample,
                composition,
                outcome,
                environment,
                lung.neighbourhood(coords, composition),
            )
        rows = []
        for held_out in samples:
            row = lung.validate(
                [data[sample] for sample in samples if sample != held_out],
                [data[held_out]],
            )
            row["held_out_sample"] = held_out
            rows.append(row)
        deltas = np.asarray([row["context_minus_focal_r2"] for row in rows])
        results[str(radius)] = {
            "radius_um": radius,
            "samples": rows,
            "mean_delta_r2": float(deltas.mean()),
            "median_delta_r2": float(np.median(deltas)),
            "positive_samples": int(np.sum(deltas > 0)),
        }
    OUTPUT.write_text(json.dumps({
        "accession": "GSE189487",
        "radii_um": RADII,
        "results": results,
        "interpretation": (
            "Multiscale lung-adapted held-out diagnostic; radius was not "
            "selected as a confirmatory result from these same samples."
        ),
    }, indent=2) + "\n")
    for radius, result in results.items():
        print(radius, result["mean_delta_r2"], result["positive_samples"])


if __name__ == "__main__":
    main()
