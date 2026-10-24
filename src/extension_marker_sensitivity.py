#!/usr/bin/env python3
"""Assess marker-set sensitivity on the public paired NSCLC extension."""

import json
from pathlib import Path

import numpy as np

from extension_pilot_baselines import GENE_SPEC, block_cv, load_sample


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "extension_pilot"
SAMPLES = ("P10_B1", "P10_T1", "P11_B1", "P11_T1",
           "P15_B1", "P15_T1", "P16_B1", "P16_T1",
           "P17_B1", "P17_T1", "P19_B1", "P19_T1",
           "P24_B1", "P24_T1", "P25_B1", "P25_T1")

VARIANTS = {
    "reference": GENE_SPEC,
    "remove_vim_from_resistant": {
        **GENE_SPEC,
        "T_R": ["AXL", "CDH2", "MMP2"],
    },
    "expanded_epithelial_sensitive": {
        **GENE_SPEC,
        "T_S": ["EPCAM", "KRT8", "KRT18", "CDH1"],
    },
    "expanded_caf": {
        **GENE_SPEC,
        "S_CAF": ["ACTA2", "COL1A1", "FAP", "PDGFRB", "DCN"],
    },
}


def design(coords, frequencies):
    from extension_pilot_baselines import neighborhood_frequencies
    local = neighborhood_frequencies(coords, frequencies)
    game = np.column_stack([
        frequencies[:, row, None] * local
        for row in range(frequencies.shape[1])
    ])
    composition = np.column_stack([frequencies, local])
    return composition, game


def main():
    output = {}
    for name, spec in VARIANTS.items():
        rows = []
        for sample in SAMPLES:
            coords, frequencies, outcome = load_sample(sample, spec)
            composition, game = design(coords, frequencies)
            rows.append({
                "sample": sample,
                "composition_block_cv_r2": block_cv(
                    composition, outcome, coords
                )[0],
                "game_block_cv_r2": block_cv(game, outcome, coords)[0],
                "game_in_sample_r2": float(
                    1.0 - np.sum((outcome - game @ np.linalg.lstsq(
                        game, outcome, rcond=None
                    )[0]) ** 2) / np.sum((outcome - outcome.mean()) ** 2)
                ),
            })
        output[name] = {
            "rows": rows,
            "mean_composition_block_cv_r2": float(
                np.mean([row["composition_block_cv_r2"] for row in rows])
            ),
            "mean_game_block_cv_r2": float(
                np.mean([row["game_block_cv_r2"] for row in rows])
            ),
            "positive_game_advantage_count": int(sum(
                row["game_block_cv_r2"] > row["composition_block_cv_r2"]
                for row in rows
            )),
        }
    destination = PILOT / "marker_sensitivity_results.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(destination)
    for name, values in output.items():
        print(name, values["mean_game_block_cv_r2"],
              values["positive_game_advantage_count"])


if __name__ == "__main__":
    main()
