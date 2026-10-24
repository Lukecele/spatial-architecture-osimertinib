#!/usr/bin/env python3
"""Audit whether published molecular anchors are present in current inputs.

This prevents importing a microRNA, methylation marker, or pathway claim when
the bundled spatial matrices do not measure it.
"""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "clinical_pan_cancer_cache.npz"
OUTPUT = ROOT / "results" / "MOLECULAR_ANCHOR_AUDIT.json"

ANCHORS = {
    "melanoma_MITF_APAF1_resistance": {
        "cohort": "melanoma",
        "required_transcripts": ["MITF", "APAF1"],
        "source": "Carotenuto et al.; MITF/APAF-1 resistance axis",
        "status_if_complete": "candidate spatial anchor, not causal validation",
    },
    "PDAC_miR1307_CLIC5_FOLFIRINOX": {
        "cohort": None,
        "required_transcripts": ["MIR1307", "CLIC5"],
        "source": "Carotenuto et al.; miR-1307 and FOLFIRINOX response",
        "status_if_complete": "candidate molecular anchor, not present in current cohorts",
    },
    "biliary_miR181_SIRT1_drug_response": {
        "cohort": None,
        "required_transcripts": ["MIR181A1", "SIRT1"],
        "source": "Carotenuto et al.; miR-181/SIRT1 drug response",
        "status_if_complete": "candidate molecular anchor, not present in current cohorts",
    },
}


def available_genes(cache, cohort):
    prefix = f"{cohort}_gene_"
    return {key[len(prefix):] for key in cache.files if key.startswith(prefix)}


def main():
    with np.load(CACHE) as cache:
        cohorts = {
            cohort: sorted(available_genes(cache, cohort))
            for cohort in ("breast", "melanoma", "lung")
        }

    results = {}
    for name, anchor in ANCHORS.items():
        cohort = anchor["cohort"]
        measured = cohorts[cohort] if cohort else []
        present = [gene for gene in anchor["required_transcripts"] if gene in measured]
        missing = [gene for gene in anchor["required_transcripts"] if gene not in measured]
        results[name] = {
            "cohort": cohort,
            "source": anchor["source"],
            "required_features": anchor["required_transcripts"],
            "present_features": present,
            "missing_features": missing,
            "status": (
                anchor["status_if_complete"]
                if not missing
                else "not identifiable from current cached spatial matrix"
            ),
        }

    output = {
        "input": str(CACHE),
        "measured_data_type": (
            "cached gene-level spatial expression and coordinates; "
            "no microRNA counts and no methylome/CpG matrix"
        ),
        "cohort_feature_counts": {cohort: len(genes) for cohort, genes in cohorts.items()},
        "anchors": results,
        "current_action": (
            "Do not claim microRNA or methylome validation. Use the available "
            "melanoma transcript markers only as exploratory state descriptors; "
            "obtain matched miRNA/methylome or perturbation data before using "
            "the published anchors mechanistically."
        ),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
