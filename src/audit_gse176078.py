#!/usr/bin/env python3
"""Audit public GSE176078 files before model validation."""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "public_validation" / "GSE176078" / "extracted"
OUTPUT = ROOT / "data" / "public_validation" / "GSE176078_audit.json"

MARKERS = {
    "tumour_state": ["EPCAM", "KRT8", "KRT18", "KRT5", "TP63", "VIM", "AXL", "CDH2"],
    "outcome": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "hypoxia": ["HIF1A", "VEGFA", "CA9"],
    "vascular": ["PECAM1", "VWF"],
}


def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def read_features(path):
    with path.open() as handle:
        return [line.rstrip("\n").split("\t") for line in handle]


def read_positions(path):
    with path.open() as handle:
        return [line.rstrip("\n").split(",") for line in handle]


def inspect_sample(sample):
    metadata = read_rows(BASE / "metadata" / f"{sample}_metadata.csv")
    positions = read_positions(
        BASE / "spatial" / f"{sample}_spatial" / "tissue_positions_list.csv"
    )
    features = read_features(
        BASE / "filtered_count_matrices" / f"{sample}_filtered_count_matrix" / "features.tsv.gz"
    )
    names = {row[1] if len(row) > 1 else row[0] for row in features if row}
    matrix_path = (
        BASE / "filtered_count_matrices" / f"{sample}_filtered_count_matrix" / "matrix.mtx.gz"
    )
    with matrix_path.open("rb") as matrix_handle:
        shape = mmread(matrix_handle).shape
    metadata_barcodes = {row.get("", "") for row in metadata}
    position_barcodes = {row[0] for row in positions if row}
    marker_presence = {
        group: [gene for gene in genes if gene in names]
        for group, genes in MARKERS.items()
    }
    return {
        "sample": sample,
        "patient_ids": sorted({row.get("patientid", "") for row in metadata}),
        "subtypes": sorted({row.get("subtype", "") for row in metadata}),
        "metadata_rows": len(metadata),
        "spatial_rows": len(positions),
        "metadata_spatial_barcode_overlap": len(metadata_barcodes & position_barcodes),
        "matrix_shape_genes_by_spots": list(shape),
        "marker_presence": marker_presence,
        "classification_counts": {
            label: sum(row.get("Classification", "") == label for row in metadata)
            for label in sorted({row.get("Classification", "") for row in metadata})
        },
    }


def main():
    samples = sorted(
        path.name.removesuffix("_metadata.csv")
        for path in (BASE / "metadata").glob("*_metadata.csv")
    )
    result = {
        "accession": "GSE176078",
        "samples": [inspect_sample(sample) for sample in samples],
        "unique_patient_ids": samples,
        "limitations": [
            "GEO raw human sequencing files are controlled; this audit uses public processed matrices.",
            "Patient-level validation requires grouping all sections from the same donor.",
            "Marker presence does not validate biological state or causality.",
        ],
    }
    result["unique_patient_ids"] = sorted({
        patient
        for sample in result["samples"]
        for patient in sample["patient_ids"]
    })
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(OUTPUT)
    print(json.dumps({
        "samples": len(result["samples"]),
        "unique_patient_ids": len(result["unique_patient_ids"]),
        "patients": result["unique_patient_ids"],
    }, indent=2))


if __name__ == "__main__":
    main()
