#!/usr/bin/env python3
"""Audit spatial proxies relevant to EGFR-TKI resistance in public LUAD ST.

GSE189487 is not a treated/response cohort. This script therefore reports
candidate transcriptional resistance proxies only and never labels them as
osimertinib response or resistance measurements.
"""

import gzip
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "public_validation" / "GSE189487" / "extracted"
OUTPUT = ROOT / "data" / "public_validation" / "GSE189487_OSIMERTINIB_RELEVANCE.json"
RADII = [55.0, 100.0, 150.0, 200.0, 300.0, 500.0]

MODULES = {
    "egfr_axis": ["EGFR", "ERBB2", "ERBB3", "MET"],
    "bypass_resistance": ["MET", "AXL", "FGFR1", "IGF1R", "ERBB2"],
    "emt_resistance": ["VIM", "AXL", "CDH2", "ZEB1", "SNAI1", "SNAI2"],
    "stem_like": ["ALDH1A1", "SOX2", "PROM1", "CD44"],
    "proliferation": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "epithelial_alveolar": ["EPCAM", "KRT8", "SCGB1A1", "SFTPC"],
}


def sample_names():
    return sorted(
        path.name.removesuffix("_features.tsv.gz")
        for path in BASE.glob("*_features.tsv.gz")
    )


def read_gzip_lines(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def load_sample(sample):
    features = [line.split("\t")[1] for line in read_gzip_lines(BASE / f"{sample}_features.tsv.gz")]
    feature_index = {gene: index for index, gene in enumerate(features)}
    barcodes = read_gzip_lines(BASE / f"{sample}_barcodes.tsv.gz")
    positions = {
        fields[0]: [float(fields[4]), float(fields[5])]
        for fields in (
            line.split(",")
            for line in read_gzip_lines(BASE / f"{sample}_tissue_positions_list.csv.gz")
        )
    }
    kept = [barcode for barcode in barcodes if barcode in positions]
    columns = np.array([barcodes.index(barcode) for barcode in kept])
    selected = sorted({
        gene for genes in MODULES.values() for gene in genes
        if gene in feature_index
    })
    with gzip.open(BASE / f"{sample}_matrix.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsr()
    rows = [feature_index[gene] for gene in selected]
    values = matrix[rows][:, columns].toarray()
    totals = np.asarray(matrix[:, columns].sum(axis=0)).ravel()
    normalized = {
        gene: np.log1p(values[index] / np.maximum(totals, 1.0) * 1e4)
        for index, gene in enumerate(selected)
    }
    modules = {
        name: np.mean([normalized[gene] for gene in genes if gene in normalized], axis=0)
        for name, genes in MODULES.items()
    }
    return np.asarray([positions[barcode] for barcode in kept]), modules


def moran(values, neighbours):
    centered = values - np.mean(values)
    denominator = np.sum(centered**2)
    weights = sum(len(row) for row in neighbours)
    numerator = sum(
        centered[index] * centered[j]
        for index, row in enumerate(neighbours)
        for j in row
    )
    return float(len(values) / weights * numerator / denominator) if weights and denominator else 0.0


def radius_neighbours(coords, radius):
    scale = 100.0 / np.median(cKDTree(coords).query(coords, k=2)[0][:, 1])
    tree = cKDTree(coords * scale)
    return [
        [j for j in indices if j != i]
        for i, indices in enumerate(tree.query_ball_point(tree.data, radius))
    ]


def main():
    sample_results = []
    for sample in sample_names():
        coords, modules = load_sample(sample)
        sample_result = {"sample": sample, "spots": len(coords), "radii": {}}
        for radius in RADII:
            neighbours = radius_neighbours(coords, radius)
            sample_result["radii"][str(radius)] = {
                "moran_emt_resistance": moran(modules["emt_resistance"], neighbours),
                "moran_bypass_resistance": moran(modules["bypass_resistance"], neighbours),
                "moran_proliferation": moran(modules["proliferation"], neighbours),
                "mean_neighbours": float(np.mean([len(row) for row in neighbours])),
            }
        sample_result["module_genes"] = MODULES
        sample_results.append(sample_result)
    OUTPUT.write_text(json.dumps({
        "accession": "GSE189487",
        "cohort_type": "untreated LUAD progression spatial transcriptomics",
        "osimertinib_response_measured": False,
        "interpretation": (
            "Candidate resistance-related transcriptional proxy audit only; "
            "not a drug-response or resistance validation."
        ),
        "samples": sample_results,
    }, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
