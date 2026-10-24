#!/usr/bin/env python3
"""Small, non-inferential program and neighbourhood audit for GSE301973."""

import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "public_validation" / "GSE301973" / "extracted"
OUTPUT = ROOT / "data" / "public_validation" / "GSE301973_small_analysis.json"

MODULES = {
    "egfr_bypass": ["EGFR", "ERBB2", "ERBB3", "MET", "AXL", "FGFR1", "IGF1R"],
    "emt": ["VIM", "CDH2", "ZEB1", "SNAI1", "SNAI2"],
    "stem_like": ["ALDH1A1", "SOX2", "PROM1", "CD44"],
    "spp1_cd44": ["SPP1", "CD44"],
    "proliferation": ["MKI67", "PCNA", "TOP2A", "CDK1"],
}


def decode(values):
    return [value.decode() if isinstance(value, bytes) else value for value in values]


def load_slide(slide):
    path = BASE / slide / f"{slide}_stardist_integrated.h5ad"
    with h5py.File(path, "r") as handle:
        shape = [int(value) for value in handle["X"].attrs["shape"]]
        coords = handle["obsm"]["spatial"][:]
        cell_types = decode(handle["obs"]["cell_type"]["categories"][:])
        codes = handle["obs"]["cell_type"]["codes"][:]
        labels = np.asarray([cell_types[int(code)] for code in codes])
        genes = decode(handle["var"]["_index"][:])
        gene_index = {gene: index for index, gene in enumerate(genes)}
        selected = sorted({
            gene for values in MODULES.values() for gene in values if gene in gene_index
        })
        selected_indices = [gene_index[gene] for gene in selected]
        matrix = np.zeros((shape[0], len(selected)), dtype=np.float32)
        data = handle["X"]["data"][:]
        indices = handle["X"]["indices"][:]
        indptr = handle["X"]["indptr"][:]
        for row in range(shape[0]):
            start, stop = int(indptr[row]), int(indptr[row + 1])
            if stop > start:
                row_indices = indices[start:stop]
                values = data[start:stop]
                positions = np.searchsorted(row_indices, selected_indices)
                valid = (positions < len(row_indices))
                valid &= row_indices[np.minimum(positions, len(row_indices) - 1)] == np.asarray(selected_indices)
                matrix[row, np.where(valid)[0]] = values[positions[valid]]
    totals = matrix.sum(axis=1)
    normalized = np.log1p(matrix / np.maximum(totals[:, None], 1.0) * 1e4)
    scores = {
        name: normalized[:, [selected.index(gene) for gene in genes if gene in selected]].mean(axis=1)
        for name, genes in MODULES.items()
    }
    region_ids = KMeans(n_clusters=2, random_state=42, n_init=10).fit_predict(coords)
    result = {"slide": slide, "cells": int(len(coords)), "regions": []}
    for region in [0, 1]:
        mask = region_ids == region
        tree = cKDTree(coords[mask])
        distances, indices = tree.query(coords[mask], k=min(7, int(mask.sum())))
        nearest_types = labels[mask][indices[:, 1:]]
        result["regions"].append({
            "candidate_region": int(region),
            "cells": int(mask.sum()),
            "centroid": coords[mask].mean(axis=0).tolist(),
            "cell_type_counts": {
                key: int(value)
                for key, value in Counter(labels[mask]).items()
            },
            "module_means": {
                key: float(values[mask].mean())
                for key, values in scores.items()
            },
            "nearest_neighbour_type_fractions": {
                key: float(np.mean(nearest_types == key))
                for key in sorted(set(labels))
            },
        })
    return result


def main():
    output = {
        "accession": "GSE301973",
        "analysis_type": "descriptive pretreatment spatial program audit",
        "patient_assignment_verified": False,
        "warning": (
            "Regions are coordinate-derived candidate partitions. They are not "
            "assigned to long/short PFS until the sample-position mapping is verified."
        ),
        "modules": MODULES,
        "slides": [load_slide(slide) for slide in ["A", "C"]],
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
