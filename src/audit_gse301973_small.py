#!/usr/bin/env python3
"""Minimal integrity and biological-availability audit for GSE301973."""

import json
from collections import Counter
from pathlib import Path

import h5py


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "public_validation" / "GSE301973" / "extracted"
OUTPUT = ROOT / "data" / "public_validation" / "GSE301973_small_audit.json"


def values(group, key):
    node = group[key]
    if "categories" in node:
        categories = [x.decode() for x in node["categories"][:]]
        codes = node["codes"][:]
        return [categories[int(code)] for code in codes]
    return [x.decode() if isinstance(x, bytes) else x for x in node[:]]


def inspect_slide(slide):
    path = BASE / slide / f"{slide}_stardist_integrated.h5ad"
    with h5py.File(path, "r") as handle:
        shape = [int(value) for value in handle["X"].attrs["shape"]]
        spatial = handle["obsm"]["spatial"][:]
        genes = [x.decode() for x in handle["var"]["_index"][:]]
        cell_types = values(handle["obs"], "cell_type")
        return {
            "slide": slide,
            "cells": shape[0],
            "genes": shape[1],
            "spatial_dimensions": [int(value) for value in spatial.shape],
            "spatial_range": {
                "min": spatial.min(axis=0).tolist(),
                "max": spatial.max(axis=0).tolist(),
            },
            "cell_type_counts": {
                key: int(value) for key, value in Counter(cell_types).items()
            },
            "marker_presence": {
                gene: gene in genes
                for gene in [
                    "EGFR", "MET", "AXL", "FGFR1", "IGF1R", "VIM", "CDH2",
                    "ZEB1", "SNAI1", "SNAI2", "ALDH1A1", "SOX2", "PROM1",
                    "CD44", "SPP1", "MKI67", "PCNA", "TOP2A", "CDK1",
                    "EPCAM", "KRT8", "SCGB1A1", "SFTPC", "COL1A1", "FAP",
                ]
            },
        }


def main():
    result = {
        "accession": "GSE301973",
        "cohort": "four pretreatment EGFR-mutant NSCLC patients",
        "osimertinib_response_measured": False,
        "slides": [inspect_slide(slide) for slide in ["A", "C"]],
        "metadata_boundary": (
            "PFS groups are available at patient level; treatment is pretreatment. "
            "PFS association is exploratory because n=4."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "slides": len(result["slides"]),
        "cells": [slide["cells"] for slide in result["slides"]],
        "all_markers_present": all(
            all(slide["marker_presence"].values()) for slide in result["slides"]
        ),
    }, indent=2))


if __name__ == "__main__":
    main()
