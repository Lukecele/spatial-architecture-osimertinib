#!/usr/bin/env python3
"""Tumour-adapted patient-held-out validation for GSE189487 lung ST."""

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "public_validation" / "GSE189487" / "extracted"
OUTPUT = ROOT / "data" / "public_validation" / "GSE189487_validation.json"
RADIUS_UM = 150.0

GENE_SPEC = {
    "epithelial_alveolar": ["EPCAM", "KRT8", "SCGB1A1", "SFTPC"],
    "invasive_mesenchymal": ["VIM", "AXL", "CDH2", "KRT19"],
    "caf": ["COL1A1", "COL3A1", "ACTA2", "FAP"],
    "immune": ["CD68", "CD3D", "CD8A", "MS4A1", "NKG7"],
    "outcome": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "hypoxia": ["HIF1A", "VEGFA", "CA9"],
    "vascular": ["PECAM1", "VWF"],
}


def samples():
    return sorted(
        path.name.removesuffix("_features.tsv.gz")
        for path in BASE.glob("*_features.tsv.gz")
    )


def read_lines(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def load_sample(sample):
    features = [line.split("\t")[1] for line in read_lines(BASE / f"{sample}_features.tsv.gz")]
    feature_index = {gene: index for index, gene in enumerate(features)}
    barcodes = read_lines(BASE / f"{sample}_barcodes.tsv.gz")
    positions = {}
    for line in read_lines(BASE / f"{sample}_tissue_positions_list.csv.gz"):
        fields = line.split(",")
        positions[fields[0]] = [float(fields[4]), float(fields[5])]
    keep_barcodes = [barcode for barcode in barcodes if barcode in positions]
    columns = np.asarray([barcodes.index(barcode) for barcode in keep_barcodes])
    selected = sorted({
        gene for genes in GENE_SPEC.values() for gene in genes
        if gene in feature_index
    })
    row_indices = [feature_index[gene] for gene in selected]
    matrix_path = BASE / f"{sample}_matrix.mtx.gz"
    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr()
    counts = matrix[row_indices][:, columns].toarray()
    total = np.asarray(matrix[:, columns].sum(axis=0)).ravel()
    normalized = {
        gene: np.log1p((counts[index] / np.maximum(total, 1.0)) * 1e4)
        for index, gene in enumerate(selected)
    }

    def module(name):
        return np.mean([normalized[gene] for gene in GENE_SPEC[name]], axis=0)

    scores = {name: module(name) for name in GENE_SPEC}
    composition = np.vstack([
        scores["epithelial_alveolar"],
        scores["invasive_mesenchymal"],
        scores["caf"],
        scores["immune"],
    ])
    composition -= composition.max(axis=0, keepdims=True)
    composition = np.exp(composition)
    composition /= composition.sum(axis=0, keepdims=True)
    coords = np.asarray([positions[barcode] for barcode in keep_barcodes])
    return coords, composition.T, scores["outcome"], np.column_stack(
        [scores["hypoxia"], scores["vascular"]]
    )


def neighbourhood(coords, composition):
    scale = 100.0 / np.median(cKDTree(coords).query(coords, k=2)[0][:, 1])
    tree = cKDTree(coords * scale)
    local = np.empty_like(composition)
    for index, neighbours in enumerate(tree.query_ball_point(tree.data, RADIUS_UM)):
        neighbours = [value for value in neighbours if value != index]
        local[index] = composition[neighbours].mean(axis=0) if neighbours else composition[index]
    return local


def validate(train, test):
    focal_train = np.vstack([item[1] for item in train])
    context_train = np.vstack([
        np.column_stack([item[1], item[4]]) for item in train
    ])
    outcome_train = np.concatenate([item[2] for item in train])
    environment_train = np.vstack([item[3] for item in train])
    focal_test = np.vstack([item[1] for item in test])
    context_test = np.vstack([
        np.column_stack([item[1], item[4]]) for item in test
    ])
    outcome_test = np.concatenate([item[2] for item in test])
    environment_test = np.vstack([item[3] for item in test])
    environment_model = LinearRegression().fit(environment_train, outcome_train)
    train_residual = outcome_train - environment_model.predict(environment_train)
    test_residual = outcome_test - environment_model.predict(environment_test)
    focal_model = LinearRegression(fit_intercept=False).fit(focal_train, train_residual)
    context_model = LinearRegression(fit_intercept=False).fit(context_train, train_residual)
    focal_r2 = r2_score(test_residual, focal_model.predict(focal_test))
    context_r2 = r2_score(test_residual, context_model.predict(context_test))
    return {
        "focal_r2": float(focal_r2),
        "focal_plus_neighbourhood_r2": float(context_r2),
        "context_minus_focal_r2": float(context_r2 - focal_r2),
        "spots": int(len(outcome_test)),
    }


def bootstrap(values, seed=42, replicates=10000):
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(replicates, len(values)), replace=True).mean(axis=1)
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))]


def main():
    data = {}
    for sample in samples():
        coords, composition, outcome, environment = load_sample(sample)
        data[sample] = (
            sample, composition, outcome, environment,
            neighbourhood(coords, composition),
        )
    rows = []
    for held_out in samples():
        result = validate(
            [data[sample] for sample in samples() if sample != held_out],
            [data[held_out]],
        )
        result["held_out_sample"] = held_out
        rows.append(result)
    deltas = np.asarray([row["context_minus_focal_r2"] for row in rows])
    output = {
        "accession": "GSE189487",
        "radius_um": RADIUS_UM,
        "samples": samples(),
        "marker_spec": GENE_SPEC,
        "patients_assumed_one_per_spatial_sample": True,
        "patients": rows,
        "mean_patient_delta_r2": float(deltas.mean()),
        "median_patient_delta_r2": float(np.median(deltas)),
        "samples_with_positive_delta": int(np.sum(deltas > 0)),
        "bootstrap_95ci_mean_delta_r2": bootstrap(deltas),
        "interpretation": (
            "Lung-adapted patient-held-out validation using public processed "
            "matrices; not causal or clinical validation."
        ),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({
        "mean_patient_delta_r2": output["mean_patient_delta_r2"],
        "samples_with_positive_delta": output["samples_with_positive_delta"],
        "bootstrap_95ci_mean_delta_r2": output["bootstrap_95ci_mean_delta_r2"],
    }, indent=2))


if __name__ == "__main__":
    main()
