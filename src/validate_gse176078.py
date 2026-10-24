#!/usr/bin/env python3
"""Patient-held-out validation on the public GSE176078 breast cohort."""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "public_validation" / "GSE176078" / "extracted"
OUTPUT = ROOT / "data" / "public_validation" / "GSE176078_validation.json"
RADIUS_UM = 150.0

GENE_SPEC = {
    "positive": ["EPCAM", "KRT8", "CDH1"],
    "resistant": ["VIM", "CDH2", "AXL"],
    "caf": ["ACTA2", "COL1A1", "FAP"],
    "outcome": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "hypoxia": ["HIF1A", "VEGFA", "CA9"],
    "vascular": ["PECAM1", "VWF"],
}
EXCLUDED_CLASSIFICATIONS = {"Artefact", "Necrosis"}


def read_metadata(sample):
    path = BASE / "metadata" / f"{sample}_metadata.csv"
    with path.open() as handle:
        return {row[""]: row for row in csv.DictReader(handle)}


def read_coordinates(sample):
    path = BASE / "spatial" / f"{sample}_spatial" / "tissue_positions_list.csv"
    coords = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split(",")
            if fields:
                coords[fields[0]] = [float(fields[4]), float(fields[5])]
    return coords


def read_features(sample):
    path = BASE / "filtered_count_matrices" / f"{sample}_filtered_count_matrix" / "features.tsv.gz"
    with path.open() as handle:
        return [line.rstrip("\n").split("\t")[0] for line in handle]


def load_sample(sample):
    metadata = read_metadata(sample)
    coordinates = read_coordinates(sample)
    features = read_features(sample)
    genes = {gene: index for index, gene in enumerate(features)}
    selected = sorted({
        gene for group in GENE_SPEC.values() for gene in group if gene in genes
    })
    rows = np.array([genes[gene] for gene in selected])
    matrix_path = (
        BASE / "filtered_count_matrices" / f"{sample}_filtered_count_matrix" / "matrix.mtx.gz"
    )
    with matrix_path.open("rb") as handle:
        matrix = mmread(handle).tocsr()[rows].toarray()
    counts = {gene: matrix[index] for index, gene in enumerate(selected)}
    barcodes = sorted(set(metadata) & set(coordinates))
    barcode_indices = {
        barcode: index for index, barcode in enumerate(
            _read_barcodes(sample)
        )
    }
    columns = np.array([barcode_indices[barcode] for barcode in barcodes])
    total = _read_total_counts(sample, matrix.shape[1])
    normalized = {
        gene: np.log1p((values[columns] / np.maximum(total[columns], 1.0)) * 1e4)
        for gene, values in counts.items()
    }

    def module(group):
        present = [normalized[gene] for gene in GENE_SPEC[group] if gene in normalized]
        return np.mean(present, axis=0)

    valid = np.array([
        metadata[barcode].get("Classification", "") not in EXCLUDED_CLASSIFICATIONS
        for barcode in barcodes
    ])
    coords = np.asarray([coordinates[barcode] for barcode in barcodes])[valid]
    scores = {group: module(group)[valid] for group in GENE_SPEC}
    expression = np.vstack([scores["positive"], scores["resistant"], scores["caf"]])
    expression -= expression.max(axis=0, keepdims=True)
    frequencies = np.exp(expression)
    frequencies /= frequencies.sum(axis=0, keepdims=True)
    environment = np.column_stack([scores["hypoxia"], scores["vascular"]])
    return coords, frequencies.T, scores["outcome"], environment, int(valid.sum())


def _read_barcodes(sample):
    path = BASE / "filtered_count_matrices" / f"{sample}_filtered_count_matrix" / "barcodes.tsv.gz"
    with path.open() as handle:
        return [line.rstrip("\n") for line in handle]


def _read_total_counts(sample, length):
    path = BASE / "filtered_count_matrices" / f"{sample}_filtered_count_matrix" / "matrix.mtx.gz"
    with path.open("rb") as handle:
        matrix = mmread(handle).tocsr()
    return np.asarray(matrix.sum(axis=0)).ravel()[:length]


def neighbourhood(coords, frequencies):
    scale = 100.0 / np.median(cKDTree(coords).query(coords, k=2)[0][:, 1])
    tree = cKDTree(coords * scale)
    local = np.empty_like(frequencies)
    for index, neighbours in enumerate(tree.query_ball_point(tree.data, RADIUS_UM)):
        neighbours = [value for value in neighbours if value != index]
        local[index] = frequencies[neighbours].mean(axis=0) if neighbours else frequencies[index]
    return local


def fit_predict(train, test):
    train_focal = np.vstack([row[1] for row in train])
    train_context = np.vstack([np.column_stack([row[1], row[5]]) for row in train])
    train_outcome = np.concatenate([row[2] for row in train])
    train_environment = np.vstack([row[3] for row in train])
    test_focal = np.vstack([row[1] for row in test])
    test_context = np.vstack([np.column_stack([row[1], row[5]]) for row in test])
    test_outcome = np.concatenate([row[2] for row in test])
    env_model = LinearRegression().fit(train_environment, train_outcome)
    train_residual = train_outcome - env_model.predict(train_environment)
    test_residual = test_outcome - env_model.predict(np.vstack([row[3] for row in test]))
    focal_model = LinearRegression(fit_intercept=False).fit(train_focal, train_residual)
    context_model = LinearRegression(fit_intercept=False).fit(train_context, train_residual)
    return {
        "focal_r2": float(r2_score(test_residual, focal_model.predict(test_focal))),
        "focal_plus_neighbourhood_r2": float(
            r2_score(test_residual, context_model.predict(test_context))
        ),
    }


def bootstrap(values, seed=42, replicates=10000):
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(replicates, len(values)), replace=True).mean(axis=1)
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))]


def main():
    samples = sorted(
        path.name.removesuffix("_metadata.csv")
        for path in (BASE / "metadata").glob("*_metadata.csv")
    )
    data = {}
    for sample in samples:
        coords, frequencies, outcome, environment, spots = load_sample(sample)
        data[sample] = (sample, frequencies, outcome, environment, spots, neighbourhood(coords, frequencies))

    rows = []
    for held_out in samples:
        result = fit_predict(
            [data[sample] for sample in samples if sample != held_out],
            [data[held_out]],
        )
        result.update({"held_out_patient": held_out, "spots": data[held_out][4]})
        result["context_minus_focal_r2"] = (
            result["focal_plus_neighbourhood_r2"] - result["focal_r2"]
        )
        rows.append(result)
    deltas = np.asarray([row["context_minus_focal_r2"] for row in rows])
    output = {
        "accession": "GSE176078",
        "radius_um": RADIUS_UM,
        "marker_spec": GENE_SPEC,
        "excluded_classifications": sorted(EXCLUDED_CLASSIFICATIONS),
        "patients": rows,
        "mean_patient_delta_r2": float(deltas.mean()),
        "median_patient_delta_r2": float(np.median(deltas)),
        "patients_with_positive_delta": int(np.sum(deltas > 0)),
        "bootstrap_95ci_mean_delta_r2": bootstrap(deltas),
        "interpretation": (
            "Patient-held-out exploratory validation using public processed "
            "matrices; not causal or clinical validation."
        ),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({
        "mean_patient_delta_r2": output["mean_patient_delta_r2"],
        "patients_with_positive_delta": output["patients_with_positive_delta"],
        "bootstrap_95ci_mean_delta_r2": output["bootstrap_95ci_mean_delta_r2"],
    }, indent=2))


if __name__ == "__main__":
    main()
