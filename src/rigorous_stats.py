"""Small, dependency-light statistical primitives for the spatial audits."""

from __future__ import annotations

import math

import numpy as np


def boolean_mask(
    *,
    total_counts: np.ndarray,
    detected_features: np.ndarray,
    tissue: np.ndarray | None = None,
    min_counts: float = 500.0,
    min_features: float = 200.0,
    max_counts: float | None = None,
    max_features: float | None = None,
) -> np.ndarray:
    """Return an explicit QC mask; no sample is silently rescued."""
    total_counts = np.asarray(total_counts, dtype=float)
    detected_features = np.asarray(detected_features, dtype=float)
    if total_counts.shape != detected_features.shape:
        raise ValueError("total_counts and detected_features must have the same shape")
    mask = (total_counts >= min_counts) & (detected_features >= min_features)
    if max_counts is not None:
        mask &= total_counts <= max_counts
    if max_features is not None:
        mask &= detected_features <= max_features
    if tissue is not None:
        tissue = np.asarray(tissue, dtype=bool)
        if tissue.shape != mask.shape:
            raise ValueError("tissue mask must have the same shape as QC arrays")
        mask &= tissue
    return mask


def log_cpm(counts: np.ndarray, scale: float = 10_000.0) -> np.ndarray:
    """Library-size normalisation followed by log1p."""
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2 or np.any(counts < 0):
        raise ValueError("counts must be a non-negative 2-D matrix")
    totals = counts.sum(axis=0)
    if np.any(totals <= 0):
        raise ValueError("every column must contain positive counts")
    return np.log1p(counts / totals * scale)


def zscore_matrix(matrix: np.ndarray, axis: int = 0) -> np.ndarray:
    """Z-transform a matrix, returning zero for constant features."""
    matrix = np.asarray(matrix, dtype=float)
    mean = matrix.mean(axis=axis, keepdims=True)
    scale = matrix.std(axis=axis, keepdims=True, ddof=1)
    return np.divide(matrix - mean, scale, out=np.zeros_like(matrix), where=scale > 0)


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values in the original order."""
    p_values = np.asarray(p_values, dtype=float)
    if np.any((p_values < 0) | (p_values > 1)):
        raise ValueError("p-values must be between zero and one")
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def normal_normal_posterior(
    effects: np.ndarray,
    standard_errors: np.ndarray | float,
    prior_mean: float = 0.0,
    prior_sd: float = 0.1,
) -> dict[str, float]:
    """Conjugate normal-normal posterior for a patient-level effect.

    This is a transparent Bayesian sensitivity summary, not a substitute for
    a hierarchical model when patient-level sample sizes are adequate.
    """
    effects = np.asarray(effects, dtype=float)
    errors = np.broadcast_to(np.asarray(standard_errors, dtype=float), effects.shape)
    if effects.size == 0 or np.any(errors <= 0) or prior_sd <= 0:
        raise ValueError("effects must be non-empty and all standard deviations positive")
    prior_precision = 1.0 / prior_sd**2
    likelihood_precision = np.sum(1.0 / errors**2)
    posterior_variance = 1.0 / (prior_precision + likelihood_precision)
    posterior_mean = posterior_variance * (
        prior_precision * prior_mean + np.sum(effects / errors**2)
    )
    posterior_sd = np.sqrt(posterior_variance)
    probability_positive = 0.5 * (
        1.0 + math.erf(posterior_mean / (posterior_sd * np.sqrt(2.0)))
    )
    return {
        "posterior_mean": float(posterior_mean),
        "posterior_sd": float(posterior_sd),
        "posterior_ci95_low": float(posterior_mean - 1.96 * posterior_sd),
        "posterior_ci95_high": float(posterior_mean + 1.96 * posterior_sd),
        "posterior_probability_positive": float(probability_positive),
    }
