# Contribution and novelty scope

## Purpose

This document defines the strongest defensible contribution of the merged
study after a targeted literature check. It is not a systematic review and
does not establish absolute priority over every unpublished, preprint, or
non-English result.

## What is not new by itself

The following ingredients are established or actively used in spatial
transcriptomics and spatial cancer modelling:

- nearest-neighbour or fixed-radius spatial graphs;
- local neighbourhood composition;
- spatial autocorrelation and graph summaries;
- multiscale spatial analysis;
- regression using spatial context;
- patient-held-out validation as a robustness strategy;
- evolutionary-game models of cancer populations;
- spatial statistics used to classify games in simulated agent-based systems.

Relevant comparators include:

- [SpNeigh](https://www.biorxiv.org/content/10.1101/2025.11.07.687304v1),
  which models neighbourhoods, distance weighting, boundaries and spatial
  gradients;
- [MCIST](https://pmc.ncbi.nlm.nih.gov/articles/PMC12667549/), which studies
  multiscale cell-cell interactions;
- [SpatialZoomer](https://www.nature.com/articles/s41551-026-01800-w), which
  uses multiscale spectral spatial features;
- [Leither et al.](https://www.biorxiv.org/content/10.1101/2025.07.09.664005v1),
  which uses spatial statistics to infer game classes in an agent-based cancer
  model.

## Candidate contribution of this study

The potentially original contribution is the **specific, interpretable
validation composition**:

1. a spot-level molecular outcome proxy;
2. local phenotype composition;
3. explicit radius-based neighbourhood composition;
4. a bilinear interaction model compared with a non-game spatial baseline;
5. environmental adjustment for hypoxia and vascular markers;
6. spatial-block out-of-fold testing;
7. patient-held-out testing and patient-level bootstrap;
8. marker and depth sensitivity audits;
9. explicit negative-result interpretation when the interaction model does
   not reliably improve the baseline.

The novelty claim must therefore be about an **auditable analysis protocol and
its application to human tumour spatial transcriptomics**, not about inventing
neighbourhood analysis, graph methods, multiscale modelling, or evolutionary
game theory separately.

## Recommended central claim

> We present an interpretable and falsifiable pipeline for testing whether
> multiscale neighbourhood composition adds predictive information beyond
> spot composition for a molecular proliferation proxy in tumour spatial
> transcriptomics. The pipeline separates non-game spatial context from
> bilinear interaction terms and evaluates them with spatial and patient-level
> holdouts.

## Claims to avoid

Do not claim that the study:

- proves a tumour evolutionary game;
- directly measures biological payoff matrices;
- identifies genetic clones from transcriptomic states;
- reconstructs temporal replicator dynamics from one snapshot;
- validates an adaptive treatment policy;
- is the first use of any spatial graph, neighbourhood, or multiscale method;
- has established absolute methodological priority without a systematic review.

## Current empirical conclusion

The protocol is technically reproducible and the data contain spatial
structure. The public patient-held-out multiscale analysis shows small positive
mean game-minus-baseline effects from 150 to 300 µm, but all patient bootstrap
intervals include zero. Thus the current data support the protocol and an
exploratory scale-localized signal, not a validated incremental game effect.

This negative or inconclusive result is part of the contribution: the protocol
tests whether a game interpretation earns incremental predictive value rather
than assuming that interaction coefficients are biological payoffs.

## Publication positioning

The most defensible positioning is a methods/benchmarking paper on
interpretable spatial-context modelling and validation. A game-theoretic
interpretation can remain a prespecified candidate model or secondary
analysis, provided the non-game baselines and uncertainty results remain
primary.
