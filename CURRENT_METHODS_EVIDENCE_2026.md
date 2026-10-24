# Current methods evidence register (October 2026)

This register records current method families found in the online literature
search. A method being published or popular does not make it appropriate for
the present data.

## Evidence relevant to the protocol

### Patient-level inference and pseudoreplication

The recent comparative literature continues to warn that treating spots as
independent replicates inflates evidence. The protocol therefore uses grouped
patient-held-out validation, patient-level bootstrap, and explicit null models.
The [PLOS Computational Biology comparative study](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013956)
is recorded as methodological context, not as a reason to adopt one package
without checking its assumptions.

### Graph and AI models

Recent spatial methods include graph neural networks and multimodal graph
learning, including [SpaInGNN](https://www.sciencedirect.com/science/article/pii/S1046202324002433),
[STING](https://www.biorxiv.org/content/biorxiv/early/2025/04/26/2025.04.23.650236.full.pdf),
and [AugGCL](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013912).
These are candidate P2/P3 comparisons. They require nested grouped validation,
leakage checks, and a sample-size justification before they can be used for a
primary result.

### Foundation models

[Novae](https://www.nature.com/articles/s41592-025-02899-6) is recorded as a
recent graph-based foundation-model direction. Access, licensing, input
compatibility, pretrained-data provenance and patient leakage must be checked
before use. A foundation model cannot replace a simple baseline or create
biological validation.

### Spatial DE and post-imputation inference

The current protocol avoids default pooled spot-level differential testing.
Pseudobulk, mixed-effects/GEE approaches, and uncertainty-aware post-imputation
testing are possible P1/P2 modules. [TIDEST](https://www.biorxiv.org/content/10.64898/2026.06.19.733432v1)
is recorded as a research direction, not a dependency of the current assay.

## Status categories

| Category | Meaning |
|---|---|
| Implemented | Code exists locally and targeted tests pass |
| Ready candidate | Dependencies/data are known and can be added after P0 |
| Research-only | Interesting method, but prerequisites or validation are currently missing |
| Excluded | Cannot be identified from the available measurements |

## Current status

- Implemented: explicit boolean QC masks, log-CPM, z-score matrices,
  Benjamini-Hochberg adjustment, transparent normal-normal Bayesian effect
  summary, grouped validation and patient bootstrap, and patient-level Visium HD
  subcellular spatial neighborhood and PFS program analysis on GSE301973 (33/33 tests passing).
- Ready candidate: Scanpy/Squidpy, cell2location/RCTD/Tangram, PyMC/NumPyro,
  graph methods, pathology integration, and mixed-effects/GEE.
- Research-only: large GNNs, foundation models, ligand-receptor prioritisation,
  and causal/dynamic models.
- Excluded without new measurements: miRNA mechanisms, methylome mechanisms,
  protein validation, genetic clonality, causal communication, and clinical
  efficacy.
