# Method inventory and filtering plan

## Objective

The project will consider a broad set of computational and bioinformatic tools,
but will retain only methods that answer a defined question with available data
and measurable validation. A method is not included in the final claim merely
because it produces a figure, network, score, or p-value.

## Universal analysis contract

Every retained method must specify:

1. input data and data provenance;
2. biological or geometric question;
3. independent validation unit;
4. null or baseline model;
5. uncertainty estimate;
6. failure criterion;
7. interpretation boundary.

The independent unit is the patient/donor whenever available. Spots are
observations within a patient, not independent biological replicates.

## Broad inventory

| Domain | Candidate methods | What it can answer | Priority | Keep only if |
|---|---|---|---|---|
| Quality control | library size, detected genes, mitochondrial/ribosomal fractions, tissue masks, outlier detection | whether measurements are technically usable | P0 | QC is reported before modelling |
| Spatial geometry | coordinate audit, radius/kNN graphs, connected components, edge distance | whether the spatial graph is valid | P0 | coordinates and units are verified |
| Spatial statistics | Moran's I, Geary's C, LISA, variograms, distance decay | whether a feature is spatially structured | P0 | null preserves relevant geometry |
| Predictive baselines | mean-only, focal composition, smoothing, tissue-region covariates | whether added complexity is useful | P0 | grouped held-out performance is measured |
| Multiscale analysis | fixed radii, distance weighting, scale-response curves | at which scale context matters | P0 | scale selection is nested or pre-specified |
| Graph signal analysis | Laplacian energy, graph smoothing, graph regularization | local smoothness and heterogeneity | P1 | not interpreted as biological mechanism |
| Cell composition | marker modules, deconvolution, reference mapping | which cell states or compartments co-localise | P0/P1 | reference and uncertainty are documented |
| Spatial domains | BayesSpace, SpaGCN, STAGATE, Giotto/Seurat domains | whether anatomical/molecular regions exist | P1 | domains are not treated as clones without genomics |
| Cell-cell neighbourhoods | Squidpy, Giotto, neighborhood enrichment, permutation tests | which cell-type pairs co-occur | P1 | composition and tissue abundance are controlled |
| Communication | ligand-receptor scoring, NicheNet-like prioritisation | candidate associations between source and target cells | P2 | reported as hypothesis generation, not causality |
| Deconvolution validation | matched scRNA references, cell2location, stereoscope, RCTD-like methods | whether spot mixtures are plausible | P1 | reference tissue and cross-platform limits are explicit |
| Hierarchical statistics | mixed-effects models, patient random effects, cluster bootstrap | patient-level uncertainty and generalisation | P0 | patient IDs and nested structure are available |
| Confounding control | QC covariates, batch, tissue region, density, depth, sensitivity analyses | whether signal survives technical alternatives | P0 | covariates are measured, not invented |
| Differential spatial analysis | spatially variable genes, region-aware differential expression | where genes or programs vary spatially | P1 | multiple testing and patient grouping are controlled |
| Histology integration | pathology masks, H&E features, morphology embeddings | whether molecular context follows anatomy/morphology | P1/P2 | images and registration are available |
| Image models | Cellpose/StarDist, segmentation, feature extraction | cell-level morphology and density | P2 | segmentation quality is independently checked |
| Multi-omics | matched scRNA, ATAC, CNV, methylation, miRNA | molecular anchoring beyond RNA | P2 | matched data actually exist |
| Genomics | CNV, mutations, clonality, lineage inference | whether transcriptional states map to genomic populations | P2 | DNA/lineage measurements are available |
| Dynamics | longitudinal samples, treatment response, state transitions | temporal or treatment-associated change | P2 | time/order is measured |
| Survival/clinical | survival models, response prediction, risk stratification | clinical utility | P2 | independent clinical endpoints exist |
| Information theory | entropy, mutual information, transfer entropy | dependence, diversity, or information flow | P1/P2 | estimator uncertainty and sample structure are handled |
| Physics-inspired | graph diffusion, reaction-diffusion, energy landscapes | phenomenological spatial dynamics | P2 | equations are identifiable from observed data |
| Game theory | interaction/payoff models, Nash or evolutionary simulations | whether game-like interactions are supported | P3 | independent dynamic or perturbational evidence exists |
| Machine learning | elastic net, random forest, gradient boosting, representation learning | nonlinear prediction and feature ranking | P1 | grouped nested CV beats simple baselines |
| Deep learning | GNNs, CNNs, transformers, multimodal models | high-capacity spatial prediction | P2 | sample size, leakage control, and external validation are adequate |

## Filtering gates

### Gate 1: data availability

Exclude methods requiring measurements absent from the cohort. In particular,
do not infer miRNA, methylation, DNA clones, protein, treatment response, or
time from ordinary spatial RNA counts.

### Gate 2: identifiability

Exclude methods whose parameters cannot be estimated from the available
replication structure. A large number of spots does not solve a shortage of
patients.

### Gate 3: baseline competition

Exclude a complex method from the primary claim if it does not beat a relevant
simple baseline under the same grouped split.

### Gate 4: null and uncertainty

Exclude results without an appropriate null, patient-level uncertainty, or a
clear negative-result rule.

### Gate 5: biological interpretation

Downgrade to exploratory when a method produces association without an
independent biological anchor. Never rename a statistical coefficient as
fitness, payoff, clone, communication, or causal effect.

### Gate 6: reproducibility

Require public code, pinned parameters, provenance, version recording, and
machine-readable output for methods used in the primary result.

## Execution tiers

### Tier P0: necessary now

- data and coordinate audit;
- QC and filtering;
- marker modules adapted to tumour;
- focal versus context baselines;
- spatial autocorrelation;
- multiscale fixed-grid analysis;
- grouped patient-held-out validation;
- patient bootstrap and restricted nulls;
- confounding and platform sensitivity.

### Tier P1: add after P0

- deconvolution or reference mapping;
- spatial domains and neighbourhood enrichment;
- graph signal diagnostics;
- mixed-effects or hierarchical models;
- differential spatial programs;
- benchmark against established spatial methods.

### Tier P2: add only when data support it

- ligand-receptor prioritisation;
- histology and morphology;
- matched single-cell integration;
- CNV/mutation anchoring;
- information-theoretic summaries;
- interpretable machine learning.

### Tier P3: optional and high-risk

- game-theoretic interpretation;
- deep graph neural networks;
- causal inference;
- dynamic models;
- treatment-response or survival models.

## Final reporting filters

The final study must separate:

- verified computation;
- exploratory association;
- negative or inconclusive results;
- hypotheses requiring laboratory data;
- methods not run because prerequisites were absent.

The primary claim should be based on the smallest model that survives grouped
validation. Additional methods can provide triangulation, diagnostics, or
hypothesis generation, but cannot rescue a failed primary endpoint by selective
reporting.
