# InSituGame integrated study copy

This is an independent study copy that combines:

1. the original InSituGame project and its existing results;
2. the public-data metadata audit;
3. the paired E-MTAB-13530 pilot for Patient 10.

The original project at `/home/luca/spatial-game-assay` remains unchanged.
This directory is not a replacement for the original repository and has not
been pushed to GitHub.

## What "merged" means here

The materials are merged into one reproducible study package, but the
biological evidence is **not pooled into one homogeneous estimate**. The
original cohorts and the public pilot differ in cohort design, tissue
selection, depth and provenance. The merged study therefore reports:

- original results as the primary proof-of-concept layer;
- the paired public pilot as an exploratory external evidence layer;
- candidate public cohorts as future transfer/replication layers.

This prevents spot counts from being mistaken for independent patients and
prevents a cross-domain comparison from being reported as a direct
replication.

## Evidence layers

| Layer | Material | Role | Current status |
|---|---|---|---|
| Original study | Existing project data, assay and results | Primary proof-of-concept | Existing, unchanged |
| Paired public pilot | E-MTAB-13530, Patient 10, `P10_B1` and `P10_T1` | Technical and intra-patient exploratory check | Completed |
| External cohorts | GSE300445, GSE176078, E-MTAB-14560, HTAN candidates | Future replication or transfer analysis | Not yet analysed |

## Integrated pilot result

| Sample | Tissue | Regime from current assay | In-sample R² |
|---|---|---|---:|
| `P10_B1` | adjacent | Resistant Dominance | 0.0151 |
| `P10_T1` | tumour | Stag Hunt (Bistable), `x_sep = 0.206` | 0.0201 |

These results are exploratory only. They do not resolve the interpretation of
the coefficients as biological payoffs, establish causal mechanisms, or prove
clinical generalization. The low R² and unequal library depth must be
considered before any biological interpretation.

The first baseline audit found only small in-sample gains for the game design
over composition-plus-neighborhood regression: +0.0012 R² for `P10_B1` and
+0.0034 R² for `P10_T1`. This supports testing the game terms against simpler
models rather than treating their presence as evidence by itself.

Five-fold spatial block validation did not show generalization: game block-CV
R² was 0.0061 for `P10_B1` and -0.0101 for `P10_T1`, while the composition
baseline was 0.0089 and -0.0048. The pilot therefore provides a useful
negative/limiting result: its small in-sample gain does not transfer to held
out spatial blocks.

The extension now includes paired samples from Patients 10, 11 and 15 with
leave-one-patient-out evaluation. Game R² exceeded the composition baseline
for five of six samples, but not for `P10_T1`; several scores remained
negative. This is evidence of heterogeneous transfer, not proof of a
generalizable game mechanism.

The extension was subsequently completed for eight patients with paired
samples available in the same E-MTAB-13530 record. The mean patient-level
game-minus-baseline R² was +0.00185, with 5/8 patients positive. The
10,000-replicate patient bootstrap 95% interval was [-0.00883, +0.01050],
including zero. Thus the full same-study audit does not establish reliable
superiority of the game model.

Marker sensitivity reinforces this limitation. Across four alternative marker
panels and 16 samples, mean held-out game R² remained below the
composition-plus-neighborhood baseline in every panel; the game model won in
4/16, 5/16, 7/16 and 4/16 samples respectively. This is a robustness
failure for an incremental game claim, not evidence that the biological
question is meaningless.

A depth audit found correlations of -0.426 between median UMI and the
game-minus-baseline R² difference and -0.446 between detected features and that
difference. These are unadjusted diagnostics, but they reinforce the need to
model measurement quality before interpreting the interaction terms.

A dependency-free QC-adjusted patient-held-out audit added log median UMI and
log detected features as sample-level covariates. Its mean patient-level
game-minus-baseline R² was +0.00395, with bootstrap 95% CI
[-0.01289, +0.01994]. The interval includes zero; QC adjustment does not
establish incremental biological value.

## Reproducibility

The integrated pilot artifacts are in
[`data/extension_pilot/`](./data/extension_pilot/). The report is in
[`EXTENSION_PILOT_REPORT.md`](./EXTENSION_PILOT_REPORT.md), and the candidate
cohort review is in [`PUBLIC_COHORT_CATALOG.md`](./PUBLIC_COHORT_CATALOG.md).

The original test suite was run in this copy:

```bash
python3 -m unittest discover -s tests -v
```

Result: 19 tests passed.

The interdisciplinary roadmap is documented in
[`INTERDISCIPLINARY_FRAMEWORK.md`](./INTERDISCIPLINARY_FRAMEWORK.md). It
prioritizes measurement validation, spatial/hierarchical statistics and
uncertainty before causal, dynamical or clinical extensions. It explicitly
rejects quantum, relativistic and Millennium Prize analogies as direct
validation methods because they add no testable prediction to these data.

The mathematical applicability audit is documented in
[`MATHEMATICAL_EXTENSION_REPORT.md`](./MATHEMATICAL_EXTENSION_REPORT.md). Its
graph-Laplacian and truncated spectral-zeta calculations are descriptive
spatial diagnostics only; they do not support a Riemann-zeta mechanism or
rescue the uncertain game-model advantage.

The proposed fluid-dynamical extension is assessed in
[`VORTICITY_APPLICABILITY_REPORT.md`](./VORTICITY_APPLICABILITY_REPORT.md).
No velocity field or time series is available, so Navier–Stokes vorticity
cannot be estimated from the current snapshot without inventing a physical
quantity.

The most suitable additional analysis currently supported by the data is the
multiscale spatial audit in
[`MULTISCALE_SPATIAL_AUDIT.md`](./MULTISCALE_SPATIAL_AUDIT.md). It finds an
intermediate-scale OOF improvement in the original lung and breast sections,
but the scale was selected from the same data and is not independent patient
validation.

The pre-specified 150--300 µm window was subsequently tested with
patient-held-out training across eight public patients; results are in
[`MULTISCALE_PATIENT_HOLDOUT_REPORT.md`](./MULTISCALE_PATIENT_HOLDOUT_REPORT.md).
Mean game-minus-baseline R² was positive at all four radii
(+0.00185 to +0.00531), but every patient bootstrap interval included zero.
This is a useful scale-localized exploratory result, not confirmatory evidence
of a game mechanism.

The contribution and novelty boundary is documented in
[`CONTRIBUTION_AND_NOVELTY_SCOPE.md`](./CONTRIBUTION_AND_NOVELTY_SCOPE.md).
The defensible contribution is the interpretable validation protocol that
separates spot composition, multiscale neighbourhood context and bilinear
interaction terms on human tumour spatial transcriptomics. The individual
ingredients are not claimed as new, and the current empirical result does not
validate a game-specific effect.

Before adding microRNA or methylome claims, the input readiness audit in
[`MOLECULAR_ANCHOR_READINESS.md`](./MOLECULAR_ANCHOR_READINESS.md) must be
respected: the current cache has restricted gene-level spatial expression but
no microRNA counts, CpG matrix, or matched PDAC/biliary cohorts.

The primary defensible specification is now summarized in
[`CORE_VALIDATED_MODEL.md`](./CORE_VALIDATED_MODEL.md). It defines the
focal-composition and neighbourhood baselines, the bilinear interaction as an
exploratory extension, and the gates required before promoting it to a
biological claim.

The fixed 150 µm non-game audit is reported in
[`CORE_SPATIAL_CONTEXT_AUDIT.md`](./CORE_SPATIAL_CONTEXT_AUDIT.md). It improves
spatial-block OOF R² over focal composition in all three cached cohorts
(Breast +0.0498, Melanoma +0.0787, Lung +0.0249). This supports spatial
context as the primary descriptive result, without implying causality or
validating the game interpretation.

The same core was tested patient-held-out at 150 µm across eight paired public
patients; see [`CORE_PATIENT_HOLDOUT_REPORT.md`](./CORE_PATIENT_HOLDOUT_REPORT.md)
and [`data/extension_pilot/core_patient_holdout.json`](./data/extension_pilot/core_patient_holdout.json).
The mean patient-level increment was only +0.00091 R² (6/8 patients positive),
with bootstrap 95% CI [-0.00163, +0.00316], so the patient-held-out result
remains uncertain and does not confirm a reproducible cross-patient gain.

An independent breast cohort, [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078),
was audited and tested using the public processed matrices; see
[`data/public_validation/GSE176078_VALIDATION.md`](./data/public_validation/GSE176078_VALIDATION.md).
Its six-patient held-out result was negative/uncertain (mean increment
-0.06454 R²; 3/6 positive; bootstrap 95% CI [-0.39027, +0.18983]).

The same skeleton was then applied with a lung-specific marker representation
to [GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487);
see [`data/public_validation/GSE189487_VALIDATION.md`](./data/public_validation/GSE189487_VALIDATION.md).
The six-sample held-out result was also negative/uncertain (mean increment
-0.04793 R²; 4/6 positive; bootstrap 95% CI [-0.29986, +0.14579]).

The project architecture has been broadened beyond game theory in
[`EXPANDED_METHOD_ARCHITECTURE.md`](./EXPANDED_METHOD_ARCHITECTURE.md).
Spatial statistics, graph signal models, multiscale analysis, composition-aware
regression, hierarchical validation, pathology integration and exploratory
communication analysis are now first-class alternatives. The game-theoretic
module remains optional and requires independent evidence for dynamics or
perturbation before its coefficients can be interpreted as biological payoffs.

The complete method inventory and filtering policy are in
[`METHOD_INVENTORY_AND_FILTER.md`](./METHOD_INVENTORY_AND_FILTER.md). It
includes bioinformatics, spatial statistics, graph methods, deconvolution,
single-cell integration, imaging, multi-omics, information theory, machine
learning and clinical modelling, but assigns them tiers and hard prerequisites
so that only necessary, identifiable and validated methods enter the primary
claim.

The October 2026 rigorous protocol, including boolean QC filters, matrix
orientation and z-score conventions, patient-level inference, Bayesian
sensitivity summaries, FDR control, negative controls and AI-tool gates is in
[`RESEARCH_PROTOCOL_2026.md`](./RESEARCH_PROTOCOL_2026.md).
The current literature and method-status register is in
[`CURRENT_METHODS_EVIDENCE_2026.md`](./CURRENT_METHODS_EVIDENCE_2026.md).

For the original lung/osimertinib direction, the GSE189487 cohort has been
audited separately in
[`data/public_validation/GSE189487_OSIMERTINIB_RELEVANCE.md`](./data/public_validation/GSE189487_OSIMERTINIB_RELEVANCE.md).
It is an untreated LUAD progression cohort, so the analysis can test only
candidate resistance-associated transcriptional programs and their spatial
organisation, not osimertinib response or resistance itself.
The lung-specific multiscale sensitivity grid is reported in
[`data/public_validation/GSE189487_MULTISCALE.md`](./data/public_validation/GSE189487_MULTISCALE.md).

The focused EGFR-mutant NSCLC Visium HD cohort GSE301973 has now passed a
complete integrity and patient mapping verification; see
[`data/public_validation/GSE301973_SMALL_AUDIT.md`](./data/public_validation/GSE301973_SMALL_AUDIT.md) and
[`data/public_validation/GSE301973_MAPPING_AUDIT.md`](./data/public_validation/GSE301973_MAPPING_AUDIT.md).
It contains four pretreatment patients with PFS metadata and cell-level
annotations (Slide A: `s1_PFS_69` and `s2_PFS_22`, long PFS; Slide C: `s3_PFS_15` and `s4_PFS_12`, short PFS).
The sample positions documented in `biological_sample_metadata.xlsx` (`Left`/`Right`) were verified against
the bimodal coordinate distributions of the two slides.
The patient-level progression-free survival (PFS) spatial analysis is implemented in
[`src/analyze_gse301973_patient_pfs.py`](./src/analyze_gse301973_patient_pfs.py) and reported in
[`data/public_validation/GSE301973_PFS_ANALYSIS.md`](./data/public_validation/GSE301973_PFS_ANALYSIS.md).
Key findings reveal that short-PFS tumors exhibit markedly higher cellular confluence and homotypic
tumor-tumor contact (mean 75.9% vs 64.3%), driven by patient `s3_PFS_15` (77.5% tumor purity, 87.9% homotypic contact),
while long-PFS tumors exhibit richer heterotypic integration with the surrounding stroma and vasculature before treatment.

## Required interpretation

The merged package may support a manuscript section titled, for example,
"Public-data paired exploratory analysis". It must not yet be described as a
multi-patient validation, because the analysed public pilot contains one
patient and two paired samples.
