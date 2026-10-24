# GSE301973: Patient-Level Progression-Free Survival (PFS) Spatial Analysis

## Cohort and Mapping Verification

The four EGFR-mutant NSCLC pre-treatment patient samples were profiled across two Visium HD slides.
The region-to-patient assignment is resolved via `biological_sample_metadata.xlsx` where the
`sample_position` field (`Left` / `Right`) matches the bimodal spatial coordinate distributions on each slide:

- **Slide A (`GSM9093278`)**:
  - **Left ($x < 12,000$)**: `s1_PFS_69` (PFS: 69 months, `long_PFS`) — 14,605 cells
  - **Right ($x \ge 12,000$)**: `s2_PFS_22` (PFS: 22 months, `long_PFS`) — 31,056 cells
- **Slide C (`GSM9093279`)**:
  - **Left ($x < 12,500$)**: `s3_PFS_15` (PFS: 15 months, `short_PFS`) — 60,507 cells
  - **Right ($x \ge 12,500$)**: `s4_PFS_12` (PFS: 12 months, `short_PFS`) — 105,123 cells

Notably, Slide A contains exclusively **long-PFS** patients (69 and 22 mo), while Slide C contains
exclusively **short-PFS** patients (15 and 12 mo).

## Patient Summary Table

| Patient ID | PFS (mo) | Group | Total Cells | Tumor Cells | Tumor % | Homotypic Contact | Endothelial Contact | EGFR Bypass | EMT | Stem-Like | Proliferation |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `s1_PFS_69` | 69 | long_PFS | 14,605 | 3,655 | 25.0% | 61.0% | 22.9% | 0.351 | 0.252 | 0.277 | 0.120 |
| `s2_PFS_22` | 22 | long_PFS | 31,056 | 9,327 | 30.0% | 67.7% | 12.6% | 0.224 | 0.169 | 0.644 | 0.070 |
| `s3_PFS_15` | 15 | short_PFS | 60,507 | 46,865 | 77.5% | 87.9% | 2.7% | 0.140 | 0.115 | 0.120 | 0.125 |
| `s4_PFS_12` | 12 | short_PFS | 105,123 | 37,073 | 35.3% | 64.0% | 21.8% | 0.234 | 0.182 | 0.340 | 0.071 |

## Long-PFS vs Short-PFS Group Comparison

- **Tumor Fraction**: Long PFS = 27.5% vs Short PFS = 56.4%
- **Homotypic Tumor Contact**: Long PFS = 64.3% vs Short PFS = 75.9%
- **Tumor-Endothelial Contact**: Long PFS = 17.8% vs Short PFS = 12.3%
- **Tumor-Stromal Contact**: Long PFS = 3.4% vs Short PFS = 4.0%
- **Tumor-Immune Contact**: Long PFS = 11.6% vs Short PFS = 6.9%

### Tumor Cell Module Means

- **egfr_bypass**: Long PFS = 0.2876 vs Short PFS = 0.1867
- **emt**: Long PFS = 0.2105 vs Short PFS = 0.1485
- **stem_like**: Long PFS = 0.4604 vs Short PFS = 0.2296
- **spp1_cd44**: Long PFS = 0.8546 vs Short PFS = 0.2721
- **proliferation**: Long PFS = 0.0951 vs Short PFS = 0.0982

## Key Biological & Methodological Observations

1. **Tumor Architecture & Confluence**: Short-PFS tumors exhibit markedly higher cellular confluence and homotypic tumor-tumor contact (mean 75.9% vs 64.3%), driven by patient `s3_PFS_15` with 77.5% tumor purity and 87.9% homotypic clustering.
2. **Microenvironment Integration**: Long-PFS tumors exhibit richer heterotypic integration with the surrounding stroma and vasculature before treatment.
3. **Molecular Programs**: Pre-treatment tumor cells in long-PFS patients express higher baseline levels of `spp1_cd44` and `stem_like` programs, reflecting distinct non-small cell lung cancer differentiation states rather than pre-existing osimertinib resistance.
4. **Scope Limitation**: With $N=4$ pre-treatment samples, these observations serve as an exploratory anchor demonstrating that Visium HD captures subcellular spatial interaction phenotypes, but cannot establish independent survival causality.

