# Public validation cohorts: first priority

This is the first-pass registry of public cohorts selected for validation of
the minimal spatial-context assay. The selection prioritizes downloadable
expression/spatial data, independent donors or patients, and metadata that can
support patient-held-out testing.

## Priority 1: start here

| Priority | Accession | Disease | Platform | Public material | Why it is first | Main caveat |
|---|---|---|---|---|---|---|
| 1 | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) / [Zenodo mirror](https://zenodo.org/records/4739739) | Breast carcinoma | 10x Visium | Raw-count and filtered matrices, spatial files, images and annotations are reported in the public package | Six primary cancers/donors; same broad platform as the current pilot; suitable for patient-held-out replication | Six patients is still a small validation cohort; subtypes and tissue composition vary |
| 2 | [GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487) | Lung adenocarcinoma progression | 10x Visium | Public processed spatial data and study metadata | Spatial samples span AIS, MIA and invasive disease, giving a biologically ordered test | Stage and patient can be partly confounded; patient mapping must be checked before analysis |
| 3 | [GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318) / [10x code](https://github.com/10XGenomics/HumanColonCancer_VisiumHD) | Colorectal cancer | Visium HD, FFPE | Public count/processed data and analysis code are reported | Largest and highest-resolution candidate in this first wave; useful for fine-scale neighbourhood tests | Reported sample count is not automatically patient count; donor/section mapping is mandatory |

## Priority 2: useful independent or complementary cohorts

| Accession/source | Disease | Platform | Intended use | Limitation |
|---|---|---|---|---|
| [GSE246518](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246518) | BRCA2-associated breast tissue and controls | 10x Visium | Genotype/risk-group and breast-context replication | Not a pure tumour-only cohort; privacy restrictions affect raw-data availability |
| [GSE226997](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226997) | Colorectal cancer | 10x Visium | Tumour–immune–stromal neighbourhood transfer test | Four patients; section and donor independence must be verified |
| [GSE108290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108290) | Cutaneous malignant melanoma | Original spatial transcriptomics array | Historical melanoma replication and qualitative spatial organisation | Only two patients and pre-Visium technology; not sufficient for a definitive melanoma claim |
| [Mendeley 29ntw7sh4r](https://www.spatialresearch.org/resources-published-datasets/10-1038-s41551-020-0578-x/) | Breast cancer | Earlier spatial transcriptomics | Morphology-oriented complementary replication | Older platform and donor count require direct metadata audit |

## Discovery sources, not yet validation cohorts

- [HTAN](https://humantumoratlas.org/) and its [data updates](https://humantumoratlas.org/data-updates):
  rich cancer atlases and multiple assay types, but each release needs a
  study-level audit of download status, donor identifiers, assay compatibility,
  and licensing.
- [CROST](https://ngdc.cncb.ac.cn/crost/), [STOmicsDB](https://db.cngb.org/stomics/)
  and [OncoMap](https://www.oncomap.org/): useful catalogues for finding
  additional cohorts, not automatically harmonized datasets.
- [10x Genomics datasets](https://www.10xgenomics.com/datasets): excellent
  technical controls, but many entries are single-donor demonstrations and
  therefore weak for patient-level generalization.
- [HEST-1k](https://arxiv.org/html/2406.16192v2): large histology/spatial
  benchmark resource; suitability depends on whether the specific cancer
  subset has expression, coordinates, donor metadata and compatible licensing.

## Validation rules

1. The independent unit is the patient/donor, not the spot.
2. Sections from one donor are grouped in the same held-out fold.
3. A cohort is not called independent until accession, donor mapping, platform,
   and preprocessing are verified from the repository files.
4. Processed matrices can support this assay; FASTQ availability is not
   required, but raw-count provenance must be recorded.
5. A positive result is not expected in advance. Effect size, patient-level
   uncertainty, and null-model behaviour will be reported even if the result
   is negative.

## What is confirmed versus pending

**Confirmed at the catalogue level:** these accessions and source pages exist,
and the listed study/platform descriptions are publicly reported.

**Pending before analysis:** exact unique-patient counts, sample-to-donor
mapping, marker availability, coordinate units, tissue filtering, and whether
the downloadable files contain the required count and spatial-coordinate
layers.

**First completed cohort:** GSE176078 passed the file/metadata audit and was
then evaluated with patient-held-out validation. The mean increment was
-0.06454 R² (3/6 patients positive; bootstrap 95% CI [-0.39027, +0.18983]),
so it is currently a negative/uncertain replication rather than a confirmation.

**Second completed cohort:** GSE189487 was evaluated with a lung-specific
marker representation and the same universal skeleton. The mean increment was
-0.04793 R² (4/6 samples positive; bootstrap 95% CI [-0.29986, +0.14579]),
also negative/uncertain at this sample size.

**Third completed cohort:** GSE301973 (Visium HD NSCLC pre-treatment osimertinib cohort)
passed full mapping verification via `biological_sample_metadata.xlsx` (Left/Right sample
positions matching coordinate bimodal distributions). Patient-level analysis across all
4 patients (`s1_PFS_69`, `s2_PFS_22`, `s3_PFS_15`, `s4_PFS_12`) demonstrated that short-PFS
tumors exhibit higher cellular confluence and homotypic contact (75.9% vs 64.3%), while
long-PFS tumors exhibit richer heterotypic integration with the vascular/stromal stroma.
Reported in `data/public_validation/GSE301973_PFS_ANALYSIS.md`.

**Not inferable from public metadata alone:** causal interaction, evolutionary
payoff, clinical utility, miRNA or methylome mechanisms.
