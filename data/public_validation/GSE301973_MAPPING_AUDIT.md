# GSE301973 region-to-patient mapping audit

## Question

Can the public files assign each spatial region on slide A or C to one of the
four biological patients and therefore to a PFS group?

## Sources checked

- GEO series record [GSE301973](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301973);
- sample records `GSM9093278` and `GSM9093279`;
- public biological sample metadata workbook (`biological_sample_metadata.xlsx`);
- 2 µm tissue-position parquet files;
- StarDist-integrated H5AD objects;
- related BioProject [PRJNA1218489](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1218489);
- related SRA run metadata.

## Verified and Resolved

1. **Workbook metadata**: `biological_sample_metadata.xlsx` explicitly provides the four-patient mapping:
   - `s1_PFS_69`: Slide A, `GSM9093278`, `VisiumHD_1`, Pretreatment, PFS 69 mo, `long_PFS`, `sample_position: Left`
   - `s2_PFS_22`: Slide A, `GSM9093278`, `VisiumHD_1`, Pretreatment, PFS 22 mo, `long_PFS`, `sample_position: Right`
   - `s3_PFS_15`: Slide C, `GSM9093279`, `VisiumHD_2`, Pretreatment, PFS 15 mo, `short_PFS`, `sample_position: Left`
   - `s4_PFS_12`: Slide C, `GSM9093279`, `VisiumHD_2`, Pretreatment, PFS 12 mo, `short_PFS`, `sample_position: Right`
2. **Coordinate verification**:
   - `spatial_x` coordinates correlate perfectly ($r = 1.0000$) with `array_col`, establishing that low $x$ corresponds to `Left` and high $x$ corresponds to `Right`.
   - **Slide A**: Has a prominent physical coordinate gap at $x \approx 11,800 - 12,200$, separating the Left specimen ($x < 12,000$, $N=14,605$, `s1_PFS_69`) from the Right specimen ($x \ge 12,000$, $N=31,056$, `s2_PFS_22`).
   - **Slide C**: Separates into two distinct histological specimens at $x \approx 12,500$: the Left specimen ($x < 12,500$, $N=60,507$, `s3_PFS_15`) is a dense confluent tumor (77.5% tumor cells, 87.9% homotypic contact), while the Right specimen ($x \ge 12,500$, $N=105,123$, `s4_PFS_12`) is a mixed microenvironment (35.3% tumor cells, 28.4% endothelial cells).
3. **Cohort level observation**:
   - Slide A contains exclusively long-PFS patients (69 and 22 months).
   - Slide C contains exclusively short-PFS patients (15 and 12 months).

## Decision

The region-to-patient mapping is **verified from the public files and metadata**.
The patient-level spatial and module comparison is implemented in:
- [`src/analyze_gse301973_patient_pfs.py`](../../src/analyze_gse301973_patient_pfs.py)
- [`GSE301973_PFS_ANALYSIS.md`](./GSE301973_PFS_ANALYSIS.md)
- [`GSE301973_patient_pfs_analysis.json`](./GSE301973_patient_pfs_analysis.json)

## Osimertinib boundary and scope

These four samples are pre-treatment clinical biopsies with known progression-free
survival under subsequent osimertinib therapy. While they provide an authentic clinical
molecular anchor for Visium HD spatial analysis, $N=4$ remains hypothesis-generating
and exploratory. It demonstrates that spatial tumor architecture, cellular confluence,
and microenvironmental contact are quantifiable at subcellular resolution in clinical
specimens, but does not claim predictive clinical validation.
