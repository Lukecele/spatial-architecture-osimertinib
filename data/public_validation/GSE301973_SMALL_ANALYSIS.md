# GSE301973 small spatial analysis

This is the first descriptive analysis after the integrity audit. It does not
test PFS or claim osimertinib resistance.

## Design

- two candidate coordinate-derived regions per Visium HD slide;
- cell-type composition from the supplied StarDist-integrated annotations;
- nearest-neighbour composition;
- module scores for:
  - EGFR/bypass signalling;
  - EMT;
  - stem-like state;
  - SPP1-CD44;
  - proliferation.

The code is [`src/analyze_gse301973_small.py`](../../src/analyze_gse301973_small.py)
and the output is [`GSE301973_small_analysis.json`](./GSE301973_small_analysis.json).

## Critical mapping rule

The integrated objects identify the slide and contain coordinates, but do not
directly identify each cell with the biological patient ID. Therefore the two
coordinate-derived regions are explicitly labelled `candidate_region`.
They are not assigned to the long-PFS or short-PFS metadata groups.

Until the left/right or sample-position mapping is verified, any comparison
between these regions and PFS would be invalid. This analysis is therefore
descriptive and suitable for checking whether the spatial programs are
measurable, not for clinical inference.

The official 2 µm tissue-position parquet files were also inspected. They
contain barcode, tissue flag, array coordinates and image coordinates, but no
patient identifier or left/right sample label. The integrated objects likewise
identify only slide A or C. Therefore the region-to-patient assignment remains
unverified from the downloaded public files.

## Descriptive observation

The candidate regions differ in cell composition and in EGFR/bypass, EMT,
stem-like, SPP1-CD44 and proliferation module means. This confirms that the
programs are measurable, but it does **not** distinguish a patient/PFS effect
from anatomical composition or tissue-region structure.

## Expansion gate

The study will expand only after:

1. sample-to-patient spatial mapping is verified;
2. per-patient QC is reported;
3. region-level estimates are stable to the partition rule;
4. leave-one-patient-out sensitivity is feasible;
5. the resulting interpretation remains explicitly exploratory at n=4.
