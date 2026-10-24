# Molecular-anchor readiness audit

The proposed integration of Carotenuto-related microRNA and methylation
evidence starts with an input audit rather than an inferred biological claim.
The implementation is
[`src/molecular_anchor_audit.py`](./src/molecular_anchor_audit.py), and its
output is
[`results/MOLECULAR_ANCHOR_AUDIT.json`](./results/MOLECULAR_ANCHOR_AUDIT.json).

## What the current cache contains

The bundled data contain gene-level spatial expression for a restricted marker
panel and coordinates. They do **not** contain:

- microRNA counts;
- methylation or CpG beta values;
- matched PDAC or biliary tumour cohorts;
- perturbation labels for miR-1307, miR-181, or MITF/APAF-1.

## Anchor status

| Published biological anchor | Current data status | Permitted use |
|---|---|---|
| MITF/APAF-1 resistance axis in melanoma | MITF present, APAF1 absent | Exploratory MITF-associated state only; not pathway validation |
| miR-1307/CLIC5 and FOLFIRINOX response | Neither measured; no PDAC cohort | Cannot be analysed from current cache |
| miR-181/SIRT1 drug response in biliary cancer | Neither measured; no biliary cohort | Cannot be analysed from current cache |

The current melanoma panel contains MITF and other melanocytic/resistance
descriptors (MLANA, PMEL, TYR, WNT5A, VIM, AXL), but that does not reproduce
the published MITF/APAF-1 mechanism. It is therefore not valid to call the
result a microRNA or methylome validation.

## What can be done now

The defensible immediate analysis is an exploratory spatial map of the
available melanoma transcriptomic state, with the existing composition,
neighbourhood, multiscale and patient-held-out audits. It must be labelled as
a transcriptomic state analysis.

## What requires new data

To test the published anchors mechanistically, obtain matched data with:

1. small-RNA sequencing or spatial microRNA measurements;
2. methylation/CpG measurements if an epigenetic claim is intended;
3. the relevant tumour type and treatment context;
4. perturbation or independent response labels;
5. registered spatial coordinates where spatial localization is claimed.

This audit is a productive first step because it prevents overclaiming and
identifies the exact missing measurements before a new analysis is designed.
