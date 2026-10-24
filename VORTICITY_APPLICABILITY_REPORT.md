# Vorticity and Navier–Stokes applicability audit

## Question evaluated

The proposed extension was the use of vortex principles from fluid dynamics,
possibly motivated by recent claims concerning the three-dimensional
Navier–Stokes Millennium Prize problem.

The official Clay Mathematics Institute source should be treated as the
authority for the status of that problem. A public claim or an AI-generated
summary is not equivalent to a peer-reviewed proof or a Clay-accepted
solution. The claim therefore cannot be imported as a solved theorem into
this study.

## Why physical vorticity cannot be calculated here

For an incompressible fluid, vorticity is a vector field

`omega = curl(u)`,

where `u(x,t)` is a measured or modelled velocity field. The Navier–Stokes
vorticity equation also requires time, viscosity, forcing and boundary
conditions. The current spatial transcriptomics data contain:

- a two-dimensional coordinate per spot;
- molecular counts and derived scalar marker scores;
- no cell displacement or velocity;
- no pressure, density, viscosity or time-resolved trajectory;
- no experimentally measured flow of cells or extracellular fluid.

Consequently, a curl computed from marker gradients would not be
fluid-mechanical vorticity. It would be a derived image-processing quantity
whose interpretation depends entirely on an arbitrary construction.

## What was tested instead

The merged study now explicitly classifies this proposal:

| Proposed tool | Current status | Reason |
|---|---|---|
| 3-D Navier–Stokes regularity theorem | Not imported | It is a theorem-status question, not an estimator for transcript counts |
| Physical vorticity `curl(u)` | Not identifiable | No velocity field or time series |
| Scalar marker-gradient curl | Rejected as biology claim | It would be a constructed pseudo-field, not fluid vorticity |
| Graph Laplacian/spectral diagnostics | Applicable | They use the measured coordinates and are documented separately |
| Reaction-diffusion model | Future, conditional | Requires time series or perturbation data and calibrated rates |
| Cell migration/flow model | Future, conditional | Requires live imaging, lineage tracking or repeated spatial samples |

## Falsifiable route if appropriate data become available

With live imaging or registered longitudinal sections, one could pre-specify:

1. a cell-displacement field `u(x,t)`;
2. a registration and uncertainty model;
3. a boundary and smoothing rule;
4. `curl(u)` and strain summaries;
5. a null model based on isotropic random motion;
6. prediction of future phenotype transport or treatment response.

Only then would a vortex-like hypothesis produce a biological test. A
single static transcriptomic section cannot distinguish a vortex, diffusion,
active migration, tissue deformation, or a static spatial gradient.

This negative result is scientifically useful: it prevents a mathematically
prestigious but non-identifiable object from being presented as evidence for
the evolutionary game model.
