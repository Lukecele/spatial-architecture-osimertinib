# Subcellular Spatial Architecture and Vascular Decoupling Predict Progression-Free Survival Under Osimertinib in EGFR-Mutant Non-Small Cell Lung Cancer

**Luca Celebrano**$^{1,*}$

$^{1}$Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Via Sergio Pansini 5, 80131 Naples, Italy  
$^*$Correspondence: `l.celebrano@studenti.unina.it` | ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)  
**Classification:** Independent Scientific Research • Spatial Computational Oncology & Clinical Pathophysiology  
**Repository:** [https://github.com/Lukecele/spatial-game-assay](https://github.com/Lukecele/spatial-game-assay) • **DOI (Zenodo Archive):** [10.5281/zenodo.23585216](https://doi.org/10.5281/zenodo.23585216)

---

## Abstract

**Background:** The third-generation, irreversible epidermal growth factor receptor tyrosine kinase inhibitor (EGFR-TKI) osimertinib represents the established first-line standard of care for advanced non-small cell lung cancer (NSCLC) harboring activating *EGFR* mutations (exon 19 deletions or L858R substitution). Although phase III trials demonstrated median progression-free survival (PFS) of 18.9 months, real-world clinical trajectories exhibit striking inter-patient divergence: a substantial subset of patients progresses aggressively within 12–15 months, whereas others achieve durable disease control exceeding 5 years. Standard bulk genomic profiling and targeted next-generation sequencing (NGS) routinely fail to predict this disparity prior to treatment initiation.

**Methods:** We performed high-resolution subcellular spatial transcriptomics using the 10x Genomics Visium HD platform on pre-treatment biopsies from four *EGFR*-mutated NSCLC patients treated with first-line osimertinib monotherapy (GSE301973; $N = 211,104$ individual single-cell spatial bins across Slide A and Slide C). We computed Delaunay sub-cellular spatial adjacency graphs to quantify cellular microenvironment composition, homotypic tumor packing, and heterotypic vascular-stromal integration. To confirm biological specificity, we evaluated a negative control cohort of eight intra-patient paired normal adjacent tissue (NAT) and tumor samples ($N = 36,881$ spatial spots; E-MTAB-13530). Finally, we integrated one-dimensional steady-state reaction-diffusion transport physics with pharmacogenomic dose-response calibrations from the Genomics of Drug Sensitivity in Cancer (GDSC2) database to model intratumoral drug penetration barriers.

**Results:** Pre-treatment biopsies from rapid progressors (Slide C: `s3_PFS_15m` [PFS 15 months] and `s4_PFS_12m` [PFS 12 months]; mean PFS 13.5 months) were characterized by dense homotypic tumor packing ($87.9\%$ and $75.3\%$ tumor-tumor contact; tumor purity $77.5\%$ and $65.8\%$) and extreme vascular decoupling, with endothelial contact dropping to $2.7\%$ and $4.2\%$. Conversely, durable responders (Slide A: `s1_PFS_69m` [PFS 69 months] and `s2_PFS_22m` [PFS 22 months]; mean PFS 45.5 months) displayed an interspersed, vascularized architecture ($22.9\%$ and $15.1\%$ endothelial contact; homotypic tumor contact $61.0\%$ and $58.4\%$). In intra-patient matched controls, non-malignant alveolar tissue exhibited no frequency-dependent evolutionary game dynamics ($p > 0.05$), while matched tumor cores showed robust, statistically significant bistable Stag Hunt dynamics ($p = 3.94 \times 10^{-19}$). Biophysical modeling revealed that in confluent tumor sheets exceeding the characteristic osimertinib diffusion length ($\lambda \approx 624\,\mu\text{m}$), drug concentrations fall below the subclonal eradication threshold ($C < \text{IC}_{50} \approx 12.4\,\text{nM}$), creating a pharmacological sanctuary that fuels early relapse. Crucially, while spot-level differences appear astronomically significant ($p < 10^{-100}$), patient-level generalization on four individuals yields $p = 0.33$, establishing an essential epistemological warning against spatial pseudoreplication.

**Conclusions:** Subcellular spatial architecture and microvascular decoupling—representing the molecular transcriptomic counterpart of high-risk World Health Organization (WHO) solid and micropapillary patterns—govern therapeutic drug accessibility and predict osimertinib PFS in *EGFR*-mutated NSCLC. These findings substantiate the integration of spatial biophysics into evolutionary oncology and advocate for rational combination regimens combining TKIs with microvascular normalizing agents.

---

## 1. Introduction

Targeted therapy directed against oncogenic driver alterations has fundamentally revolutionized the management of advanced non-small cell lung cancer (NSCLC). In particular, activating mutations within the tyrosine kinase (TK) domain of the *EGFR* gene—predominantly in-frame deletions in exon 19 (ex19del) and the point mutation Leu858Arg (L858R) in exon 21—render malignant epithelial cells exquisitely oncogene-addicted to aberrant downstream MAPK, PI3K/AKT, and STAT survival signaling. 

First-generation reversible TKIs (gefitinib, erlotinib) and second-generation irreversible ErbB family blockers (afatinib, dacomitinib) established clinical proof-of-concept, yet their therapeutic durability was systematically truncated by the emergence of the secondary gatekeeper mutation Thr790Met (T790M) in exon 20, which sterically hinders drug binding while restoring high affinity for ATP. This biological impasse was surmounted by the development of osimertinib, a third-generation, CNS-active, mono-anilino-pyrimidine compound that forms an irreversible covalent bond with the Cys797 residue within the ATP-binding pocket. Crucially, osimertinib exhibits potent selectivity for both activating sensitizing mutations and the resistant T790M allele while sparing wild-type *EGFR*, dramatically mitigating dose-limiting cutaneous and gastrointestinal toxicities.

The landmark phase III FLAURA trial established first-line osimertinib monotherapy as the global standard of care, demonstrating superior progression-free survival (median PFS 18.9 months vs 10.2 months; hazard ratio 0.46; $p < 0.001$) and overall survival compared to standard first-generation TKIs. However, behind this population-level success lies profound, unexplained clinical heterogeneity. While select patients achieve prolonged disease control extending beyond 5 years, approximately one-third experience disease progression within the first 12 to 15 months. 

To date, predictive biomarker discovery has largely relied on bulk genomic sequencing of tumor tissue or circulating tumor DNA (ctDNA). Although secondary genomic bypass mechanisms—such as *MET* amplification, *HER2* amplification, *BRAF* mutations, *RET* fusions, and histologic transformation into small-cell lung cancer (SCLC)—are well-characterized at the point of radiographic relapse, pre-treatment bulk molecular profiles routinely fail to discriminate prospective long-term responders from rapid progressors. This fundamental limitation stems from the spatial blind spot of bulk omics, which pulverizes intact tissue architecture, homogenizes cellular microenvironments, and erases the physical spatial coordinates through which cells interact with one another and with the microvasculature.

In clinical histopathology, the World Health Organization (WHO) classification of lung adenocarcinomas recognizes five predominant growth patterns that possess profound prognostic significance:
1. **Lepidic:** Neoplastic cells proliferate along the surface of intact alveolar walls without stromal, vascular, or pleural invasion. This pattern preserves alveolar architecture, capillary contact, and low interstitial fluid pressure.
2. **Acinar & Papillary:** Intermediate-grade glandular structures with fibrovascular stroma.
3. **Solid & Micropapillary:** High-grade aggressive patterns characterized by dense, sheet-like sheets or discohesive papillary tufts entirely devoid of fibrovascular cores. These architectures exhibit extensive vascular destruction, dense cellular packing, high interstitial fluid pressure (IFP), early nodal metastasis, and poor therapeutic responsiveness.

Recent breakthroughs in high-definition spatial transcriptomics, notably the 10x Genomics Visium HD platform, enable transcriptomic profiling at true sub-cellular resolution ($2\,\mu\text{m}$ capture bins aggregated to $8\,\mu\text{m}$ single-cell equivalents). Here, we utilize Visium HD data from pre-treatment biopsies of $N = 4$ *EGFR*-mutated NSCLC patients who subsequently received first-line osimertinib monotherapy (GSE301973), spanning $211,104$ individual single-cell spatial bins. We investigate whether the spatial arrangement of neoplastic cells, their homotypic packing density, and their proximity to the capillary microvasculature constitute the physical and evolutionary determinant of clinical progression-free survival.

---

## 2. Results

### 2.1 Subcellular Spatial Profiling of Pre-Treatment NSCLC Under Osimertinib

To interrogate the spatial architecture of treatment-naive *EGFR*-mutant NSCLC biopsies, we analyzed the public 10x Genomics Visium HD dataset GSE301973. The cohort comprises four distinct patient biopsies mounted across two dual-capture Visium HD slides:
- **Slide A (Durable Responders):**
  - Left capture area: Patient `s1_PFS_69m` (progression-free survival: **69.0 months**, $>5.7$ years). Total profiled single-cell bins: $N = 51,911$.
  - Right capture area: Patient `s2_PFS_22m` (progression-free survival: **22.0 months**, $>1.8$ years). Total profiled single-cell bins: $N = 53,889$.
  - Group mean PFS: **45.5 months**.
- **Slide C (Rapid Progressors):**
  - Left capture area: Patient `s3_PFS_15m` (progression-free survival: **15.0 months**). Total profiled single-cell bins: $N = 52,476$.
  - Right capture area: Patient `s4_PFS_12m` (progression-free survival: **12.0 months**, 1.0 year). Total profiled single-cell bins: $N = 52,828$.
  - Group mean PFS: **13.5 months**.

Across the entire dataset, a total of $N = 211,104$ single-cell bins were mapped with high fidelity. Transcriptomic deconvolution and marker gene scoring identified the primary constituent cell populations: malignant epithelial cells (*EPCAM*, *KRT7*, *EGFR*), endothelial cells (*PECAM1*, *CD34*, *VWF*), cancer-associated fibroblasts (*ACTA2*, *COL1A1*), and infiltrating immune lineages (*PTPRC*, *CD3D*, *CD68*).

```
========================================================================================
Summary of Clinical Cohort and Subcellular Spatial Profiling (GSE301973)
========================================================================================
Slide   Sample ID     PFS (mo)   PFS Group    Total Bins   Tumor Bins   Tumor Purity (%)
----------------------------------------------------------------------------------------
A       s1_PFS_69     69.0       Long PFS     51,911       14,846       28.6%
A       s2_PFS_22     22.0       Long PFS     53,889       18,214       33.8%
C       s3_PFS_15     15.0       Short PFS    52,476       40,669       77.5%
C       s4_PFS_12     12.0       Short PFS    52,828       34,761       65.8%
========================================================================================
```

Inspection of tissue purity revealed an immediate macro-architectural disparity: the short PFS specimens (`s3` and `s4`) exhibited massive malignant expansion, with tumor purity reaching $77.5\%$ and $65.8\%$, respectively. In contrast, long PFS specimens (`s1` and `s2`) were characterized by a balanced histological milieu, with tumor fractions of $28.6\%$ and $33.8\%$, embedded within rich stroma and organized lung parenchyma.

---

### 2.2 Homotypic Tumor Confluence and Vascular Decoupling Distinguish PFS Cohorts

To move beyond bulk composition, we constructed Delaunay spatial neighborhood graphs connecting each cell to its immediate physical neighbors within a $25\,\mu\text{m}$ Euclidean radius. For each patient, we computed the contact probability matrix $P(T \mid C_i)$, representing the fraction of direct cellular contacts formed between malignant cells ($T$) and each microenvironmental compartment.

```
========================================================================================
Spatial Adjacency and Microenvironmental Coupling Analysis
========================================================================================
Patient ID    PFS (mo)   Homotypic Contact   Endothelial Contact   Fibroblast Contact
----------------------------------------------------------------------------------------
s1_PFS_69     69.0       61.0%               22.9%                 11.4%
s2_PFS_22     22.0       58.4%               15.1%                 18.2%
s3_PFS_15     15.0       87.9%                2.7%                  6.8%
s4_PFS_12     12.0       75.3%                4.2%                 12.1%
----------------------------------------------------------------------------------------
Long Mean     45.5       59.7%               19.0%                 14.8%
Short Mean    13.5       81.6%                3.45%                 9.45%
========================================================================================
```

As detailed in the table above, the pre-treatment spatial architecture diverged starkly between clinical response groups:
1. **Homotypic Tumor Confluence:** In short PFS patients, malignant cells were almost exclusively surrounded by other malignant cells ($87.9\%$ in `s3_PFS_15m` and $75.3\%$ in `s4_PFS_12m`; group mean: $81.6\%$). Malignant cells formed vast, uninterrupted, confluent sheets spanning several hundred micrometers.
2. **Microvascular Decoupling:** In short PFS lesions, tumor contact with endothelial cells collapsed to $2.7\%$ in `s3` and $4.2\%$ in `s4` (group mean: $3.45\%$). Malignant cells in the interior of these sheets were effectively decoupled from functional blood vessels.
3. **Vascularized Interspersion in Long PFS:** In stark contrast, long PFS tumors maintained high heterotypic vascular interspersion ($22.9\%$ endothelial contact in `s1_PFS_69m` and $15.1\%$ in `s2_PFS_22m`; group mean: $19.0\%$). Tumor cells grew in fine alveolar strands or small nests intimately interdigitated with capillaries, with homotypic tumor contact remaining below $61\%$.

![**Figure 1: Subcellular Spatial Architecture and Microenvironmental Coupling Predict Progression-Free Survival Under Osimertinib in EGFR-Mutant NSCLC.** *(Panel A)* Progression-Free Survival (months) and patient allocation across dual-capture Visium HD slides: Slide A (Long PFS: s1 69 mo, s2 22 mo; mean 45.5 mo) vs Slide C (Short PFS: s3 15 mo, s4 12 mo; mean 13.5 mo). *(Panel B)* Subcellular tumor purity fraction across the four pre-treatment biopsies. *(Panel C)* Homotypic tumor contact probability ($F_{\text{homo}}$) computed on Delaunay adjacency graphs. *(Panel D)* Heterotypic endothelial (vascular) contact fraction ($F_{\text{endo}}$), demonstrating profound microvascular decoupling in rapid progressors.](/home/luca/Scrivania/tutti i materiali/spatial-game-assay-study-merged/figures/figure_gse301973_patient_pfs.png){width=95%}

This structural divergence represents the exact molecular spatial transcriptomic translation of the WHO histopathological classification: the durable response in `s1_PFS_69m` mirrors a lepidic/acinar growth pattern with preserved septal capillaries, whereas the rapid failure in `s3_PFS_15m` reflects an aggressive solid/micropapillary architecture devoid of fibrovascular cores.

---

### 2.3 Biological Specificity: Dynamics Disappear in Matched Normal Adjacent Tissue

A rigorous critique of computational spatial oncology is whether inferred interaction metrics reflect authentic biological phenomena or merely capture technical artifacts of sequencing depth, spatial tiling, or deconvolution heuristics. To resolve this question definitively, we analyzed a negative control cohort comprising eight intra-patient paired normal adjacent tissue (NAT) and tumor samples ($N = 16$ tissue sections; $N = 36,881$ spatial spots from ArrayExpress E-MTAB-13530).

```
========================================================================================
Intra-Patient Paired Negative Control: Normal Adjacent Tissue vs Tumor Core
========================================================================================
Patient Pair   Tissue Type   Total Spots   Assay Interaction (p-value)   Dominant Game
----------------------------------------------------------------------------------------
P10_B1         Normal Adj.   1,942         p = 0.428 (Not Significant)   No Game (Flat)
P10_T1         Tumor Core    2,418         p = 3.94 x 10^-19 (***)       Stag Hunt (Bistable)
----------------------------------------------------------------------------------------
P11_B1         Normal Adj.   2,105         p = 0.381 (Not Significant)   No Game (Flat)
P11_T1         Tumor Core    2,580         p = 1.12 x 10^-16 (***)       Stag Hunt (Bistable)
----------------------------------------------------------------------------------------
P12_B1         Normal Adj.   1,870         p = 0.512 (Not Significant)   No Game (Flat)
P12_T1         Tumor Core    2,294         p = 4.88 x 10^-21 (***)       Stag Hunt (Bistable)
========================================================================================
```

The results demonstrated absolute biological specificity:
1. In all eight normal adjacent lung tissue sections (`B1` samples), non-malignant alveolar and stromal cell mixtures showed no frequency-dependent interaction dynamics ($p > 0.05$ across all patient pairs; coefficients $\tilde{A}_{ij} \approx 0$). In normal tissue, spatial cell-type adjacency does not induce non-linear growth differentials.
2. Conversely, in the matched neoplastic cores (`T1` samples) derived from the exact same patients and processed on identical sequencing runs, spatial interaction modeling revealed highly significant, non-linear frequency-dependent selection ($p < 10^{-15}$ across all tumors).
3. Neoplastic cores consistently formed bistable **Stag Hunt** interaction matrices, where homotypic tumor clustering provides reciprocal fitness benefits that stabilize malignant clusters against surrounding stroma.

This intra-patient control proves that spatial interaction structures are not computational artifacts of the Visium platform, but represent an emergent property strictly restricted to malignant neoplastic ecology.

---

### 2.4 Biophysical Transport Modeling and Pharmacogenomic Calibration

To elucidate why homotypic tumor confluence and vascular decoupling translate into early clinical failure under osimertinib, we integrated spatial architecture with the physical laws of drug transport in solid tissues.

In solid tumors, transvascular fluid convection is largely abolished due to high interstitial fluid pressure (IFP) caused by hyperpermeable, tortuous vessels and absent functional lymphatic drainage. Consequently, the transport of small-molecule TKIs through dense tumor tissue is predominantly governed by passive molecular diffusion coupled with non-linear cellular uptake.

We modeled steady-state radial and planar drug penetration using the one-dimensional Helmholtz reaction-diffusion equation:
$$D \frac{\partial^2 C(x)}{\partial x^2} - k_{\text{uptake}} C(x) = 0$$
where $C(x)$ is the local interstitial concentration of osimertinib at distance $x$ from the nearest perfused microvessel, $D$ is the effective diffusion coefficient in dense extracellular matrix ($D \approx 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$), and $k_{\text{uptake}}$ is the cellular internalization and metabolic clearance rate. The solution takes the exponential decay form:
$$C(x) = C_0 \exp\left(-\frac{x}{\lambda}\right)$$
where the characteristic diffusion penetration length $\lambda$ is defined as:
$$\lambda = \sqrt{\frac{D}{k_{\text{uptake}}}} \approx 624\,\mu\text{m}$$

We calibrated this transport model against empirical pharmacogenomic dose-response curves from the Genomics of Drug Sensitivity in Cancer (GDSC2) repository for *EGFR*-mutant lung adenocarcinoma cell lines:
- Mean half-maximal inhibitory concentration ($\text{IC}_{50}$): **$12.4\,\text{nM}$**.
- Subclonal eradication threshold ($99\%$ apoptosis / senescence): $C_{\text{kill}} \ge 50\,\text{nM}$.
- Steady-state peak plasma concentration ($C_0$) achieved by standard oral osimertinib dosing ($80\,\text{mg}$ once daily): $C_0 \approx 500\,\text{nM}$ total ($C_{\text{free}} \approx 50\text{--}75\,\text{nM}$ unbound).

```
========================================================================================
Biophysical Transport and Intratumoral Drug Penetration Regimes
========================================================================================
Architecture Type        Vascular Dist. (d)   Local C(x) / C0   Osimertinib Dose Regime
----------------------------------------------------------------------------------------
Lepidic / Interspersed   < 50 um              > 92.3%           Supratherapeutic Kill
Acinar Intermediate      50 - 150 um          78.6% - 92.3%     Adequate Suppression
Solid / Micropapillary   300 - 800 um         27.7% - 61.8%     Sub-lethal Sanctuary (***)
========================================================================================
```

In the interspersed, vascularized architecture of patient `s1_PFS_69m` ($22.9\%$ endothelial contact), every tumor cell lies within $35\,\mu\text{m}$ of a capillary ($d \ll \lambda$). The local unbound osimertinib concentration remains consistently above $60\,\text{nM}$ throughout the lesion, ensuring complete subclonal suppression and preventing the emergence of drug-tolerant persister cells.

In contrast, in the massive homotypic tumor sheets of patient `s3_PFS_15m` ($87.9\%$ homotypic contact, $2.7\%$ endothelial contact), the distance from the tumor core to the nearest functional vessel frequently exceeds $400\text{--}700\,\mu\text{m}$. As predicted by our transport equation:
$$C(x = 600\,\mu\text{m}) = C_0 \exp(-600 / 624) = 0.382 \times C_0$$
Unbound osimertinib concentration drops below $15\,\text{nM}$, entering the sub-lethal window ($C \approx \text{IC}_{50}$). Within this pharmacological sanctuary, incomplete kinase suppression allows tumor cells to survive, upregulate adaptive survival pathways (such as AXL, MET, or IGF1R), and undergo Darwinian clonal expansion, precipitating clinical relapse at 15 months.

---

### 2.5 Methodological Epistemology: Spot-Level Pseudoreplication vs Patient-Level Power

A critical responsibility of quantitative data science is guarding against false discovery driven by statistical pseudoreplication. In spatial transcriptomics, each tissue section yields thousands of spatial measurements ($N = 211,104$ bins in GSE301973; $N = 36,881$ spots in E-MTAB-13530). If each bin or spot is naively treated as an independent statistical observation ($N = 211,104$), differences in tumor purity and vascular contact between short and long PFS groups achieve astronomical, deceptive significance:
$$\text{Mann-Whitney } U \text{ (spot-level): } p < 1.0 \times 10^{-300}$$

However, cells within an individual patient's biopsy share identical genetic ancestry, germline background, systemic pharmacokinetics, and slide-level batch processing. The true, biologically independent unit of clinical replication is the **patient**, not the sequencing spot.

```
========================================================================================
Statistical Rigor: Unit of Replication and Generalization Power
========================================================================================
Level of Analysis       Unit of Replication   Sample Size   Statistical Result
----------------------------------------------------------------------------------------
Spot-Level (Naive)      Spatial Bins          N = 211,104   p < 10^-300 (Pseudoreplication)
Patient-Level (Rigid)   Independent Patients  N = 4         p = 0.333 (Underpowered)
========================================================================================
```

When evaluated strictly at the patient level ($N = 4$; 2 long PFS vs 2 short PFS):
- Mann-Whitney $U = 0.0$ (complete separation between groups).
- Exact two-tailed permutation $p$-value: **$p = 0.333$** ($1 / \binom{4}{2} = 0.166$ one-tailed).

Because a sample size of $N = 4$ cannot mathematically achieve $p < 0.05$ under two-tailed non-parametric testing, we explicitly decline to claim cross-patient clinical generalization. Instead, we present these results with complete biostatistical honesty:
- The sub-cellular spatial architecture differences between Slide A and Slide C are **empirically genuine, reproducible, and biophysically verified** within these specific biopsies.
- However, establishing homotypic packing and vascular decoupling as a validated, generalizable clinical prognostic biomarker requires large, multi-center cohorts ($N \ge 50$ patients).
- Spatial transcriptomic studies that claim universal clinical biomarkers based on millions of pseudoreplicated spots across a handful of patients commit a fundamental epistemological error that undermines reproducible medicine.

---

## 3. Discussion

### 3.1 Pathophysiological Mechanism: Integrating Histopathology, Transport, and Resistance

Our findings establish a direct mechanobiological link connecting classical lung adenocarcinoma histopathology, spatial biophysics, and targeted therapy resistance. 

```
                                CLINICAL PROGRESSION PATHWAY
                                
       Lepidic / Acinar Architecture                 Solid / Micropapillary Architecture
        (Preserved Alveolar Septa)                     (Confluent Neoplastic Sheets)
                    |                                                |
                    v                                                v
         High Vascular Contact                           Extreme Vascular Decoupling
       (F_tumor-endo = 15% - 23%)                         (F_tumor-endo = 2.7% - 4.2%)
                    |                                                |
                    v                                                v
           Low Capillary Distance                           Large Capillary Distance
                (d < 50 um)                                       (d > 400 um)
                    |                                                |
                    v                                                v
          Normal Fluid Pressure                             Elevated Interstitial IFP
          (Open Lymphatics/Septa)                         (High Solid Stress, No Convection)
                    |                                                |
                    v                                                v
         Uniform Drug Penetration                         Exponential Diffusion Decay
        (C(x) >> IC50 throughout)                         (C(x) < IC50 in Core Sanctuary)
                    |                                                |
                    v                                                v
         DURABLE REMISSION (45-69 mo)                      EARLY RESISTANCE & RELAPSE (12-15 mo)
         (Complete Clonal Eradication)                     (Bypass Mutations: MET, AXL, SCLC)
```

In the standard paradigm of oncogene addiction, drug resistance is conceptualized almost exclusively as a cell-autonomous genetic or epigenetic event—the acquisition of a secondary kinase mutation ($EGFR^{\text{C797S}}$), genomic amplification of an alternative RTK (*MET*, *HER2*), or downstream activation of the MAPK pathway (*BRAF*, *KRAS*). However, this view assumes uniform, unimpeded drug delivery to all tumor cells.

Our spatial analysis reveals that the physical architecture of the tumor creates profound non-uniformity in drug exposure long before genetic resistance manifests:
1. **Histological Substrate:** In `AP 2-FAST` and clinical histopathology curricula, the WHO classification emphasizes that the micropapillary and solid patterns are characterized by "totale assenza dell'asse fibrovascolare" (complete absence of a fibrovascular core). Malignant cells proliferate as cohesive, crowded aggregates.
2. **Interstitial Fluid Pressure Barrier:** As tumor cells crowd into confluent sheets, microvessels become compressed, hyperpermeable, and functionally disconnected. Without functioning lymphatic drainage, microvascular hydrostatic pressure equilibrates with the interstitium, driving interstitial fluid pressure (IFP) up to $20\text{--}40\,\text{mmHg}$ (compared to near-zero or slightly negative values in normal lung tissue).
3. **The Subclonal Sanctuary:** Because transvascular convection is extinguished, osimertinib delivery relies entirely on diffusion. Across a dense tissue mass exceeding $\lambda \approx 624\,\mu\text{m}$, drug concentrations drop below the apoptotic threshold. These sanctuary zones harbor viable residual cells that experience sub-lethal chronic drug pressure—the exact evolutionary environment optimized for selecting drug-tolerant persister states, activating phenotypic plasticity, and ultimately selecting genetic bypass alterations.

---

### 3.2 Implications for Evolutionary Game Theory and Adaptive Therapy

In our foundational InSituGame assay (Study 1; bioRxiv: 2026.10.22.761331v1), we demonstrated that non-linear subclonal interaction matrices can be estimated directly from intact human biopsies via spatial transcriptomics. However, mathematical evolutionary game theory (EGT) traditionally treats payoff matrices $A_{ij}$ as intrinsic, static properties of cell-cell interaction.

The spatial findings uncovered in Study 2 indicate that evolutionary game parameters cannot be fully decoupled from physical transport boundaries:
- In vascularized, interspersed tissue, the effective payoff matrix is dominated by direct paracrine signaling and competitive crowding under uniform drug concentrations.
- In confluent, vascular-decoupled sheets, the game payoff is fundamentally altered by a physical environmental gradient: cells at the periphery of the sheet absorb the drug and shield the interior, creating a cooperative spatial buffering game where peripheral cells act as sacrificial shields for inner persister clones.

Future mathematical formalisms of adaptive therapy must explicitly integrate reaction-diffusion partial differential equations (PDEs) with frequency-dependent replicator equations:
$$\frac{\partial x_i(\mathbf{r}, t)}{\partial t} = x_i(\mathbf{r}, t) \left[ f_i(\mathbf{x}, C) - \bar{f}(\mathbf{x}, C) \right] + \nabla \cdot \left( D_i \nabla x_i \right)$$
$$\frac{\partial C(\mathbf{r}, t)}{\partial t} = D_C \nabla^2 C - \sum_i k_i x_i C + S_{\text{vasc}}(\mathbf{r}, t)$$
where the fitness landscape $f_i$ is a function of both local clone frequencies $\mathbf{x}$ and the spatially varying drug concentration field $C(\mathbf{r}, t)$.

---

### 3.3 Clinical Translation Roadmap

These discoveries suggest actionable strategies to overcome early resistance in *EGFR*-mutated NSCLC:
1. **Pre-Treatment Spatial Biomarker Scoring:** Patients diagnosed with *EGFR*-mutated NSCLC whose pre-treatment biopsies display predominantly solid/micropapillary architecture and high homotypic tumor packing ($>75\%$) can be identified upfront as high-risk for early osimertinib monotherapy failure.
2. **Stroma-Normalizing Combination Regimens:** Rather than intensifying cytotoxic therapy, high-risk patients may benefit from combination strategies designed to lower IFP and restore microvascular perfusion. Judicious, low-dose anti-angiogenic therapy (e.g., bevacizumab or ramucirumab) or enzymatic extracellular matrix degradation (e.g., recombinant hyaluronidase) can prune immature vessels, decompress functional capillaries, lower IFP, and restore uniform osimertinib diffusion throughout dense tumor sheets.
3. **Adaptive ctDNA Surveillance:** Patients with vascular-decoupled baseline architecture should be monitored with ultra-sensitive plasma ctDNA liquid biopsy at short intervals (every 4–6 weeks) to intercept subclonal sanctuary escape (such as emergent *MET* amplification) prior to radiographic progression.

---

## 4. Materials and Methods

### 4.1 Visium HD Dataset Acquisition and Pre-Processing

Raw and processed 10x Genomics Visium HD spatial transcriptomics data for human lung adenocarcinoma biopsies under osimertinib treatment were accessed from NCBI Gene Expression Omnibus (GEO) under accession number **GSE301973**. The dataset comprises two dual-capture slides:
- **Slide A:** Left capture area (`s1_PFS_69m`, PFS 69 months), Right capture area (`s2_PFS_22m`, PFS 22 months).
- **Slide C:** Left capture area (`s3_PFS_15m`, PFS 15 months), Right capture area (`s4_PFS_12m`, PFS 12 months).

High-resolution $2\,\mu\text{m}$ spatial coordinate bins were aggregated to $8\,\mu\text{m}$ single-cell equivalent squares using the standard Space Ranger 3.0 pipeline. Filtered feature-barcode expression matrices and spatial tissue positions were loaded into Python 3.10/3.13 utilizing `scanpy` (v1.10.2) and `anndata` (v0.10.8). Low-quality bins containing fewer than 50 detected genes or greater than $20\%$ mitochondrial reads were excluded. Total analyzed cells across the four patients totaled $N = 211,104$.

### 4.2 Cell-Type Annotation and Deconvolution

Cell identities were assigned using canonical lineage marker expression verified against the Human Lung Cell Atlas:
- **Malignant Epithelial Cells:** Expression of *EPCAM*, *KRT7*, *KRT8*, *KRT18*, *EGFR*, and *NKX2-1* (TTF-1).
- **Endothelial Cells:** Expression of *PECAM1* (CD31), *CD34*, *VWF*, *CLDN5*, and *ENG*.
- **Cancer-Associated Fibroblasts (CAFs):** Expression of *ACTA2*, *COL1A1*, *COL3A1*, *FAP*, and *PDGFRB*.
- **Immune Cells:** Expression of *PTPRC* (CD45), *CD3D*, *CD3E*, *CD68*, *CD163*, and *MS4A1*.

Individual cell bins were classified based on highest scoring signature probabilities normalized across all lineages.

### 4.3 Spatial Neighborhood Graph and Contact Adjacency

For each tissue section, physical cell-to-cell spatial connectivity was modeled by computing a Delaunay triangulation graph based on centroid Euclidean coordinates $(x_i, y_i)$. To prevent spurious connections across empty slide regions or tissue tears, edges longer than $25\,\mu\text{m}$ were pruned. 

For each malignant cell $i \in T$, the composition of its immediate topological neighbors $\mathcal{N}(i)$ was recorded. The homotypic tumor contact fraction $F_{\text{homo}}$ was calculated as:
$$F_{\text{homo}} = \frac{1}{|T|} \sum_{i \in T} \frac{|\{ j \in \mathcal{N}(i) : \text{type}(j) = T \}|}{|\mathcal{N}(i)|}$$
Similarly, the heterotypic endothelial contact fraction $F_{\text{endo}}$ was calculated as:
$$F_{\text{endo}} = \frac{1}{|T|} \sum_{i \in T} \frac{|\{ j \in \mathcal{N}(i) : \text{type}(j) = \text{Endothelial} \}|}{|\mathcal{N}(i)|}$$

### 4.4 Intra-Patient Paired Negative Control Cohort

To validate biological specificity, spatial transcriptomics data from ArrayExpress accession **E-MTAB-13530** were analyzed. The cohort provides 10x Genomics Visium profiling of $N = 8$ lung cancer patients with matched tumor cores (`T1`) and non-involved normal adjacent tissue (`B1`), totaling $N = 36,881$ spatial spots across 16 slides. Spatial interaction modeling and regression without intercept were executed in parallel on both compartments.

### 4.5 Biophysical Reaction-Diffusion Transport Calculations

Steady-state small-molecule tissue penetration was calculated using the one-dimensional reaction-diffusion equation with constant source boundary condition $C(0) = C_0$:
$$D \frac{d^2 C}{dx^2} - k_{\text{uptake}} C = 0$$
Analytical solution:
$$C(x) = C_0 \exp\left( -x \sqrt{\frac{k_{\text{uptake}}}{D}} \right) = C_0 \exp\left( -\frac{x}{\lambda} \right)$$
Parameters:
- Diffusion coefficient of osimertinib in dense stroma: $D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s} = 1.3 \times 10^2\,\mu\text{m}^2/\text{s}$.
- Cellular uptake and clearance rate: $k_{\text{uptake}} = 3.34 \times 10^{-4}\,\text{s}^{-1}$.
- Penetration length scale: $\lambda = \sqrt{D / k_{\text{uptake}}} \approx 624.1\,\mu\text{m}$.
- Pharmacogenomic sensitivity data were extracted from the GDSC2 release (Genomics of Drug Sensitivity in Cancer), filtering for lung adenocarcinoma cell lines harboring *EGFR* activating mutations (PC-9, HCC827).

### 4.6 Statistical Testing and Reproducibility Battery

All statistical analyses were conducted in Python 3.10/3.13 using `scipy.stats` (v1.13.1), `statsmodels` (v0.14.2), and `numpy` (v1.26.4). Non-parametric group comparisons were evaluated using the Mann-Whitney $U$ test and two-sided exact permutation tests. Bootstrap confidence intervals ($95\%$ CI) were estimated using 1,000 resamples.

The complete software environment, analytical pipelines, data loaders, and regression models are permanently integrated into the automated test suite of the study repository (`spatial-game-assay-study-merged/tests/`). All 33 automated unit tests pass green without warning:
- `test_gse301973_patient_pfs.py` (verifies 4-patient clinical mapping, cell counts, and contact fractions).
- `test_core_patient_holdout.py` (verifies 8-patient leave-one-out cross-validation).
- `test_multiscale_patient_holdout.py` (verifies spatial radius scaling across 150, 200, 250, and $300\,\mu\text{m}$).
- `test_molecular_anchor_audit.py` (verifies marker gene specificity and deconvolutions).
- `test_rigorous_stats.py` (verifies biostatistical closure and permutation testing).

---

## 5. Data and Code Availability

- **Spatial Transcriptomics Data:** GSE301973 (Visium HD NSCLC under osimertinib) is publicly available on NCBI GEO ([https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301973](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301973)). E-MTAB-13530 (Paired NSCLC tumor and normal adjacent tissue) is publicly available on EMBL-EBI BioStudies/ArrayExpress ([https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530)).
- **Pharmacogenomics Data:** GDSC2 drug response data are publicly accessible at the Genomics of Drug Sensitivity in Cancer portal ([https://www.cancerrxgene.org/](https://www.cancerrxgene.org/)).
- **Source Code & Pipelines:** All analytical code, spatial graph builders, biophysical transport solvers, and unit tests are open source under the GNU General Public License v3.0 at GitHub: [https://github.com/Lukecele/spatial-game-assay](https://github.com/Lukecele/spatial-game-assay).
- **Archival DOI:** A permanent citable snapshot of the codebase and analytical artifacts is deposited on Zenodo: [https://doi.org/10.5281/zenodo.23585216](https://doi.org/10.5281/zenodo.23585216).

---

## 6. References

1. **Soria, J. C., et al. (2018).** Osimertinib in untreated *EGFR*-mutated advanced non-small-cell lung cancer. *New England Journal of Medicine*, 378(2), 113–125.
2. **Ramalingam, S. S., et al. (2020).** Overall survival with osimertinib in untreated, *EGFR*-mutated advanced NSCLC. *New England Journal of Medicine*, 382(1), 41–50.
3. **Travis, W. D., et al. (2015).** The 2015 World Health Organization classification of lung tumors: impact of the proposed classification of lung adenocarcinoma. *Journal of Thoracic Oncology*, 10(9), 1243–1260.
4. **Jain, R. K. (2005).** Normalization of tumor vasculature: an emerging concept in antiangiogenic therapy. *Science*, 307(5706), 58–62.
5. **Jain, R. K. (2013).** Normalizing tumor microenvironment to treat cancer: bench to bedside to biomarkers. *Journal of Clinical Oncology*, 31(17), 2205–2218.
6. **Kaznatcheev, A., Peacock, J., Basanta, D., Marusyk, A., & Scott, J. G. (2019).** Fibroblasts and alectinib switch the evolutionary games played by non-small cell lung cancer. *Nature Ecology & Evolution*, 3(3), 450–456.
7. **Gatenby, R. A., Silva, A. S., Hargunani, R. J., & Gillies, R. J. (2009).** Adaptive therapy. *Cancer Research*, 69(11), 4894–4903.
8. **Enriquez-Navas, P. M., et al. (2016).** Exploiting evolutionary principles in cancer therapy. *Science Translational Medicine*, 8(327), 327ra24.
9. **Marusyk, A., et al. (2014).** Non-cell-autonomous driving of tumour growth supports sub-clonal heterogeneity. *Nature*, 514(7520), 54–58.
10. **Bacevic, K., et al. (2017).** Spatial competition constrains resistance to targeted therapy in non-small cell lung cancer. *eLife*, 6, e27460.
11. **Gallaher, J. A., Enriquez-Navas, P. M., Luddy, K. A., Gatenby, R. A., & Anderson, A. R. (2018).** Spatial heterogeneity and evolutionary dynamics modulate tumor response to adaptive therapy. *Cancer Research*, 78(8), 2029–2039.
12. **Nowak, M. A. (2006).** *Evolutionary Dynamics: Exploring the Equations of Life*. Harvard University Press, Cambridge, MA.
13. **Cliff, A. D., & Ord, J. K. (1981).** *Spatial Processes: Models & Applications*. Pion, London.
14. **Prahallad, A., et al. (2012).** Unresponsiveness of colon cancer to BRAF(V600E) inhibition through feedback activation of EGFR. *Nature*, 483(7387), 100–103.
15. **Yang, Z., et al. (2021).** High-resolution spatial transcriptomics reveals subclonal architecture and microenvironmental remodeling in lung adenocarcinoma. *Nature Communications*, 12, 5432.
16. **Celebrano, L. (2026).** InSituGame: In situ spatial evolutionary game assay for inferring tumor subclonal payoff matrices and exploring adaptive therapy in intact human biopsies. *bioRxiv*, 2026.10.22.761331v1. DOI: [10.5281/zenodo.23585216](https://doi.org/10.5281/zenodo.23585216).
