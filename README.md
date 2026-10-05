# hcc-singlecell

Single-cell analysis of human hepatocellular carcinoma (HCC): cellular heterogeneity across patients and tissue sites, and communication between tumour cells with progenitor/stem-like features and immunosuppressive myeloid populations.

> **Status: work in progress.** The repository is developed in the open; see [Status](#status) for what is done.

## Objectives

1. Build an integrated, annotated single-cell atlas across patients and tissue sites.
2. Define cancer stem-like and immunosuppressive myeloid states and their molecular signatures.
3. Reconstruct myeloid differentiation trajectories.
4. Infer ligand–receptor communication between these compartments and prioritise candidate mediators.
5. *(Stretch)* Map the reference onto public spatial data (Xenium/Visium) and quantify colocalisation.

## Data

| Dataset | Type | Content |
|---|---|---|
| [GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614) (Lu et al., 2022) | scRNA-seq (10x) | ~70k cells, 10 patients; primary tumour, portal vein tumour thrombus, metastatic lymph node, non-tumour liver |
| Public 10x liver cancer dataset (TBD) | Xenium / Visium | Spatial validation (stretch goal) |

The analysis starts from the author-provided **raw count matrix** and cell metadata (no FASTQ re-processing). Consequently, spliced/unspliced counts are not available and RNA velocity is out of scope; trajectories rely on `slingshot`.

Raw data are not tracked. `slurm/02_download.sh` retrieves them from GEO and records MD5 checksums.

## Analysis plan

| Step | Methods |
|---|---|
| QC | Lower-tail per-sample MAD outliers (UMIs, genes), absolute % mitochondrial cap, doublets with `scDblFinder` run per sample |
| Normalisation & integration | Log-normalisation, 3,000 HVGs (MT/ribosomal excluded), PCA, `Harmony` across **patients** |
| Integration assessment | Unintegrated vs integrated embeddings; kNN-based LISI for patient (mixing), tissue site and author cell type (should stay separated) |
| Clustering & annotation | Leiden clustering at several resolutions, UMAP, canonical markers, `presto` cluster markers, `SingleR` (HPCA) and comparison with the authors' labels |
| Malignant cells | `inferCNV` per patient; reference = the patient's immune/stromal cells plus non-tumour epithelium from the other patients; cells called on their correlation with the patient's consensus tumour CNV profile, calibrated on the epithelial reference; the patient's own non-tumour epithelium as negative check |
| Progenitor-like (stem-like candidate) cells | Malignant cells only: top quartile, within each patient, of a hepatic progenitor / biliary programme score (EPCAM, KRT19, KRT7, SOX9, CD24, PROM1); `CytoTRACE2` potency as an independent, whole-transcriptome check |
| Myeloid compartment | Re-analysis of myeloid cells alone (HVGs, PCA, Harmony, Leiden); immunosuppression programme score; tumour vs non-tumour enrichment of each state in paired patients |
| Differential expression | Pseudobulk `DESeq2` on raw counts, `~ patient + group`: progenitor-like vs other malignant cells; myeloid tumour vs non-tumour |
| Pathway enrichment | `fgsea` on Hallmark gene sets, genes ranked by the DESeq2 Wald statistic |
| Trajectories | `slingshot` on the myeloid Harmony embedding, rooted in classical monocytes; monocyte-derived states only; gene–pseudotime correlations checked per patient |
| Cell–cell communication | `LIANA` consensus (NATMI, Connectome, log2FC, SingleCellSignalR) on tumour-site cells |
| Candidate prioritisation | Progenitor-like ↔ TAM interactions ranked by the mean of two ranks, LIANA aggregate rank and specificity vs other malignant cells; DE of the tumour-side gene reported, genes of the defining programme flagged — an ordering of hypotheses for validation, not a test |
| Report | HTML report rendered as the last pipeline step (`docs/index.html`, GitHub Pages) |
| Spatial (stretch) | Reference mapping, niche detection, colocalisation statistics |

### Design decisions

- **Study design limits the comparisons.** Tumour and non-tumour tissue are paired in 8 of 10 patients, so tumour vs non-tumour tests are done within patients. Portal vein tumour thrombus (2 patients) and lymph-node metastasis (1 patient) are described, not tested.
- **QC is a re-assessment, and composition-aware.** The GEO matrix was already filtered by the authors (71,915 cells). A first pass with per-sample MAD thresholds on both tails removed 15% of cells, and not at random: upper-tail UMI outliers were 8.8% of non-tumour cells vs 0.8% of tumour cells. Liver samples mix small immune cells with large, mitochondria-rich hepatocytes, so the upper tail is a cell type rather than an artefact. Final rules: lower-tail MAD outliers (UMIs, genes) per sample, a single absolute mitochondrial cap (30%), and doublets removed with scDblFinder instead of an upper UMI cut-off. The rejected rules are still computed and tabulated by cell type (`results/tables/02_qc_flags_by_celltype.tsv`): a per-sample upper MAD on % mitochondrial would have removed 13% of hepatocytes (and 11% of endothelial cells) but 3% of T/NK cells, and the upper UMI rule 8–9% of myeloid, B and endothelial cells. With the final rules, 67,908 of 71,915 cells are kept; almost all removals are scDblFinder doublets (5.6%).
- **Log-normalisation rather than SCTransform.** Differential expression is done on raw counts at the pseudobulk level, so SCTransform would only affect the embedding, at a memory cost that matters on a cluster limited to 2 GB per CPU.
- **Integrate across patients, not tissue sites.** Tissue site is the biological signal of interest; correcting it would remove it.
- **Differential expression at the pseudobulk level** with the patient as blocking factor, so that tumour vs non-tumour comparisons are made within patients and cells are not treated as independent replicates.
- **Stem-like states are treated cautiously.** Cancer stem cells in HCC are a debated concept and marker scores are weak evidence on their own. A first definition (top quartile of a CSC marker score AND of a within-sample gene-count rank) was dropped after the second run: the two readouts were uncorrelated in every patient (Spearman −0.17 to 0.28), so their intersection (~6% of cells) was what chance alone gives, and differential expression only returned the defining markers. The state is now defined by one programme (hepatic progenitor / biliary genes) and named for it, "progenitor-like"; CytoTRACE2, which does not use the marker list, tests whether it also has higher developmental potential. The rule was fixed before the result: the state is called stem-like only where CytoTRACE2 agrees. It did not (Results §7), so the state is reported as progenitor/biliary-like.
- **No circular evidence in the candidate ranking.** Differential expression of a gene that defines the group cannot support that gene as a mediator, so DE is reported but not scored, and defining genes are flagged.
- **Myeloid-derived suppressor cells** are hard to separate from monocytes/neutrophils in scRNA-seq alone; states are named by function-associated programmes rather than by assumed identity.

## Results

### 1. Data and study design

71,915 cells from 21 samples (10 patients) were read from the author-provided count matrix; all cells matched the metadata. Tumour and non-tumour tissue are paired in 8 patients; PVTT is available for 2 patients and lymph-node metastasis for 1 ([`01_cells_per_patient_site.tsv`](results/tables/01_cells_per_patient_site.tsv)).

### 2. Quality control

The matrix had already been filtered by the authors: mitochondrial content stops sharply at 20% and UMIs have a hard lower bound in every sample. QC therefore re-assesses that filtering rather than repeating it.

![QC per sample](results/figures/02_qc_per_sample.png)

Several tumour samples show bimodal UMI distributions (small immune cells vs large tumour/hepatocyte-like cells). Per-sample MAD thresholds on the upper tail would have removed cells by type rather than by quality:

| Author cell type | Cells | Upper-MAD % mito rule | Upper-MAD UMI rule |
|---|---:|---:|---:|
| Hepatocyte (incl. malignant) | 20,782 | 13.0% | 1.8% |
| Endothelial | 3,644 | 11.2% | 8.4% |
| Fibroblast | 2,266 | 9.1% | 4.2% |
| B | 3,685 | 5.3% | 9.3% |
| Myeloid | 15,947 | 4.2% | 8.8% |
| T/NK | 25,591 | 3.0% | 1.1% |

With the final rules, **67,908 cells** are kept; removals are almost entirely scDblFinder doublets (3–7% per sample) ([`02_qc_summary.tsv`](results/tables/02_qc_summary.tsv), [`02_qc_flags_by_celltype.tsv`](results/tables/02_qc_flags_by_celltype.tsv)).

### 3. Integration across patients

Harmony (batch = patient) raised the median patient LISI from 1.1 to 3.0 (10 patients), while cell-type LISI stayed at 1.0, i.e. patients mix within cell types and cell types are not merged. Tissue-site LISI barely moved (1.0 → 1.07), as intended: site is biology, not batch ([`03_integration_lisi.tsv`](results/tables/03_integration_lisi.tsv)).

![Integration](results/figures/03_integration_umap.png)

Immune and stromal cells mix across patients after integration; hepatocyte-lineage cells remain partly patient-specific (one cluster is 67% from a single patient). That is expected for malignant cells, whose copy-number profiles differ between tumours, and is why malignant cells will be identified with inferCNV rather than by forcing further integration.

### 4. Clustering and first-pass annotation

Leiden clustering (resolution 0.5) gives 23 clusters; majority SingleR (HPCA) labels and the authors' labels agree with the markers for all major compartments ([`04_cluster_annotation.tsv`](results/tables/04_cluster_annotation.tsv), [`04_cluster_markers.tsv`](results/tables/04_cluster_markers.tsv)).

![Annotation](results/figures/04_annotation_umap.png)

![Markers](results/figures/04_marker_dotplot.png)

| Compartment | Clusters | Key markers |
|---|---|---|
| T / NK | 0, 2, 9, 12 (cycling) | CD3D, NKG7, GZMA, IL32; MKI67/TOP2A in 12 |
| Myeloid | 1, 5, 7, 18 (cycling), 22 | SPP1, CD68 (1); LYZ, HLA-DR (5); CD5L, C1QA, CD163 (7, Kupffer-like); LILRA4, IRF4 (22, pDC) |
| B / plasma | 14, 15 | CD79A; MZB1, JCHAIN, IGHG1 |
| Endothelial | 8, 16 | PECAM1, VWF; CLEC4G, FCN2, CRHBP (16, sinusoidal) |
| Fibroblast / stellate | 10 | ACTA2, TAGLN, COL1A2 |
| Hepatocyte-like | 3, 11, 13 | ALB, APOA1, TTR, ORM1, GSTA1 |
| Tumour-like | 4, 6 (cycling), 17 | GPC3, MDK, CD24, KRT8/18; SPINK1, AGR2 (17) |
| Cholangiocyte | 21 | EPCAM, KRT19, KRT7 |

Two observations guide the next steps: an **SPP1⁺ macrophage** cluster (1), a candidate immunosuppressive TAM population, and a **CD24⁺ MDK⁺ epithelial** cluster (4) with a proliferating counterpart (6), candidate stem-like tumour states. Both are hypotheses until malignant cells are confirmed by CNV and the myeloid compartment is sub-clustered. Clusters 19 (plasma-like, 76% one patient) and 20 (T cells with hepatocyte transcripts) are flagged as possible residual doublets or ambient RNA.

### 5. Malignant cells (inferCNV) — first pass and revision

inferCNV was run per patient (9 patients with ≥ 50 epithelial cells). In five patients malignant and reference cells separate cleanly (e.g. HCC02, HCC04, HCC05, HCC09: 93–100% of tumour epithelium called malignant, non-tumour epithelium of HCC04 and HCC09 at 0–2%).

![CNV first pass](results/figures/06_cnv_scores.png)

The first pass also exposed a design flaw: with immune/stromal cells as the only reference, non-tumour epithelium was called malignant in 31–58% of cells in four patients (HCC03, HCC05, HCC06, HCC10). Hepatocyte-specific genes cluster along some chromosomes, so lineage differences read as copy-number change. The reference now also includes non-tumour epithelium pooled from the other patients (leave-one-patient-out), and thresholds are calibrated on that lineage-matched reference; each patient's own non-tumour epithelium stays out and remains the negative check.

The second run fixed the false positives (non-tumour epithelium 0–0.5% malignant in HCC03, 04, 05, 06 and 09) but exposed a second problem: requiring both the CNV score and the correlation to exceed the 99th percentile of the epithelial reference left most tumour epithelium "unresolved" in some patients (HCC03 1.4% malignant, HCC05 2.6%). The epithelial reference comes from nine other livers, so its CNV score is widely spread and overlaps the tumour; the correlation with the patient's own tumour profile separates the two in every patient. Cells are now called on that correlation alone (malignant above the 99th percentile of the epithelial reference, non-malignant below the 95th).

Final calls ([`06_cnv_calls_by_sample.tsv`](results/tables/06_cnv_calls_by_sample.tsv)): **98.4–100% of tumour, PVTT and lymph-node epithelium is malignant in every patient**, and non-tumour epithelium is called malignant in 0–5% of cells where it is well represented (HCC03 0.1%, HCC04 0%, HCC05 0%, HCC09 2.9%, HCC06 4.9%). Three non-tumour samples have few epithelial cells and higher rates: HCC07 (4 of 34), HCC08 (3 of 19) and **HCC10 (35 of 70, 50%)**. HCC10 is the only clear exception; tumour cells in the adjacent sample (the patient also has a lymph-node metastasis) are a plausible explanation that the data cannot confirm. 18,475 malignant cells go forward.

### 6. Myeloid states

Re-analysis of 14,959 myeloid cells gives 12 states. Tumour enrichment is the median, over the 8 patients with paired tumour and non-tumour tissue, of log2(fraction of the patient's tumour myeloid cells in the state / fraction in non-tumour) ([`07_myeloid_states.tsv`](results/tables/07_myeloid_states.tsv)).

![Myeloid states](results/figures/07_myeloid_states.png)

| State | Identity (markers) | Cells | Tumour-site cells | log2 tumour / non-tumour | Immunosuppression score |
|---|---|---:|---:|---:|---:|
| 3 | SPP1⁺ TAM (SPP1, LGALS1, S100A10) | 2,002 | 98% | +1.6 | 0.20 |
| 2 | C1Q⁺ APOE⁺ TAM (C1QC, APOE, APOC1, TREM2) | 2,109 | 96% | +4.6 | 0.19 |
| 0 | APOE⁺ LGMN⁺ TAM (RNASE1, APOE, LGMN, A2M) | 2,658 | 97% | +3.8 | 0.19 |
| 9 | Cycling macrophages (MKI67, STMN1) | 456 | 93% | +1.9 | 0.08 |
| 5 | Stress-high monocyte-like (HSP genes, FCN1) | 1,134 | 87% | +2.2 | −0.05 |
| 1 | Kupffer cells (CD5L, MARCO, VCAM1, CD163) | 2,542 | 3% | −5.8 | 0.13 |
| 4 | Classical monocytes (FCN1, S100A4) | 1,240 | 23% | −3.3 | −0.22 |
| 6 | cDC2 / mature DC (CD1C, LAMP3, CCR7) | 1,068 | 46% | −2.2 | −0.17 |
| 8 | cDC1 (CLEC9A, IDO1) | 479 | 33% | −2.9 | −0.18 |
| 7 | MHC-II-high macrophages / DC | 926 | 52% | −1.9 | −0.03 |

The three tumour-associated macrophage states (0, 2, 3; ~6,800 cells) are almost confined to tumour sites (96–98% of their cells), enriched over paired non-tumour liver in the median patient, and carry the highest immunosuppression score; non-tumour liver is dominated by Kupffer cells and monocytes. SPP1⁺ TAMs (state 3) are the prime candidate partner for progenitor-like tumour cells in the communication analysis. Caveats: state 5 carries a dissociation-stress signature, and states 10 (Kupffer cells with hepatocyte transcripts, 264 cells) and 11 (myeloid–T doublets, 81 cells) are treated as artefacts.

Pseudobulk DE of all myeloid cells, tumour vs non-tumour within patients, agrees: tumour myeloid cells up-regulate a lipid-handling programme (APOE, LPL, FABP5; Hallmark cholesterol homeostasis, adipogenesis, fatty-acid metabolism), proliferation (E2F, G2M) and SPP1/MMP14 (EMT, angiogenesis), and down-regulate interferon-γ/α responses. NK/T cytotoxic genes (GZMB, PRF1, IL2RB) also come out lower in tumour myeloid cells, which points to some NK/T contamination of non-tumour myeloid profiles rather than myeloid biology ([`09_de_myeloid_tumour_vs_normal.tsv`](results/tables/09_de_myeloid_tumour_vs_normal.tsv), [`09_hallmark_gsea.tsv`](results/tables/09_hallmark_gsea.tsv)).

### 7. Progenitor-like tumour cells: no stem-like signal

Within the 18,475 malignant cells, the top quartile per patient of the progenitor/biliary programme (EPCAM, KRT19, KRT7, SOX9, CD24, PROM1) defines the progenitor-like group. CytoTRACE2, run independently on the whole transcriptome, does **not** give these cells higher developmental potential ([`08_progenitor_by_patient.tsv`](results/tables/08_progenitor_by_patient.tsv)):

| Patient | Malignant cells | Spearman programme vs CytoTRACE2 | Median CytoTRACE2, progenitor-like vs other |
|---|---:|---:|---|
| HCC02 | 2,073 | −0.02 | 0.077 vs 0.076 |
| HCC03 | 3,394 | −0.02 | 0.022 vs 0.020 |
| HCC04 | 2,442 | −0.08 | 0.067 vs 0.068 |
| HCC05 | 1,186 | 0.05 | 0.045 vs 0.044 |
| HCC06 | 237 | 0.33 | 0.080 vs 0.069 |
| HCC07 | 193 | 0.12 | 0.117 vs 0.113 |
| HCC08 | 4,928 | 0.10 | 0.257 vs 0.243 |
| HCC09 | 1,459 | 0.21 | 0.246 vs 0.233 |
| HCC10 | 2,563 | **−0.24** | **0.132 vs 0.249** |

Correlations are near zero in most patients, the differences in median potency are ≤ 0.014 where positive, and in HCC10 the relation is reversed. By the rule set in advance, the state is therefore not described as stem-like. What it is, is shown by the DE below: a progenitor/biliary-like malignant state, the transcriptional phenotype associated with CK19⁺ HCC.

![Progenitor-like](results/figures/08_progenitor_like.png)

### 8. Differential expression and pathways

Pseudobulk DESeq2, progenitor-like vs other malignant cells within each of 9 patients ([`09_de_progenitor_vs_other_malignant.tsv`](results/tables/09_de_progenitor_vs_other_malignant.tsv)): 19 genes at padj < 0.05. Up: SOX9, KRT7, CD24, KRT19 (defining genes, so circular) and, outside the definition, CXCL6, CCL28, SPINT1, SOX4, GP2. Down: proliferation genes (CDK1, RRM2, RAD51AP1, SPC25).

Hallmark GSEA gives a consistent picture: down E2F targets, G2M checkpoint, MYC targets and oxidative phosphorylation, and the hepatocyte metabolic programmes (bile acid, fatty acid, xenobiotic metabolism); up TNFα/NF-κB, epithelial–mesenchymal transition, hypoxia, inflammatory response, TGF-β and IL-6/JAK/STAT3 signalling. The progenitor-like cells are **less proliferative** than other malignant cells, have lost hepatocyte identity and carry an inflammatory, EMT-like programme.

![GSEA](results/figures/09_hallmark_gsea.png)

### 9. Communication and candidate mediators

LIANA consensus on tumour-site cells, then progenitor-like ↔ TAM interactions ranked by LIANA rank and specificity vs other malignant cells ([`10_prioritised_mediators.tsv`](results/tables/10_prioritised_mediators.tsv)).

![Mediators](results/figures/10_prioritised_mediators.png)

- **CD24 → SIGLEC10** (progenitor-like → C1Q⁺APOE⁺ TAM) ranks first, with by far the highest specificity (LIANA rank 0.0025 with progenitor-like cells vs 0.14 with other malignant cells). CD24–Siglec-10 is a described "don't eat me" checkpoint that inhibits macrophage phagocytosis (Barkal et al., *Nature* 2019), and SIGLEC10 is part of the immunosuppression programme used here. **Caveat:** CD24 is one of the six genes that define the group, so its specificity is inflated by construction; the interaction is a hypothesis to test with a definition that excludes CD24.
- The remaining candidates (SPON2–ITGB2, SERPINE1–PLAUR, CD59–STAB1, SPP1–CD44, TNFSF12–TNFRSF12A, CXCL12–ITGB1, APOE–SCARB1) have low specificity, and the tumour-side gene is not differentially expressed. They are ordered hypotheses, not findings.

### 10. Myeloid trajectories

*(Running.)* Slingshot from classical monocytes to the TAM states; Kupffer cells (embryonic resident macrophages) and dendritic cells are excluded.

## Computing environment

The pipeline runs on a shared SLURM cluster (CentOS 7, glibc 2.17, no root access, per-user storage quota).

- **Environment:** conda-forge/bioconda via a user-level `micromamba` binary, defined in [`environment.yml`](environment.yml) (R 4.5, Seurat 5). Packages only available on GitHub (LIANA, CellChat, presto, CytoTRACE2) are installed by `slurm/03_install_env.sh` and `slurm/06_install_cytotrace2.sh` (CytoTRACE2 depends on HiClimR, which imports ncdf4; ncdf4 cannot be built against the system netCDF of CentOS 7, so HiClimR is installed without it — it only uses ncdf4 to read climate files), which also exports the exact versions (`environment.lock.yml`) and commit SHAs (`github_packages.tsv`).
- **Pipeline:** [`targets`](https://docs.ropensci.org/targets/) inside a single SLURM job launched by `slurm/run_all.sh`, which checks and installs missing packages, validates the pipeline and runs it. Only outdated targets are re-run, and a failing target does not stop independent ones (`error = "continue"`).
- **Storage budget:** ≤ 60 GB for the whole project (data, environment and target store); intermediate objects stored with `qs2`.

### Setup

```bash
cd /work/<user>/hcc-singlecell
sbatch slurm/00_check_env.sh      # cluster diagnostics (modules, network, quota, partitions)
sbatch slurm/01_test_conda.sh     # dry-run: can the environment be solved on this system?
sbatch slurm/02_download.sh       # GEO download + checksums
sbatch slurm/03_install_env.sh    # create env/, install GitHub packages, export lock files
sbatch slurm/05_download_refs.sh  # hg38 gene positions for inferCNV
sbatch slurm/06_install_cytotrace2.sh  # CytoTRACE2 (potency check for progenitor-like cells)
sbatch slurm/run_all.sh           # check packages, run or resume the whole pipeline, render the report
```

Each script writes a plain-text report (`*_report.txt`) so that runs can be inspected without an interactive session.

## Repository structure

```
hcc-singlecell/
├── README.md
├── environment.yml       # conda specification
├── _targets.R            # pipeline definition
├── R/                    # functions used by the pipeline (ingest, QC, ...)
├── slurm/                # setup and launch scripts
├── results/
│   ├── figures/
│   └── tables/
├── report/
│   └── index.Rmd         # HTML report, rendered by the pipeline
└── docs/
    └── index.html        # rendered report (GitHub Pages)
```

## Reproducibility

- Conda environment specification plus exported lock file; GitHub packages pinned by commit SHA
- Pipeline orchestrated with `targets` (cached, re-runs only what changed)
- Fixed random seeds for clustering and UMAP
- Session info included in the report

## Status

- [x] Cluster diagnostics and environment feasibility (conda solve on glibc 2.17)
- [x] Data download with checksums
- [x] Environment installation and lock files
- [x] QC (67,908 cells kept)
- [x] Integration and first-pass annotation
- [x] Malignant-cell identification (inferCNV, lineage-matched reference)
- [x] Myeloid states
- [x] Progenitor-like tumour states (CytoTRACE2 does not support a stem-like reading)
- [x] Differential expression and pathways
- [x] Cell–cell communication and candidate ranking
- [ ] Trajectories *(running)*
- [ ] Report published on GitHub Pages *(rendered by the current run)*
- [ ] Sensitivity check: progenitor-like definition without CD24
- [ ] Spatial mapping (stretch)

## References

- Lu Y, Yang A, Quan C, et al. A single-cell atlas of the multicellular ecosystem of primary and metastatic hepatocellular carcinoma. *Nature Communications* 13, 4594 (2022).
- Barkal AA, Brewer RE, Markovic M, et al. CD24 signalling through macrophage Siglec-10 is a target for cancer immunotherapy. *Nature* 572, 392–396 (2019).

## Author

Antonio Puriel Hernández — PhD in Environmental Microbiology; bioinformatics (metagenomics, reproducible pipelines).
