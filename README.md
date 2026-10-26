# Family-wide screening of RhoGAP genes in prostate adenocarcinoma

Analysis code for the manuscript:

> **Family-Wide Transcriptomic Screening of Rho GTPase-Activating Proteins Prioritizes
> *ARHGAP11A* as a Candidate Grade-Associated Member and Associates It with a
> Th2-Skewed Immune Microenvironment in Prostate Adenocarcinoma**

All results in the manuscript can be reproduced from public TCGA-PRAD data with the
scripts in this repository. No patient-identifiable or restricted data are included.

---

## 1. Data sources

All raw data are openly available; nothing in this repository needs to be requested.

| Layer | Source | Access |
|---|---|---|
| RNA-seq (STAR-FPKM, 60,660 genes x 554 samples) | UCSC Xena, GDC TCGA-PRAD hub | <https://xenabrowser.net/datapages/> |
| DNA methylation (Illumina HM450, 486,427 probes x 553 samples) | UCSC Xena, GDC TCGA-PRAD hub | same |
| Clinical / Gleason / pathologic T / survival | UCSC Xena, GDC TCGA-PRAD hub | same |
| REACTOME gene sets | MSigDB `c2.cp.reactome.v2024.1.Hs.symbols.gmt` | <https://www.gsea-msigdb.org/> |
| ESTIMATE signatures (141 immune + 141 stromal genes) | Yoshihara et al., *Nat Commun* 2013 | included as `data/estimate_signatures.csv` |
| KEGG pathway mapping | rest.kegg.jp (snapshot) | included as `data/kegg_link_pathway.txt` |

Scope used throughout: **483 primary tumours + 51 adjacent normals** (the 6 metastatic
samples present in the Xena release were excluded).

Two large files are **not** shipped and must be downloaded into `data/` before running
the R pipeline from raw counts:

```
data/TCGA-PRAD.star_fpkm.tsv.gz      # Xena: "TCGA-PRAD.star_fpkm"
data/c2.cp.reactome.gmt              # MSigDB C2 CP REACTOME, symbols
```

Everything downstream of expression quantification (the derived per-sample matrices in
`data/` and `data/r_analysis/data/`) is included, so the whole statistical pipeline runs
out of the box.

---

## 2. Repository layout

```
data/
  top_tier_permutation.py      # family screen, max-statistic permutation, Williams' t, bootstrap
  estimate_sensitivity.py      # ESTIMATE scores + partial-correlation sensitivity analysis
  export_for_R.py              # exports the derived matrices consumed by the R suite
  redraw_fig1.py ... fig5.py   # figure generation
  redraw_methylation.py        # Figure 3 (CpG-level expression-methylation correlation)
  redraw_volcano.py            # Supplementary Figure S1
  rhogap_expr_all.npz          # 45 RhoGAP x 554 samples, log2(FPKM+1)
  rhogap_family_*.csv          # per-gene family statistics (DE, Cox, clinic)
  rhogap_meth_analysis.csv     # CpG-level methylation vs expression
  estimate_scores.csv          # per-sample ESTIMATE immune/stromal scores
  estimate_signatures.csv      # ESTIMATE gene lists
  kegg_link_pathway.txt        # KEGG snapshot (offline enrichment)
  kegg_list_pathway.txt
  r_analysis/
    01_DE.R                    # limma differential expression
    02_family_screen.R         # independent re-implementation of the family screen
    03_clinical_OR.R           # grade-stratified expression, logistic ORs
    04_survival_ROC.R          # Kaplan-Meier, Cox, ROC/AUC
    05_immune.R                # ssGSEA immune correlations, partial Spearman
    06_methylation.R           # methylation vs expression
    07_GSEA.R                  # fgsea on REACTOME
    08_GO_KEGG.R               # clusterProfiler GO / KEGG
    08b_KEGG_local.R           # KEGG via offline TERM2GENE (no network)
    data/                      # derived matrices consumed by the scripts above
results/
  top_tier_stats.json          # family screen + permutation + bootstrap output
  estimate_sensitivity.json    # partial-correlation output
```

---

## 3. How to run

### Python (3.10+)

```bash
pip install -r requirements.txt

python data/top_tier_permutation.py     # -> results/top_tier_stats.json
python data/estimate_sensitivity.py     # -> results/estimate_sensitivity.json
python data/redraw_fig2.py              # figures 1-6 and S1
```

### R (4.4+)

Run the scripts **from the repository root**, in order:

```r
source("data/r_analysis/01_DE.R")
source("data/r_analysis/02_family_screen.R")
source("data/r_analysis/03_clinical_OR.R")
source("data/r_analysis/04_survival_ROC.R")
source("data/r_analysis/05_immune.R")
source("data/r_analysis/06_methylation.R")
source("data/r_analysis/07_GSEA.R")
source("data/r_analysis/08_GO_KEGG.R")      # needs internet
source("data/r_analysis/08b_KEGG_local.R")  # offline fallback
```

Required R packages: `limma`, `edgeR`, `fgsea`, `clusterProfiler`, `org.Hs.eg.db`,
`GO.db`, `DOSE`, `survival`, `pROC`, `GSVA`, `pheatmap`.
`08_GO_KEGG.R` queries `rest.kegg.jp` and may time out behind a firewall;
`08b_KEGG_local.R` reproduces the same enrichment from the bundled KEGG snapshot and is
the version used for the reported results.

---

## 4. Reproducibility notes

- **Fixed seed** `20261026` for every permutation and bootstrap procedure.
- **Family screen**: 45 RhoGAP genes were identified from the HGNC-approved list, 44 were
  evaluable in TCGA-PRAD (ARHGAP11B is not expressed/annotated in this release).
  Spearman rho was computed between each gene and ISUP grade.
- **Max-statistic permutation**: grade labels were shuffled 10,000 times and the top-*k*
  gene was re-selected in every iteration (*k* = 1, 2, 3); the empirical P value uses the
  add-one correction, `P = (b + 1) / (N + 1)`, so the minimum attainable value is 0.0001.
- **Modified Z score** follows Iglewicz & Hoaglin,
  `0.6745 x (rho - median) / MAD`, with family median rho = 0.046 and MAD = 0.0987.
- **Dependent correlations** (top-1 vs top-2) were compared with Williams' *t* test and a
  1,000-replicate bootstrap of the rho difference.
- **Multiple testing**: Benjamini-Hochberg where applicable. The exploratory immune-cell
  correlation analysis (24 cell types) is reported **unadjusted**, as stated in the paper.
- **Survival** is restricted to univariable Kaplan-Meier and Cox estimation because only
  11 death events are recorded in TCGA-PRAD.

---

## 5. Key results

| Analysis | Result |
|---|---|
| Family screen | *ARHGAP11A* rho = 0.387 (FDR = 4.3e-17); *SRGAP2* rho = 0.364; *ARHGAP33* rho = 0.315 |
| Max-statistic permutation | P < 0.0001 for *k* = 1, 2, 3 (0/10,000 permutations exceeded the observed statistic) |
| Modified Z | 2.33 / 2.17 / 1.84 (none reached the conventional 3.5 outlier threshold) |
| Top-1 vs top-2 | delta rho = 0.023, bootstrap 95% CI -0.074 to 0.111, P = 0.64; Williams' t = 0.53, P = 0.60 |
| Diagnostic AUC | 0.736 (95% CI 0.663-0.804) |
| *ARHGAP11A* - Th2 | rho = 0.702; partial rho 0.46 / 0.68 / 0.71 after adjusting for MKI67 / ISUP grade / ESTIMATE scores |
| *SRGAP2* parallel analysis | DE FDR = 0.191; Th2 rho = 0.228; OS HR = 1.28 (P = 0.68) |

---

## 6. Licence and citation

Code is released under the MIT Licence (see `LICENSE`). Derived tables in `data/` and
`results/` are released under CC BY 4.0. TCGA data remain subject to the NIH Genomic Data
Sharing Policy and the originating GDC terms.

If you use this code, please cite the manuscript above. A archived snapshot with a DOI
will be deposited in Zenodo upon acceptance.

## 7. Contact

Corresponding author: Yankun Chong — `cyk1525@yzu.edu.cn`
