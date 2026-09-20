# T1D-sceQTL

Code repository for the study:

**"Cell type–resolved eQTL mapping reveals immune regulatory effects of type 1 diabetes risk variants in Chinese children and adolescents"**

This repository contains the analysis code used to investigate cell type-specific
genetic regulation in type 1 diabetes (T1D) by integrating single-cell RNA
sequencing, genotype, HLA typing, and T1D GWAS data.

## Overview

The study includes:

- Single-cell RNA-seq profiling of peripheral blood mononuclear cells (PBMCs)
  from individuals with T1D and healthy controls
- Cell type-specific cis-eQTL mapping
- Genotype-by-disease interaction eQTL analysis
- Integration of T1D GWAS variants with cell type-specific eQTLs
- GWAS–eQTL colocalization analysis
- Classical HLA association and HLA allele-based eQTL analysis
- Clinical phenotype analyses
- Development and validation of a genetic risk score (GRS)

## Repository structure

```text
T1D-single-cell-eQTL/
├── README.md
├── LICENSE
├── environment.yml
├── config/
│   ├── config.yaml
│   └── celltype_names.tsv
├── scripts/
│   ├── 01_scRNA/        # scRNA-seq preprocessing and cell proportion analysis
│   ├── 02_eQTL/         # cis-eQTL and genotype-by-disease interaction eQTL analysis
│   ├── 03_GWAS_coloc/   # GWAS integration and colocalization analysis
│   ├── 04_HLA/          # HLA association, HLA-eQTL and conditional analyses
│   ├── 05_clinical/     # clinical phenotype association analyses
│   └── 06_GRS/          # genetic risk score development and validation
└── demo/
    ├── B_GE.txt
    ├── Covariates.txt
    ├── SNP.txt
    ├── geneloc.txt
    └── snpsloc.txt
```

## Demo

The `demo/` directory provides a small example dataset for demonstrating the use of `scripts/02_eQTL/eQTL_linear.r`.

The example files include gene expression, genotype, covariate, gene-location, and SNP-location data.
