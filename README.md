# T1D single-cell eQTL analysis

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
└── scripts/
    ├── 01_scRNA/
    ├── 02_eQTL/
    ├── 03_GWAS_coloc/
    ├── 04_HLA/
    ├── 05_clinical/
    └── 06_GRS/
