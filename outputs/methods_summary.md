# Methods summary

This repository contains a demo run of a PharmCAT-like pharmacogenetic VCF interpreter.

## Inputs

- `input/Pharma_subset.vcf`: small VCF example for repository demonstration.
- `input/dataset_full.csv`: manually curated pharmacogenetic rules.
- `input/PharmaVEP_final.csv`: VEP-derived annotation table used to map rsID to genomic coordinates.

## Method

1. Normalize pharmacogenetic rules.
2. Map rule rsIDs to `CHROM:POS:REF:ALT` coordinates.
3. Parse target variants from VCF.
4. Convert patient genotypes into marker calls.
5. Infer marker-based star alleles and simplified diplotypes.
6. Translate diplotypes/statuses into phenotypes.
7. Attach drug recommendations and actionability labels.
8. Save reports, coverage tables, run metadata and figures.

## Limitations

This is a research prototype. It does not perform full haplotype phasing, full CNV calling, structural variant calling or direct HLA typing. HLA-B*58:01 is handled by low-confidence proxy SNP logic.
