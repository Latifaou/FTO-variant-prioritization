# FTO Variant Prioritization

## Overview

This project aims to develop a reproducible bioinformatics workflow for
the functional annotation and prioritization of non-coding variants at
the FTO obesity-associated locus.

The project is based on my Master's research in Medical Biotechnology
and Bioinformatics.

## Biological Background

Genome-wide association studies (GWAS) have identified multiple variants
within the FTO locus associated with obesity-related traits. Because many
of these variants are located in non-coding regions, functional annotation
is necessary to identify variants with potential regulatory effects.

This project reconstructs and extends my Master's thesis workflow using
reproducible Python-based analyses.

## Current Workflow

1. Load the initial SNP dataset.
2. Inspect the number of variants, columns, missing values, and data types.
3. Correct the decimal-comma formatting of P-values during data import.
4. Export a standardized dataset for downstream analysis.

## Project Structure

    FTO-variant-prioritization/
    ├── data/
    │   ├── raw/
    │   │   └── fuma_input_snps.txt
    │   └── processed/
    │       └── fuma_input_snps_clean.csv
    ├── scripts/
    │   └── 01_inspect_variants.py
    ├── results/
    ├── docs/
    ├── .gitignore
    └── README.md

## Current Status

The initial data inspection and preprocessing stage is complete.

Further stages will integrate functional genomic annotations and
regulatory evidence for variant prioritization.
