# FTO Variant Prioritization

## Overview

This project implements a reproducible bioinformatics workflow for the
functional annotation and prioritization of non-coding variants at the
FTO obesity-associated locus.

The workflow is based on my Master's research in Medical Biotechnology
and Bioinformatics and reconstructs the analysis using Python to integrate
GWAS association signals with functional and regulatory genomic evidence.

The objective is to identify candidate variants supported by multiple
complementary annotation layers rather than relying on a single score
or annotation.

## Biological Background

Genome-wide association studies (GWAS) have identified multiple variants
within the FTO locus associated with body mass index (BMI) and obesity.

Many obesity-associated variants at this locus are located in non-coding
regions, making their functional interpretation challenging. Regulatory
annotations can therefore help identify variants that may influence
gene regulation rather than protein sequence.

The FTO region is particularly interesting because regulatory elements
within the locus have been linked to genes involved in adipocyte biology,
including IRX3 and IRX5.

This project integrates several types of evidence to prioritize variants
for further functional investigation.

## Workflow

The analysis consists of the following stages:

1. **Input variant preprocessing**
   - Load the initial GWAS-associated SNP dataset.
   - Standardize P-value formatting.
   - Inspect missing values and data types.

2. **FUMA SNP annotation processing**
   - Process FUMA-expanded SNP annotations.
   - Examine CADD scores, RegulomeDB categories, chromatin states,
     linkage disequilibrium information, and mapping evidence.

3. **eQTL analysis**
   - Process FUMA eQTL results.
   - Identify variants associated with gene-expression differences
     across available tissues.

4. **Chromatin-interaction analysis**
   - Process FUMA chromatin-interaction results.
   - Map Ensembl gene identifiers to gene symbols.
   - Examine chromatin contacts involving IRX3 and IRX5.
   - Expand regional SNP lists for variant-level integration.

5. **Integrated variant prioritization**
   - Combine GWAS, LD, CADD, RegulomeDB, eQTL, chromatin-interaction,
     and ChromHMM evidence.
   - Identify variants supported by convergent regulatory annotations.
   - Examine adipose-relevant chromatin states.

## Variant Prioritization Strategy

Rather than constructing an arbitrary weighted score, the workflow uses
a transparent evidence-integration approach.

An exploratory regulatory shortlist was defined using:

- CADD score >= 15
- RegulomeDB category 2a, 2b, or 3a

Variants were required to satisfy **both criteria** to enter the final
regulatory candidate set.

These thresholds are exploratory prioritization criteria and should not
be interpreted as validated causal thresholds.

Additional evidence, including GWAS association, linkage disequilibrium,
eQTL results, chromatin interactions, and adipose-relevant ChromHMM states,
is retained to support biological interpretation.

## Results

FUMA annotation produced a set of **223 expanded variants**, including
**14 independent significant SNPs**.

Of the 223 variants:

- 30 had direct GWAS P-values available in the processed FUMA table.
- 7 had eQTL evidence in the available FUMA results.
- 5 satisfied the convergent CADD and RegulomeDB prioritization criteria.

The five prioritized regulatory candidates were:

- rs77276051
- rs16952522
- rs73612011
- rs3751812
- rs7202116

Among these candidates, **rs16952522** was itself an independent
significant SNP, had a direct GWAS association (P = 8 × 10^-22),
a CADD score of 17.38, and a RegulomeDB category of 2b.

Both rs16952522 and rs77276051 were annotated as enhancer-state variants
across the three adipose-relevant epigenomic datasets examined in this
workflow.

## Chromatin Interaction Context

Chromatin-interaction mapping identified regions contacting IRX3 and IRX5.

All 223 FUMA-expanded variants occurred within regions carrying
IRX3/IRX5 chromatin-interaction evidence in the analyzed dataset.
Consequently, this evidence provides biological context but does not
discriminate between individual variants in the current variant set.

Chromatin-interaction evidence should therefore not be interpreted as
proof that an individual SNP directly regulates IRX3 or IRX5.

## eQTL Context

Seven variants had eQTL evidence in the available FUMA results, involving
RBL2 and AKTIP.

None of the five final regulatory candidates had eQTL evidence in this
dataset.

The absence of an eQTL association should not be interpreted as evidence
that a variant lacks regulatory activity, because eQTL detection depends
on tissue, cell state, sample size, and the datasets available for analysis.

## Project Structure

    FTO-variant-prioritization/
    ├── data/
    │   ├── raw/
    │   │   └── fuma_input_snps.txt
    │   ├── processed/
    │   │   ├── fuma_input_snps_clean.csv
    │   │   ├── fuma_snps_clean.csv
    │   │   ├── fuma_eqtl_clean.csv
    │   │   ├── integrated_eqtl_annotations.csv
    │   │   ├── fuma_chromatin_interactions.csv
    │   │   └── prioritized_candidates.csv
    │   └── fuma/
    │       └── FUMA output files (not tracked by Git)
    │
    ├── scripts/
    │   ├── 01_inspect_variants.py
    │   ├── 02_process_fuma_snps.py
    │   ├── 03_analyze_eqtl.py
    │   ├── 04_integrate_eqtl_annotations.py
    │   ├── 05_analyze_chromatin_interaction.py
    │   └── 06_prioritize_variants.py
    │
    ├── results/
    ├── docs/
    ├── .gitignore
    └── README.md

## Reproducibility

The workflow is implemented as sequential Python scripts.

Raw and large FUMA-generated files are excluded from version control,
while processed tables required to document the analytical results are
stored in the repository.

Python scripts use project-relative paths so the workflow is not dependent
on a specific local computer directory.

## Limitations

This workflow is intended for variant prioritization rather than causal
variant identification.

Important limitations include:

- CADD and RegulomeDB provide computational or annotation-based evidence
  and do not establish biological causality.
- Chromatin-interaction evidence operates at the genomic-region level and
  cannot by itself assign a regulatory interaction to an individual SNP.
- eQTL evidence is dependent on the tissues and datasets available.
- Linkage disequilibrium can cause multiple variants to share association
  and regulatory evidence.
- The prioritization thresholds used here are exploratory.

Experimental validation would therefore be required to establish the
functional effects of prioritized variants.

## Current Status

The reproducible computational workflow from initial variant preprocessing
through integrated regulatory prioritization is complete.

Future extensions may include additional functional annotations,
fine-mapping approaches, visualization of prioritized evidence, and
experimental validation strategies.

## Requirements and Usage

The workflow was developed using Python 3 and requires the following
Python package:

- pandas

Install the required dependency with:

    pip install pandas

Run the scripts sequentially from the project root:

    python scripts/01_inspect_variants.py
    python scripts/02_process_fuma_snps.py
    python scripts/03_analyze_eqtl.py
    python scripts/04_integrate_eqtl_annotations.py
    python scripts/05_analyze_chromatin_interaction.py
    python scripts/06_prioritize_variants.py

The FUMA-generated files required by the workflow are not included in
the repository and must be placed in the corresponding `data/fuma/`
directory before running the FUMA-dependent stages.
