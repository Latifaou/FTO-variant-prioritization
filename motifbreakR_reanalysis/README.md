
# FTO Noncoding Variant Reanalysis

## Overview

This repository documents a computational reanalysis
of noncoding genetic variants at the FTO locus,
associated with body mass index and obesity.

The project aims to provide a transparent,
reproducible workflow for variant prioritization,
transcription-factor binding motif analysis,
and regulatory annotation.

## Research objectives

1. Prioritize candidate noncoding FTO variants.
2. Evaluate allele-dependent TF motif-score changes.
3. Integrate variant-to-gene mapping evidence.
4. Interpret functional annotation results.

## Prioritized variants

- rs1421085
- rs16952522
- rs7187250
- rs8063946
- rs4784323
- rs111357538

Reference genome: GRCh37/hg19.

## Motif analysis

MotifbreakR was used with 2,078 human motif models
from JASPAR 2024 and HOCOMOCO v13.

The analysis produced 10,655 motif predictions
across the six variants.

Exploratory filtering retained strong-effect
predictions with a maximum relative motif-match
score of at least 0.85.

Threshold sensitivity was evaluated at 0.80,
0.85 and 0.90.

Detailed methods, results and limitations:

reports/motifbreakR_methods_results.md

## Reproducibility

The motif-scanning workflow can be rerun using:

scripts/01_run_motifbreakR.R

Filtering and summary generation:

scripts/02_motifbreakR_analysis.R

The regenerated motifbreakR results passed
the all.equal() comparison against the original
results. Byte-for-byte identity is not claimed.

The R session information is available in:

reports/motifbreakR_session_info.txt

## Project status

- FUMA data collection: completed
- Multicriteria variant prioritization: completed
- MotifbreakR analysis: completed
- Threshold-sensitivity analysis: completed
- Variant-to-gene mapping interpretation: pending
- FUMA GENE2FUNC interpretation: pending

The upstream FUMA prioritization workflow and its
reproducibility documentation are maintained
separately and are not fully reconstructed by
the two motifbreakR scripts in this repository.

## Scientific limitations

Predicted motif changes do not demonstrate
transcription-factor binding, regulatory activity,
or causal effects on gene expression.

Motif models can be redundant and may represent
overlapping sequence contexts.

Computational predictions require independent
functional evidence and experimental validation.
