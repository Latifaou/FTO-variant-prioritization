
# MotifbreakR Analysis of Six Prioritized FTO Variants

## 1. Objective

To systematically evaluate the predicted effects of six
prioritized noncoding variants at the FTO locus on
transcription-factor binding motifs.

This analysis identifies computational motif changes
and does not establish experimentally validated
transcription-factor binding or regulatory activity.

## 2. Input variants

- rs1421085
- rs16952522
- rs7187250
- rs8063946
- rs4784323
- rs111357538

Reference genome: GRCh37/hg19.

SNP coordinates and alleles were retrieved using
SNPlocs.Hsapiens.dbSNP144.GRCh37 and validated
against BSgenome.Hsapiens.UCSC.hg19.

## 3. Motif database

Human motif models were retrieved from MotifDb.

Selected sources:
- JASPAR 2024
- HOCOMOCO v13

Total motif models: 2,078.

## 4. MotifbreakR parameters

- Method: information content ("ic")
- Initial threshold: 0.50
- filterp: FALSE
- show.neutral: TRUE
- Execution: BiocParallel::SerialParam()

Total motif predictions: 10,655.

## 5. Exploratory filtering

Predictions were retained when:

- motifbreakR effect classification was "strong"
- max(pctRef, pctAlt) >= 0.85

The 0.85 threshold is an exploratory relative
motif-match cutoff, not a statistical significance
threshold.

Total retained predictions: 285.

## 6. Motif-transition classification

Potential creation:
pctRef < 0.85 and pctAlt >= 0.85

Potential disruption:
pctRef >= 0.85 and pctAlt < 0.85

Both alleles match:
pctRef >= 0.85 and pctAlt >= 0.85

Results:
- Potential creations: 201
- Potential disruptions: 69
- Both alleles match: 15

## 7. Reproducibility

The original and regenerated motifbreakR analyses
produced 10,655 predictions, with matching counts
for all six variants.

The complete GRanges objects passed the
all.equal() comparison.

The identical() comparison returned FALSE;
therefore, byte-for-byte object identity is
not claimed.

## 8. Limitations

- Motif predictions do not demonstrate TF occupancy.
- Related TFs may recognize highly similar motifs.
- Multiple models can represent overlapping sequence
  windows rather than independent binding events.
- TF expression and chromatin accessibility were not
  evaluated in this motif-scanning analysis.
- The 0.85 filtering threshold requires sensitivity
  assessment.
- Experimental validation is necessary to establish
  functional regulatory consequences.

## 9. Variant-level motif predictions

Following exploratory filtering, 285 strong motif
predictions were retained across six prioritized
FTO variants.

| Variant | Creation | Disruption | Both match | Total |
|---|---:|---:|---:|---:|
| rs111357538 | 78 | 13 | 4 | 95 |
| rs1421085 | 21 | 10 | 4 | 35 |
| rs16952522 | 53 | 11 | 1 | 65 |
| rs4784323 | 26 | 14 | 3 | 43 |
| rs7187250 | 19 | 7 | 3 | 29 |
| rs8063946 | 4 | 14 | 0 | 18 |
| **Total** | **201** | **69** | **15** | **285** |

These counts represent motif-model predictions,
not independent transcription-factor binding sites.

## 10. Exploratory TF-family observations

### rs111357538

The 78 predicted motif creations corresponded to
70 distinct TF labels and 10 distinct reported
sequence strings.

The highest-scoring models included multiple
homeodomain-associated TFs, particularly HOX
family members.

The overlapping sequence windows suggest motif
redundancy rather than 78 independent regulatory
events.

### rs16952522

Eleven predicted motif disruptions involved:

- RUNX family: 4 models
- T-box family: 5 models
- ZEB2: 1 model
- MESP1: 1 model

Normalized motif-score differences (ALT - REF)
ranged from approximately -0.202 to -0.073.

The tested motif models and filtering criteria
did not provide convincing support for strong
CTCF motif disruption by rs16952522.

Alternative TF-binding hypotheses remain
computational and require independent evidence.

### rs8063946

Fourteen predicted motif disruptions were retained.

Seven models corresponded to bZIP-related TFs:
BACH1, BATF3, FOSL1, FOSL2, JDP2, JUND and NFE2.

These models showed reduced alternative-allele
motif-match scores and overlapping sequence windows.

The seven TF labels were not observed among the
retained disruption predictions for the other
five variants under the same filtering criteria.

This observation does not establish altered
TF occupancy or regulatory activity.

### Other variants

rs1421085, rs4784323 and rs7187250 also showed
predicted motif changes involving several
developmental transcription-factor families.

These observations are descriptive and have not
been independently validated.

## 11. Interpretation

The six variants exhibited distinct computational
motif-change profiles.

However, motif counts, relative scores and
TF-label diversity cannot independently establish
variant functionality or biological importance.

Further interpretation requires integration with
regulatory annotations, chromatin context,
TF expression and relevant experimental evidence.

No causal regulatory mechanism is established
by this motifbreakR analysis alone.

## 12. Threshold-sensitivity analysis

To evaluate the dependence of motif predictions on
the exploratory relative-score cutoff, the filtering
and motif-transition classification were repeated
using thresholds of 0.80, 0.85 and 0.90.

Only predictions classified as "strong" by
motifbreakR were considered.

### 12.1. Threshold-dependent predictions

| Threshold | Creation | Disruption | Both match | Total |
|---|---:|---:|---:|---:|
| 0.80 | 356 | 146 | 41 | 543 |
| 0.85 | 201 | 69 | 15 | 285 |
| 0.90 | 110 | 34 | 1 | 145 |


The number and classification of retained predictions
varied substantially with the selected cutoff.

These differences demonstrate sensitivity to the
filtering criterion and do not independently establish
biological relevance.

### 12.2. rs16952522

At the 0.90 cutoff, four motif models were classified
as potential disruptions:

- MESP1.H13CORE.0.S.C (HOCOMOCO v13)
- MA0684.3 / RUNX3 (JASPAR 2024)
- TBX2.H13CORE.0.SM.B (HOCOMOCO v13)
- TBX3.H13CORE.0.PSM.A (HOCOMOCO v13)

The MESP1 and RUNX3 MA0684.3 models satisfied the
disruption criteria at all three tested thresholds.

TBX3 satisfied the criteria at 0.85 and 0.90.
TBX2 satisfied them only at 0.90.

Multiple motif models involved overlapping sequence
contexts. They must not be interpreted as independent
experimental evidence of TF binding.

### 12.3. rs8063946

At the 0.85 cutoff, 14 potential disruptions were
identified, including seven bZIP-related TF models.

At the 0.90 cutoff, only one potential disruption
remained: a SOX2 motif model.

This indicates that the bZIP-related disruption
classifications were sensitive to the cutoff.

### 12.4. Interpretation

Threshold stability reflects consistency under
the tested computational filtering criteria.

It does not establish transcription-factor occupancy,
chromatin accessibility, regulatory function or
causal effects on gene expression.

The complete threshold-sensitivity counts are
available in:

data/processed/six_variants_threshold_sensitivity.csv
