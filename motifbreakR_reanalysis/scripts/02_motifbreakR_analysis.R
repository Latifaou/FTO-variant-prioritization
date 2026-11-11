
# ============================================================
# FTO locus: motifbreakR analysis of six prioritized variants
# Genome assembly: GRCh37 / hg19
# ============================================================

# 01. Load required packages
library(dplyr)
library(readr)
library(motifbreakR)
library(S4Vectors)

# 02. Define project paths
raw_dir <- "data/raw"
processed_dir <- "data/processed"
reports_dir <- "reports"

# 03. Load validated motifbreakR results
all_six_motifs <- readRDS(
  file.path(processed_dir,
            "all_six_motifbreakR_results.rds")
)

# 04. Convert motif metadata to a data frame
six_motif_table <- as.data.frame(
  S4Vectors::mcols(all_six_motifs)
)

# 05. Verify input dimensions
cat("Total motif predictions:",
    nrow(six_motif_table), "\n")

cat("Number of variants:",
    dplyr::n_distinct(six_motif_table$SNP_id),
    "\n")

# ============================================================
# 06. Filter motif predictions
# ============================================================

motif_threshold <- 0.85

six_filtered <- six_motif_table |>
  dplyr::filter(
    effect == "strong",
    pmax(pctRef, pctAlt) >= motif_threshold
  ) |>
  dplyr::mutate(
    motif_transition = dplyr::case_when(
      pctRef < motif_threshold &
        pctAlt >= motif_threshold ~ "Potential creation",
      
      pctRef >= motif_threshold &
        pctAlt < motif_threshold ~ "Potential disruption",
      
      pctRef >= motif_threshold &
        pctAlt >= motif_threshold ~ "Both alleles match"
    ),
    delta_pct = pctAlt - pctRef
  )

# 07. Verify filtered results
cat("Retained motif predictions:",
    nrow(six_filtered), "\n")

print(
  table(
    six_filtered$SNP_id,
    six_filtered$motif_transition
  )
)

# ============================================================
# 08. Export filtered predictions and summary
# ============================================================

# Create output directory if necessary
dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Export all 285 filtered predictions
readr::write_csv(
  six_filtered,
  file.path(
    processed_dir,
    "six_variants_filtered_motifs_085.csv"
  )
)

# Generate summary by variant and transition
motif_summary <- six_filtered |>
  dplyr::group_by(SNP_id, motif_transition) |>
  dplyr::summarise(
    motif_models = dplyr::n(),
    unique_tf_labels = dplyr::n_distinct(geneSymbol),
    unique_seq_matches = dplyr::n_distinct(seqMatch),
    .groups = "drop"
  )

# Export summary
readr::write_csv(
  motif_summary,
  file.path(
    processed_dir,
    "six_variants_motif_transition_summary.csv"
  )
)

# Verify exported files
file.exists(
  file.path(
    processed_dir,
    "six_variants_filtered_motifs_085.csv"
  )
)

file.exists(
  file.path(
    processed_dir,
    "six_variants_motif_transition_summary.csv"
  )
)

# ============================================================
# 09. Automated validation checks
# ============================================================

# Verify original and filtered prediction counts
stopifnot(
  nrow(six_motif_table) == 10655,
  nrow(six_filtered) == 285,
  dplyr::n_distinct(six_filtered$SNP_id) == 6
)

# Verify transition classifications
stopifnot(
  sum(six_filtered$motif_transition ==
        "Potential creation") == 201,
  
  sum(six_filtered$motif_transition ==
        "Potential disruption") == 69,
  
  sum(six_filtered$motif_transition ==
        "Both alleles match") == 15
)

# Verify that every retained prediction meets the filter
stopifnot(
  all(six_filtered$effect == "strong"),
  all(
    pmax(six_filtered$pctRef,
         six_filtered$pctAlt) >= motif_threshold
  ),
  !anyNA(six_filtered$motif_transition)
)

cat("All motifbreakR validation checks passed.\n")


# ============================================================
# Threshold-sensitivity analysis
# ============================================================

sensitivity_thresholds <- c(0.80, 0.85, 0.90)

threshold_sensitivity <- lapply(
  sensitivity_thresholds,
  function(cutoff) {
    
    six_motif_table |>
      dplyr::filter(
        effect == "strong",
        pmax(pctRef, pctAlt) >= cutoff
      ) |>
      dplyr::mutate(
        threshold = cutoff,
        transition = dplyr::case_when(
          pctRef < cutoff & pctAlt >= cutoff ~ "Creation",
          pctRef >= cutoff & pctAlt < cutoff ~ "Disruption",
          pctRef >= cutoff & pctAlt >= cutoff ~ "Both match"
        )
      ) |>
      dplyr::count(threshold, SNP_id, transition)
  }
) |>
  dplyr::bind_rows()

# Export results
readr::write_csv(
  threshold_sensitivity,
  "data/processed/six_variants_threshold_sensitivity.csv"
)

# Validate expected totals
sensitivity_totals <- threshold_sensitivity |>
  dplyr::group_by(threshold) |>
  dplyr::summarise(
    total = sum(n),
    .groups = "drop"
  )

stopifnot(
  identical(
    sensitivity_totals$total,
    c(543L, 285L, 145L)
  )
)

cat("Threshold-sensitivity analysis verified.\n")

