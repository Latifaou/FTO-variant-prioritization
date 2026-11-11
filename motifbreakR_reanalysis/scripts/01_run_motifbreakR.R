
# ============================================================
# FTO locus: motifbreakR motif scanning
# Six prioritized noncoding variants
# Reference genome: GRCh37 / hg19
# ============================================================

# 01. Load required packages
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BiocParallel)

# 02. Define the six prioritized variants
priority_rsids <- c(
  "rs1421085",
  "rs16952522",
  "rs7187250",
  "rs8063946",
  "rs4784323",
  "rs111357538"
)

# 03. Retrieve SNP coordinates and alleles
all_six_snps <- motifbreakR::snps.from.rsid(
  rsid = priority_rsids,
  dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37,
  search.genome = BSgenome.Hsapiens.UCSC.hg19
)

# 04. Retrieve human motif models
human_motifs <- MotifDb::query(
  MotifDb::MotifDb,
  "Hsapiens"
)

# 05. Select JASPAR 2024 and HOCOMOCO v13
primary_motifs <- human_motifs[
  S4Vectors::mcols(human_motifs)$dataSource %in%
    c("jaspar2024", "HOCOMOCOv13")
]

# 06. Verify inputs
cat("Variants:", length(all_six_snps), "\n")
cat("Motif models:", length(primary_motifs), "\n")

# ============================================================
# 07. Run motifbreakR
# ============================================================

all_six_motifs <- motifbreakR::motifbreakR(
  snpList = all_six_snps,
  pwmList = primary_motifs,
  threshold = 0.50,
  filterp = FALSE,
  show.neutral = TRUE,
  method = "ic",
  BPPARAM = BiocParallel::SerialParam()
)

# 08. Verify result dimensions
cat(
  "Total motif predictions:",
  length(all_six_motifs),
  "\n"
)

cat(
  "Variants represented:",
  length(unique(all_six_motifs$SNP_id)),
  "\n"
)

# ============================================================
# 09. Save regenerated motifbreakR results
# ============================================================

saveRDS(
  all_six_motifs,
  file = "data/processed/all_six_motifbreakR_reproduced.rds"
)

# Verify successful saving
file.exists(
  "data/processed/all_six_motifbreakR_reproduced.rds"
)


# ============================================================
# 10. Verify reproducibility against original results
# ============================================================

original_results <- readRDS(
  "data/processed/all_six_motifbreakR_results.rds"
)

reproduced_results <- readRDS(
  "data/processed/all_six_motifbreakR_reproduced.rds"
)

# Verify scientific equivalence
stopifnot(
  isTRUE(
    all.equal(original_results, reproduced_results)
  )
)

cat("MotifbreakR reproducibility verification passed.\n")
