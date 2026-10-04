from pathlib import Path
import pandas as pd


# ============================================================
# 1. Project paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FUMA_DIR = PROJECT_ROOT / "data" / "fuma" / "FUMA_job721962"


# ============================================================
# 2. Load processed datasets
# ============================================================

snps = pd.read_csv(
    PROCESSED_DIR / "fuma_snps_clean.csv"
)

eqtl = pd.read_csv(
    PROCESSED_DIR / "fuma_eqtl_clean.csv"
)

ci = pd.read_csv(
    PROCESSED_DIR / "fuma_chromatin_interactions.csv"
)

annot = pd.read_csv(
    FUMA_DIR / "annot.txt",
    sep="\t"
)


# ============================================================
# 3. Build the base variant evidence table
# ============================================================

priority_columns = [
    "uniqID",
    "rsID",
    "chr",
    "pos",
    "effect_allele",
    "MAF",
    "gwasP",
    "r2",
    "IndSigSNP",
    "nearestGene",
    "CADD",
    "RDB",
    "minChrState",
    "commonChrState",
    "eqtlMapFilt",
    "ciMapFilt"
]

priority = snps[priority_columns].copy()


# ============================================================
# 4. Summarize eQTL evidence
# ============================================================

eqtl_summary = (
    eqtl.groupby("uniqID")
    .agg(
        eqtl_genes=(
            "symbol",
            lambda x: ";".join(sorted(x.dropna().unique()))
        ),
        eqtl_tissues=(
            "tissue",
            lambda x: ";".join(sorted(x.dropna().unique()))
        ),
        eqtl_min_p=("p", "min"),
        eqtl_min_FDR=("FDR", "min"),
        eqtl_n_associations=("symbol", "count")
    )
    .reset_index()
)

priority = pd.merge(
    priority,
    eqtl_summary,
    on="uniqID",
    how="left"
)


# ============================================================
# 5. Summarize chromatin-interaction evidence
# ============================================================

# FUMA chromatin-interaction rows can contain multiple SNPs.
# Expand the SNP lists so each row represents one SNP-region
# interaction before summarizing the evidence per variant.

ci_snp = ci.copy()

ci_snp["SNPs"] = ci_snp["SNPs"].str.split(";")

ci_snp = ci_snp.explode("SNPs")

ci_snp = ci_snp.rename(
    columns={"SNPs": "rsID"}
)

ci_summary = (
    ci_snp.groupby("rsID")
    .agg(
        ci_genes=(
            "symbol",
            lambda x: ";".join(sorted(x.dropna().unique()))
        ),
        ci_tissues=(
            "tissue/cell",
            lambda x: ";".join(sorted(x.dropna().unique()))
        ),
        ci_min_FDR=("FDR", "min"),
        ci_n_interactions=("symbol", "count"),
        ci_n_genes=("symbol", "nunique"),
        ci_n_tissues=("tissue/cell", "nunique")
    )
    .reset_index()
)

priority = pd.merge(
    priority,
    ci_summary,
    on="rsID",
    how="left"
)


# ============================================================
# 6. Add GWAS and chromatin-interaction indicators
# ============================================================

# True when a direct GWAS P-value is available for the variant
# in the FUMA SNP table.
priority["has_direct_gwas_p"] = priority["gwasP"].notna()

# True when the variant itself is the independent significant
# SNP representing its FUMA association signal.
priority["is_independent_sig"] = (
    priority["rsID"] == priority["IndSigSNP"]
)

# Chromatin-interaction evidence is regional rather than proof
# of a direct SNP-to-gene regulatory relationship.
priority["ci_IRX3"] = (
    priority["ci_genes"]
    .fillna("")
    .str.split(";")
    .apply(lambda genes: "IRX3" in genes)
)

priority["ci_IRX5"] = (
    priority["ci_genes"]
    .fillna("")
    .str.split(";")
    .apply(lambda genes: "IRX5" in genes)
)


# ============================================================
# 7. Select variants with convergent regulatory evidence
# ============================================================

# Exploratory prioritization criterion:
#   CADD >= 15
#   AND
#   RegulomeDB category 2a, 2b, or 3a
#
# These criteria identify candidates supported by two
# complementary regulatory annotation approaches. They do not
# establish causality.

convergent_regulatory = priority[
    (priority["CADD"] >= 15)
    & (priority["RDB"].isin(["2a", "2b", "3a"]))
].copy()


# ============================================================
# 8. Extract adipose-relevant ChromHMM annotations
# ============================================================

# Roadmap Epigenomics samples:
# E023 = mesenchymal stem cell-derived adipocyte cultured cells
# E025 = adipose-derived mesenchymal stem cell cultured cells
# E063 = adipose nuclei

adipose_states = annot[
    ["uniqID", "E023", "E025", "E063"]
].copy()

adipose_states = pd.merge(
    priority[["uniqID", "rsID"]],
    adipose_states,
    on="uniqID",
    how="left"
)


# Roadmap Epigenomics core 15-state ChromHMM labels
chromhmm_labels = {
    1: "Active TSS",
    2: "Flanking Active TSS",
    3: "Transcription at gene 5'/3'",
    4: "Strong transcription",
    5: "Weak transcription",
    6: "Genic enhancer",
    7: "Enhancer",
    8: "ZNF genes & repeats",
    9: "Heterochromatin",
    10: "Bivalent/Poised TSS",
    11: "Flanking Bivalent TSS/Enhancer",
    12: "Bivalent Enhancer",
    13: "Repressed PolyComb",
    14: "Weak Repressed PolyComb",
    15: "Quiescent/Low"
}

for column in ["E023", "E025", "E063"]:
    adipose_states[f"{column}_label"] = (
        adipose_states[column].map(chromhmm_labels)
    )


# ============================================================
# 9. Build the final candidate table
# ============================================================

adipose_for_merge = adipose_states[
    [
        "rsID",
        "E023_label",
        "E025_label",
        "E063_label"
    ]
].copy()

final_candidates = pd.merge(
    convergent_regulatory,
    adipose_for_merge,
    on="rsID",
    how="left"
)

final_columns = [
    "rsID",
    "uniqID",
    "gwasP",
    "r2",
    "IndSigSNP",
    "is_independent_sig",
    "CADD",
    "RDB",
    "eqtl_genes",
    "eqtl_tissues",
    "eqtl_min_FDR",
    "ci_IRX3",
    "ci_IRX5",
    "E023_label",
    "E025_label",
    "E063_label"
]

final_candidates = final_candidates[final_columns]


# ============================================================
# 10. Save results
# ============================================================

output_file = (
    PROCESSED_DIR / "prioritized_candidates.csv"
)

final_candidates.to_csv(
    output_file,
    index=False
)


# ============================================================
# 11. Report workflow summary
# ============================================================

print("\nVariant prioritization completed successfully.")

print(f"\nTotal FUMA variants: {len(priority)}")

print(
    "Independent significant SNPs:",
    priority["is_independent_sig"].sum()
)

print(
    "Variants with direct GWAS P-values:",
    priority["has_direct_gwas_p"].sum()
)

print(
    "Variants with eQTL evidence:",
    priority["eqtl_genes"].notna().sum()
)

print(
    "Convergent regulatory candidates:",
    len(final_candidates)
)

print("\nFinal prioritized candidate table:")
print(final_candidates.to_string(index=False))

print("\nResults saved to:")
print(output_file)