from pathlib import Path 
import pandas as pd

# Determine the project root from the location of this script
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Load FUMA chromatin interaction results
ci = pd.read_csv(
    PROJECT_ROOT / "data" / "fuma" / "FUMA_job721962" / "ci.txt",
    sep="\t"
)

# Inspect the dataset
print("Chromatin interaction table dimensions:")
print(ci.shape)

print("\nColumns:")
print(ci.columns.tolist())

print("\nFirst five interactions:")
print(ci.head())

# Count missing values
print("\nMissing values:")
print(ci.isnull().sum())

# Count chromatin interactions by tissue/cell type
print("\nInteractions by tissue/cell type:")
print(ci["tissue/cell"].value_counts())

# Count mapped vs unmapped interactions
print("\nChromatin interaction mapping filter:")
print(ci["ciMapFilt"].value_counts())

# Count unique target gene IDs
print("\nNumber of unique mapped gene IDs:")
print(ci["genes"].dropna().nunique())

# Load FUMA gene mapping results
genes = pd.read_csv(
    PROJECT_ROOT / "data" / "fuma" / "FUMA_job721962" / "genes.txt",
    sep="\t"
)

print("\nGene table dimensions:")
print(genes.shape)

print("\nGene table columns:")
print(genes.columns.tolist())

print("\nFirst five genes:")
print(genes.head())

# Keep only the gene information needed for annotation
gene_lookup = genes[
    ["ensg", "symbol"]
]

# Add gene symbols to the chromatin interaction table
ci_annotated = pd.merge(
    ci,
    gene_lookup,
    left_on="genes",
    right_on="ensg",
    how="left"
)
print("\nAnnotated chromatin interaction table:")
print(
    ci_annotated[
        [
            "region1",
            "region2",
            "tissue/cell",
            "SNPs",
            "genes",
            "symbol",
            "FDR",
            "ciMapFilt"
        ]
    ].head(10)
)

# Keep only chromatin interactions mapped to genes
mapped_ci = ci_annotated[
    ci_annotated["ciMapFilt"] == 1
].copy()

print("\nMapped chromatin interactions:")
print(mapped_ci.shape)

# Count interactions for each mapped gene
print("\nChromatin interactions by target gene:")
print(
    mapped_ci["symbol"]
    .value_counts()
)

# Select chromatin interactions mapped to IRX3 or IRX5
irx_interactions = mapped_ci[
    mapped_ci["symbol"].isin(["IRX3", "IRX5"])
].copy()

print("\nIRX3 and IRX5 chromatin interactions:")
print(irx_interactions.shape)

print("\nIRX3/IRX5 interactions by gene:")
print(irx_interactions["symbol"].value_counts())

print("\nIRX3/IRX5 interactions by tissue/cell type:")
print(
    irx_interactions.groupby(
        ["symbol", "tissue/cell"]
    ).size()
)

# Expand the semicolon-separated SNP lists
irx_snps = irx_interactions.copy()

irx_snps["SNPs"] = irx_snps["SNPs"].str.split(";")

irx_snps = irx_snps.explode("SNPs")

print("\nExpanded IRX3/IRX5 SNP-interaction table:")
print(irx_snps.shape)

print("\nNumber of unique SNPs linked to IRX3/IRX5 interactions:")
print(irx_snps["SNPs"].nunique())

# Select IRX3/IRX5 interactions in Mesenchymal Stem Cells
msc_irx = irx_snps[
    irx_snps["tissue/cell"] == "Mesenchymal_Stem_Cell"
].copy()

print("\nIRX3/IRX5 interactions in Mesenchymal Stem Cells:")
print(msc_irx.shape)

print("\nUnique SNPs in MSC interactions:")
print(msc_irx["SNPs"].nunique())

print("\nMSC interactions by target gene:")
print(msc_irx["symbol"].value_counts())

candidate_snps = [
    "rs1421085",
    "rs16952522",
    "rs7187250",
    "rs8063946",
    "rs4784323",
    "rs111357538"
]

candidate_msc = msc_irx[
    msc_irx["SNPs"].isin(candidate_snps)
]

print("\nCandidate variants in MSC IRX3/IRX5 interactions:")
print(
    candidate_msc[
        ["SNPs", "symbol", "region1", "region2", "FDR"]
    ].to_string(index=False)
)

# Save mapped chromatin interactions with gene annotations
mapped_ci.to_csv(
    PROJECT_ROOT / "data" / "processed" / "fuma_chromatin_interactions.csv",
    index=False
)

print("\nProcessed chromatin interaction table saved successfully.")