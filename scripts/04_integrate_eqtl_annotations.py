from pathlib import Path
import pandas as pd

# Determine the project root from the location of this script
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Load the processed FUMA SNP annotations
snps = pd.read_csv(
    PROJECT_ROOT / "data" / "processed" / "fuma_snps_clean.csv"
)

# Load the processed eQTL results
eqtl = pd.read_csv(
    PROJECT_ROOT / "data" / "processed" / "fuma_eqtl_clean.csv"
)

print("SNP table:", snps.shape)
print("eQTL table:", eqtl.shape)

# Merge eQTL results with SNP annotations
merged = pd.merge(
    eqtl,
    snps,
    on="uniqID",
    how="left",
    suffixes=("_eqtl", "_snp")
)

print("\nMerged table dimensions:")
print(merged.shape)

# Select key columns from the integrated dataset
key_columns = [
    "rsID",
    "uniqID",
    "tissue",
    "symbol",
    "testedAllele",
    "p",
    "signed_stats",
    "FDR",
    "MAF",
    "r2",
    "CADD",
    "RDB"
]

integrated = merged[key_columns]

print("\nIntegrated eQTL annotation table:")
print(integrated.to_string(index=False))

# Save the integrated eQTL annotation table
integrated.to_csv(
    PROJECT_ROOT / "data" / "processed" / "integrated_eqtl_annotations.csv",
    index=False
)

print("\nIntegrated eQTL annotation table saved successfully.")