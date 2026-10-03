import pandas as pd

# Load the FUMA SNP annotation table
df = pd.read_csv(
    "data/fuma/FUMA_job721962/snps.txt",
    sep="\t"
)

print(df.head())
# Display dataset dimensions
print("\nDataset dimensions:")
print(df.shape)

# Display column names
print("\nColumns:")
print(df.columns.tolist())

# Check data types
print("\nData types:")
print(df.dtypes)

# Check missing values
print("\nMissing values:")
print(df.isnull().sum())

# Inspect functional annotation categories
print("\nFunctional categories:")
print(df["func"].value_counts())

# Inspect RegulomeDB categories
print("\nRegulomeDB categories:")
print(df["RDB"].value_counts(dropna=False))

# inspect mapping filters
print("\nPositional mapping:")
print(df["posMapFilt"].value_counts())

print("\neQTL mapping:")
print(df["eqtlMapFilt"].value_counts())

print("\nChromatin interaction mapping:")
print(df["ciMapFilt"].value_counts())

# Identify variants passing the eQTL mapping filter
eqtl_variants = df[df["eqtlMapFilt"] == 1]

print("\nVariants passing the eQTL mapping filter:")
print(
    eqtl_variants[
        ["rsID", "nearestGene", "CADD", "RDB", "r2"]
    ]
)

# Save the processed FUMA SNP table
df.to_csv(
    "data/processed/fuma_snps_clean.csv",
    index=False
)

print("\nProcessed FUMA SNP table saved successfully.")