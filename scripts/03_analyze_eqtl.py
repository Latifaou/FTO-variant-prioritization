import pandas as pd

# Load FUMA eQTL mapping results
eqtl = pd.read_csv(
    "data/fuma/FUMA_job721962/eqtl.txt",
    sep="\t"
)

# Preview the dataset
print("First five eQTL associations:")
print(eqtl.head())

# Display dataset dimensions
print("\nDataset dimensions:")
print(eqtl.shape)

# Display column names
print("\nColumns:")
print(eqtl.columns.tolist())

# Check data types
print("\nData types:")
print(eqtl.dtypes)

# Check missing values
print("\nMissing values:")
print(eqtl.isnull().sum())

# Display the key eQTL associations
print("\nKey eQTL associations:")

print(
    eqtl[
        [
            "uniqID",
            "tissue",
            "symbol",
            "testedAllele",
            "p",
            "signed_stats",
            "FDR"
        ]
    ].to_string(index=False)
)

# Count eQTL associations by target gene
print("\nAssociations by target gene:")
print(eqtl["symbol"].value_counts())

# Count eQTL associations by tissue
print("\nAssociations by tissue:")
print(eqtl["tissue"].value_counts())

# Count unique variants in the eQTL results
print("\nNumber of unique eQTL variants:")
print(eqtl["uniqID"].nunique())

# Count associations per variant
print("\nAssociations per variant:")
print(eqtl["uniqID"].value_counts())

# Remove columns that contain only missing values
eqtl_clean = eqtl.dropna(axis=1, how="all")

# Save the cleaned eQTL table
eqtl_clean.to_csv(
    "data/processed/fuma_eqtl_clean.csv",
    index=False
)

print("\nCleaned eQTL table saved successfully.")
print("Cleaned dimensions:", eqtl_clean.shape)