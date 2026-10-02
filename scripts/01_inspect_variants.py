import pandas as pd 
df = pd.read_csv("../data/raw/fuma_input_snps.txt", sep="\t")
print(df.head())


import os
os.chdir(r"C:\Users\pc\FTO-variant-prioritization")


import pandas as pd

df = pd.read_csv("data/raw/fuma_input_snps.txt", sep="\t")

print(df.head())

import pandas as pd

# Save the cleaned dataset
df.to_csv(
    "data/processed/fuma_input_snps_clean.csv",
    index=False
)

print("\nCleaned dataset saved successfully.")


import pandas as pd

# Load FUMA input variants
# The original file uses a comma as the decimal separator for P-values.
df = pd.read_csv(
    "data/raw/fuma_input_snps.txt",
    sep="\t",
    decimal=","
)

# Preview the dataset
print("First five variants:")
print(df.head())

# Count the number of variants
print("\nNumber of variants:", len(df))

# Display column names
print("\nColumns:")
print(df.columns.tolist())

# Check for missing values
print("\nMissing values:")
print(df.isnull().sum())

# Check data types
print("\nData types:")
print(df.dtypes)

# Save the cleaned dataset
df.to_csv(
    "data/processed/fuma_input_snps_clean.csv",
    index=False
)

print("\nCleaned dataset saved successfully.")