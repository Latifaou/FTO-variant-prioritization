import pandas as pd 
df = pd.read_csv("../data/raw/fuma_input_snps.txt", sep="\t")
print(df.head())


import os
os.chdir(r"C:\Users\pc\FTO-variant-prioritization")


import pandas as pd

df = pd.read_csv("data/raw/fuma_input_snps.txt", sep="\t")

print(df.head())