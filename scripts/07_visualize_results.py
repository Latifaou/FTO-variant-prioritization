from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


# Determine project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent


# Load prioritized candidate variants
candidates = pd.read_csv(
    PROJECT_ROOT / "data" / "processed" / "prioritized_candidates.csv"
)


# Inspect the dataset
print("Candidate table dimensions:")
print(candidates.shape)

print("\nCandidate variants:")
print(candidates[["rsID", "CADD", "RDB"]])

# Sort variants by CADD score for visualization
plot_data = candidates.sort_values(
    "CADD",
    ascending=True
)

# Create horizontal bar plot
plt.figure(figsize=(8, 5))

plt.barh(
    plot_data["rsID"],
    plot_data["CADD"]
)

# Add exploratory CADD threshold
plt.axvline(
    x=15,
    linestyle="--",
    label="Exploratory CADD threshold (15)"
)

# Add labels and title
plt.xlabel("CADD score")
plt.ylabel("Variant")
plt.title("CADD Scores of Prioritized FTO Variants")

# Add CADD values next to each bar
for index, row in plot_data.iterrows():
    plt.text(
        row["CADD"] + 0.15,
        index,
        f'{row["CADD"]:.2f}',
        va="center"
    )

# Add some space for the value labels
plt.xlim(0, 19)

plt.legend()
plt.tight_layout()


# Create results directory if it does not already exist
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)


# Save figure
plt.savefig(
    RESULTS_DIR / "cadd_prioritized_variants.png",
    dpi=300,
    bbox_inches="tight"
)

# Display figure
plt.show()

# ============================================================
# Figure 2: Integrated evidence matrix
# ============================================================

evidence_matrix = candidates[
    [
        "rsID",
        "is_independent_sig",
        "RDB",
        "ci_IRX3",
        "ci_IRX5",
        "E023_label",
        "E025_label",
        "E063_label"
    ]
].copy()

print("\nEvidence matrix:")
print(evidence_matrix.to_string(index=False))

# Prepare categorical evidence for visualization
display_matrix = evidence_matrix[
    [
        "rsID",
        "is_independent_sig",
        "RDB",
        "E023_label",
        "E025_label",
        "E063_label"
    ]
].copy()

# Make the independent-SNP indicator easier to read
display_matrix["is_independent_sig"] = (
    display_matrix["is_independent_sig"]
    .map({
        True: "Yes",
        False: "No"
    })
)

# Rename columns for the figure
display_matrix = display_matrix.rename(
    columns={
        "is_independent_sig": "Independent\nGWAS signal",
        "RDB": "RegulomeDB",
        "E023_label": "E023\nMSC-derived adipocytes",
        "E025_label": "E025\nAdipose-derived MSC",
        "E063_label": "E063\nAdipose nuclei"
    }
)

print("\nDisplay matrix:")
print(display_matrix.to_string(index=False))

# ============================================================
# Create integrated evidence matrix figure
# ============================================================

# Use rsID as row labels
figure_data = display_matrix.set_index("rsID")

# Create figure
fig, ax = plt.subplots(figsize=(12, 3))

# Hide normal plot axes
ax.axis("off")

# Create table
table = ax.table(
    cellText=figure_data.values,
    rowLabels=figure_data.index,
    colLabels=figure_data.columns,
    cellLoc="center",
    rowLoc="center",
    loc="center"
)

# Adjust table appearance
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 2.2)

# Make column headers bold
for col in range(len(figure_data.columns)):
    table[(0, col)].set_text_props(weight="bold")

# Add title
plt.title(
    "Integrated Regulatory Evidence for Prioritized FTO Variants",
    fontsize=12,
    fontweight="bold",
    pad=15
)

# Save figure
plt.savefig(
    RESULTS_DIR / "integrated_evidence_matrix.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()