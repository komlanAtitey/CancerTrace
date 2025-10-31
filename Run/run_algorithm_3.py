#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from cancertrace_algorithm_3 import cancertrace_algorithm_3

# --- Load data ---
df1 = pd.read_csv("data/epithelial.level.time1.csv")
df2 = pd.read_csv("data/epithelial.level.time2.csv")
df3 = pd.read_csv("data/epithelial.level.time3.csv")

# --- Run the algorithm ---
result = cancertrace_algorithm_3(
    data_vector_1 = df1["level_1"],
    data_vector_2 = df2["level_2"],
    data_vector_3 = df3["level_3"],
    gene_vector   = df1["id.time1"],
    driver_genes  = ["TP53INP1", "CCNL2", "VPS37D", "ATP11AUN", "FBXO6"]
)

print("AUC mean:", result["auc_mean"])
print("CIS_matrix shape:", result["CIS_matrix"].shape)
print("Drivers:", list(result["top_influencers"].keys()))

# --- Visualization ---
df = result["knockout_results"]["table"].dropna(subset=["logp_orig", "logp_knock"])
pairs = df["non_driver"] + "→" + df["driver"]

plt.figure(figsize=(12, 5))
plt.bar(pairs, df["logp_orig"], width=0.4, label="Original", color="blue", alpha=0.7)
plt.bar(pairs, df["logp_knock"], width=0.4, label="Knockout", color="red", alpha=0.7, bottom=0)
plt.xticks(rotation=65, ha="right")
plt.ylabel("-log10(p)")
plt.title("Effect of Non-driver Knockout on Granger Causality")
plt.legend()
plt.tight_layout()

# Display the plot
#plt.show()
plt.savefig("knockout_visualization.png", dpi=300)






