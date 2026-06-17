import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Resolve the CSV one folder up from this script (data/reference/integrated_i_he_abundance.csv)
csv_path = Path(__file__).resolve().parent.parent / "integrated_i_he_abundance.csv"
if not csv_path.exists():
	raise FileNotFoundError(f"Expected CSV at {csv_path!s} not found. Put 'integrated_i_he_abundance.csv' in the parent folder of this script.")

df = pd.read_csv(csv_path)

fig, ax = plt.subplots(figsize=(8, 4))

ax.bar(df["n"], df["ionPercent"])

ax.set_ylabel("plain ion abundance / %")
ax.set_xlabel("cluster")
ax.set_title("Integrated I$^+$He$_n$ abundance")

ax.set_xlim(-0.5, df["n"].max() + 0.5)
ax.set_xticks(df["n"])
ax.set_xticklabels(df["label"], rotation=45, ha="right")

fig.tight_layout()
plt.show()