import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Resolve the CSV one folder up from this script (data/reference/integrated_i_he_abundance.csv)
csv_path = Path(__file__).resolve().parent.parent / "integrated_i_he_abundance.csv"
if not csv_path.exists():
	raise FileNotFoundError(f"Expected CSV at {csv_path!s} not found. Put 'integrated_i_he_abundance.csv' in the parent folder of this script.")

df = pd.read_csv(csv_path)

# Original (simple) histogram
fig, ax = plt.subplots(figsize=(8, 4))

ax.bar(df["n"], df["ionPercent"])

ax.set_ylabel("plain ion abundance / %")
ax.set_xlabel("cluster")
ax.set_title("Integrated I$^+$He$_n$ abundance")

ax.set_xlim(-0.5, 20 + 0.5)
ax.set_xticks(list(range(0, 21)))
ax.set_xticklabels([str(n) for n in range(0, 21)], rotation=45, ha="right")

fig.tight_layout()

# Second, nicer histogram per request, excluding n==0 and reweighing to 100%
df2 = df[df["n"] != 0].copy()
if df2.empty:
	raise ValueError("No entries with n != 0 found in the CSV; cannot create reweighted histogram.")

# compute reweighted fraction so the remaining bins sum to 100%
total = df2["ionPercent"].sum()
if total <= 0:
	raise ValueError("Sum of ionPercent for n!=0 is zero or negative; cannot reweight.")
df2["fractionPercent"] = df2["ionPercent"] / total * 100.0

fig2, ax2 = plt.subplots(figsize=(8, 4))

# Use a carmine/crimson-like color for the bars
carmine_red = "#960018"  # deep carmine red
ax2.bar(df2["n"], df2["fractionPercent"], color=carmine_red, edgecolor="black")
# ax2.set_yscale('log')
# ax2.set_ylim(1e-1)

ax2.set_ylabel("fraction / %")
ax2.set_xlabel(r"n in I$^+$He$_{\mathrm{n}}$")
#ax2.set_title("Integrated I$^+He_n$ abundance (n>0, reweighted to 100%)")

# nicer grid with low alpha
ax2.grid(True, alpha=0.3, linestyle="--")

# force x-ticks from 0..20 as requested
ax2.set_xlim(0.5, 20 + 0.5)
ax2.set_xticks(list(range(1, 21)))
ax2.set_xticklabels([str(n) for n in range(1, 21)])

fig2.tight_layout()

plt.show()


