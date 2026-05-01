import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def verify_vmi_data(data_dir: str = "data/reference"):
    """Loads and plots the exported VMI reference data to ensure parity with MATLAB."""
    base_path = Path(data_dir)

    # 1. Load the pre-processed data
    try:
        df_he = pd.read_csv(base_path / "vmi_iplus_he.csv")
        df_gas = pd.read_csv(base_path / "vmi_iplus_gas.csv")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Make sure you moved the CSVs into the data/reference/ folder.")
        return

    # 2. Recreate the MATLAB plot
    fig, ax = plt.subplots(figsize=(8, 5))

    # Gas phase normalization (MATLAB: b_v = res_Iplus_gas.r*vf_single/100>4)
    mask_gas = df_gas['v_Aps'] > 4.0
    max_gas_signal = df_gas.loc[mask_gas, 'signal_arb'].max()

    ax.plot(
        df_gas['v_Aps'],
        df_gas['signal_arb'] / max_gas_signal,
        label="I$_2$:I$^+$ (Gas Phase)",
        color="#2c7fb8",  # Using a nice blue
        linewidth=2
    )

    # He Droplet normalization
    max_he_signal = df_he['signal_arb'].max()

    ax.plot(
        df_he['v_Aps'],
        df_he['signal_arb'] / max_he_signal,
        label="I$_2$He$_N$:I$^+$He (Droplet)",
        linestyle=":",
        color="#f03b20",  # Using a nice red
        linewidth=2
    )

    # Styling to match your publication/thesis format
    ax.set_xlim([0, 28])
    ax.set_ylim([0, 1.1])
    ax.set_xlabel('v / Å/ps', fontsize=12)
    ax.set_ylabel('Normalized Signal / arb. units', fontsize=12)
    ax.legend(frameon=False, fontsize=11)
    ax.set_title("VMI Reference Data Verification", fontsize=14)

    # Remove top and right spines for a cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Adjust this path if your script is in a different directory relative to data/
    # 1. Dynamically find the project root based on this file's location
    # compare.py is in scripts/post_processing_comparison/
    # .parents[2] goes up exactly three levels to the root 'i2_helium_md' folder
    PROJECT_ROOT = Path(__file__).resolve().parents[2]

    # 2. Build the absolute path to the reference data
    DATA_REF_DIR = PROJECT_ROOT / "data" / "reference"

    verify_vmi_data(data_dir=DATA_REF_DIR)