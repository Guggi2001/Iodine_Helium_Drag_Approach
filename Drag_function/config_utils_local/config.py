import os
from pathlib import Path

home = True
if home:
    BASE_PATH = r"C:\Users\paulg\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT"
else:
    BASE_PATH = r"T:\NextCloud_PaulGuggenbichler\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT"

PATH9A  = os.path.join(BASE_PATH, "9A", "9A_All_Data.csv")
PATH18A = os.path.join(BASE_PATH, "18A", "18A_All_Data.csv")
PATH18A_OLD = os.path.join(BASE_PATH, "18A")


# Repo paths follow this file, so they hold on every machine.
_REPO_ROOT = Path(__file__).resolve().parents[2]
TDDFT_9A_DIR = str(_REPO_ROOT / "i2_helium_md" / "data" / "reference" / "TDDFT" / "9A")

# Fit exports of drag_calculation.py: scratch output, gitignored, not reference data.
DRAG_OUTPUT_DIR = str(_REPO_ROOT / "Drag_function" / "output")
