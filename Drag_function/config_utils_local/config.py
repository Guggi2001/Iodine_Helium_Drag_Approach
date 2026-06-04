import os

home = False
if home:
    BASE_PATH = r"C:\Users\paulg\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT"
else:
    BASE_PATH = r"T:\NextCloud_PaulGuggenbichler\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT"

PATH9A  = os.path.join(BASE_PATH, "9A", "9A_All_Data.csv")
PATH18A = os.path.join(BASE_PATH, "18A", "18A_All_Data.csv")
PATH18A_OLD = os.path.join(BASE_PATH, "18A")


