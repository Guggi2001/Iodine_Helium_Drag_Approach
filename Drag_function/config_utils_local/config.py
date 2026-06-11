import os

home = False
if home:
    BASE_PATH = r"C:\Users\paulg\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT"
else:
    BASE_PATH = r"T:\NextCloud_PaulGuggenbichler\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT"

PATH9A  = os.path.join(BASE_PATH, "9A", "9A_All_Data.csv")
PATH18A = os.path.join(BASE_PATH, "18A", "18A_All_Data.csv")
PATH18A_OLD = os.path.join(BASE_PATH, "18A")


if home:
    BASE_PATH_i2_helium_md = r"C:\Users\paulg\Dokumente\GitHub\Iodine_Helium_Drag_Approach\i2_helium_md"
else:
    BASE_PATH_i2_helium_md = r"T:\github synchronized\Iodine_Helium_Drag_Approach\i2_helium_md"