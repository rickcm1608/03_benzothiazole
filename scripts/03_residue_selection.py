import os
import pandas as pd
import re
from Bio.PDB import PDBParser, PDBIO, Select

workDir = os.getcwd()
temp = os.path.join(workDir, "temp.pdb")
receptor = os.path.join(workDir, "receptor.pdb")

# ======== Create subfolder for selected residues ========
selected_dir = os.path.join(workDir, "selected_residues")
os.makedirs(selected_dir, exist_ok=True)

# ======== Load P2Rank prediction ========
df = pd.read_csv("output_p2rank/receptor.pdb_predictions.csv")

# Strip whitespace from column names
df.columns = df.columns.str.strip()

# Columns from P2Rank official CSV
residue = df["residue_ids"].tolist()
score = df["score"].tolist()
center_x = df["center_x"].tolist()
center_y = df["center_y"].tolist()
center_z = df["center_z"].tolist()

# ======== Helper functions ========
def aa(residue):
    return residue.id[0] != "W"

class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue

    def accept_chain(self, chain):
        return chain.id == self.chain.id

    def accept_residue(self, residue):
        return residue == self.residue and aa(residue)

def extract_residue(target_number):
    pdb = PDBParser().get_structure("temp", temp)
    i = 1
    for model in pdb:
        for chain in model:
            for residue in chain:
                if not aa(residue):
                    continue
                if i == target_number:
                    outname = f"res_{target_number}.pdb"
                    outpath = os.path.join(selected_dir, outname)
                    io = PDBIO()
                    io.set_structure(pdb)
                    io.save(outpath, ResidueSelect(chain, residue))
                    return outpath
                i += 1
    return None

# ======== User selection ========
Selection = "Pocket"     # Change to "Residues" to select manually
number = "0"             # Pocket index or residue list: "147,150,155"

if Selection == "Pocket":
    n = int(number)

    print(f"\n=== Pocket {n} selected ===")
    print(f"Score: {score[n]}")
    print(f"Center: {center_x[n]}, {center_y[n]}, {center_z[n]}")

    res_list = residue[n].split()
    residues_num = [int(r[2:]) for r in res_list]

else:
    residues_num = [int(x) for x in number.split(",")]

print("\nSelected residues:", residues_num)

# ======== Extract individual PDB files ========
pdb_list = []
for r in residues_num:
    out = extract_residue(r)
    pdb_list.append(out)

# ======== Merge PDB files ========
merge_file = os.path.join(workDir, "selection_merge.pdb")
with open(merge_file, "w") as out:
    for file in pdb_list:
        with open(file, "r") as f:
            for line in f:
                if not line.startswith("END"):
                    out.write(line)

print("\n✅ Individual residue files saved to:")
print(selected_dir)

print("\n✅ Merged selection file generated:")
print(merge_file)
