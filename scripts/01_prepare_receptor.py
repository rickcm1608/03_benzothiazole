import os
import warnings
warnings.filterwarnings('ignore')

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import is_aa
from biopandas.pdb import PandasPdb
from pdbfixer import PDBFixer
from openmm.app.pdbfile import PDBFile

# Working directory
workDir = os.getcwd()
os.makedirs(workDir, exist_ok=True)

# Output files
temp = os.path.join(workDir, "temp.pdb")
receptor = os.path.join(workDir, "receptor.pdb")
receptor_pdbqt = os.path.join(workDir, "receptor.pdbqt")
ligand = os.path.join(workDir, "ligand.pdbqt")
name_residues_txt = os.path.join(workDir, "name_residues.txt")
name_residues_receptor_txt = os.path.join(workDir, "name_residues_receptor.txt")

# ---------- Prepare receptor ----------
receptor_file = os.path.join(workDir, "receptor.pdb")  # existing PDB file

# Remove previous output files if they exist
for f in [name_residues_txt, name_residues_receptor_txt]:
    if os.path.exists(f):
        os.remove(f)

# Read PDB with Biopandas
ppdb = PandasPdb().read_pdb(receptor_file)
ppdb.df['ATOM'] = ppdb.df['ATOM']
ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] != 'HOH']
ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] != 'OXT']
ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['element_symbol'] != 'H']
ppdb.to_pdb(path=receptor, records=['ATOM', 'HETATM'], gz=False, append_newline=True)
ppdb.to_pdb(path=temp, records=['ATOM', 'HETATM'], gz=False, append_newline=True)

# Use PDBFixer to add missing residues, atoms, and hydrogens
fixer = PDBFixer(filename=receptor)
fixer.removeHeterogens()
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open(receptor, 'w'))

# ---------- Functions to extract residues ----------
def aa(residue):
    res = residue.id[0]
    return res != "W"

class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue
    def accept_chain(self, chain):
        return chain.id == self.chain.id
    def accept_residue(self, residue):
        return residue == self.residue and aa(residue)

def extract_residues(input_pdb, output_txt):
    pdb = PDBParser().get_structure(input_pdb, input_pdb)
    io = PDBIO()
    io.set_structure(pdb)
    i = 1
    for model in pdb:
        for chain in model:
            for residue in chain:
                if not aa(residue):
                    continue
                print(f"saving {residue}", file=open(output_txt, "a"))
                i += 1

# Extract residues from both temp and fixed receptor
extract_residues(temp, name_residues_txt)
extract_residues(receptor, name_residues_receptor_txt)

# ---------- Convert receptor to PDBQT with OpenBabel ----------
os.system(f"obabel -i pdb {receptor} -o pdbqt -O {receptor_pdbqt} -xr --partialcharge")

# ---------- Read and display residues ----------
import pandas as pd

dataset = pd.read_csv(name_residues_txt, delimiter=" ", header=None)
df = dataset.iloc[:, [2]]
print("Residue - Number:")
for i, row in enumerate(df.values, 1):
    print(f"{row[0]} - {i}")
