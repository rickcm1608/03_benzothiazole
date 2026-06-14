# -*- coding: utf-8 -*-
# =============================
# IMPORTS
# =============================
import os
import torch
import torchani
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Geometry import Point3D
from openbabel import pybel
from ase import io
from ase.optimize import BFGS
from torchani.units import HARTREE_TO_KCALMOL

# Path configuration
workDir = os.getcwd()
input_xyz = os.path.join(workDir, "ligand.xyz")

# Validate input file
if not os.path.exists(input_xyz):
    raise FileNotFoundError(f"File not found: {input_xyz}")

# =============================
# PART 1 - LOAD XYZ AND PERCEIVE BONDS
# =============================
print(f"\n[1/3] Loading {input_xyz} and perceiving connectivity...")

# OpenBabel infers bonds from atomic distances
mol_pybel = next(pybel.readfile("xyz", input_xyz))
mol_block = mol_pybel.write("mol")

# Convert to RDKit preserving hydrogens
hmol = Chem.MolFromMolBlock(mol_block, removeHs=False)
if hmol is None:
    raise ValueError("Error processing structure. Check the XYZ format.")

# Save initial reference structure
AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand_inicial.mol"))
print(" -> Connectivity perceived and saved to ligand_inicial.mol")

# Quick 2D visualization
smiles_2d = Chem.MolToSmiles(Chem.RemoveHs(hmol))
img_name = os.path.join(workDir, "ligand_2d_preview.png")
Draw.MolToFile(Chem.MolFromSmiles(smiles_2d), img_name)
print(f" -> Detected SMILES: {smiles_2d}")

# =============================
# PART 2 - GEOMETRY OPTIMIZATION (ANI-2x)
# =============================
print("\n[2/3] Starting optimization with ANI-2x (Neural Network Potential)...")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
calculator = torchani.models.ANI2x().ase()
model = torchani.models.ANI2x(periodic_table_index=True).to(device)

# Load in ASE for optimization
atoms = io.read(input_xyz)
atoms.set_calculator(calculator)

# Run BFGS optimizer
opt = BFGS(atoms)
opt.run(fmax=0.001)  # fmax is the force convergence threshold

# Save optimized coordinates to XYZ
io.write(os.path.join(workDir, "ligand_min.xyz"), atoms)

# =============================
# PART 3 - UPDATE COORDINATES AND EXPORT
# =============================
print("\n[3/3] Syncing coordinates and exporting formats...")

# Transfer new ASE coordinates to the RDKit object
new_pos = atoms.get_positions()
conf = hmol.GetConformer()
for i in range(hmol.GetNumAtoms()):
    x, y, z = new_pos[i]
    conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

# Export final files
AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand_min.mol"))
AllChem.MolToPDBFile(hmol, os.path.join(workDir, "ligand_min.pdb"))

# Generate PDBQT for docking
try:
    os.system("obabel -i mol ligand_min.mol -o pdbqt -O ligand_min.pdbqt -xh --partialcharge")
    print(" -> ligand_min.pdbqt generated successfully.")
except:
    print(" -> Error: Could not generate PDBQT. Make sure 'obabel' is in your PATH.")

# Compute final ANI-2x energy
def get_ani_energy(mol):
    pos = mol.GetConformer().GetPositions().tolist()
    atomnums = [a.GetAtomicNum() for a in mol.GetAtoms()]
    species = torch.tensor([atomnums], device=device)
    coords = torch.tensor([pos], requires_grad=True, device=device)
    return model((species, coords)).energies

energy_hartree = get_ani_energy(hmol).item()
print(f"\nFINAL RESULTS:")
print(f" - Energy: {energy_hartree:.6f} Hartree")
print(f" - Energy: {energy_hartree * HARTREE_TO_KCALMOL:.2f} kcal/mol")
print(f"\nDone. Files saved to: {workDir}")
