import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# CONFIGURATION
# ===============================
topology = "system_amber.prmtop"
trajectory = "prot_lig_prod_1-21_whole.dcd"
ref_structure = "system.pdb"   # reference for protein alignment

protein_sel = "backbone"
ligand_sel = "resname UNK and not name H*"

output_csv = "rmsd_ligand.csv"
output_plot = "rmsd_ligand.png"

print("=== Ligand RMSD without local fit (reference = frame 100) ===\n")

# ===============================
# LOAD UNIVERSES
# ===============================
u = mda.Universe(topology, trajectory)          # trajectory
ref = mda.Universe(topology, ref_structure)     # external reference for protein

ligand = u.select_atoms(ligand_sel)

print(f"Ligand selection: {ligand_sel}")
print(f"Ligand atoms: {len(ligand)}\n")

# ===============================
# ALIGN PROTEIN TO REFERENCE
# ===============================
print("Aligning protein backbone to reference structure...")

aligner = align.AlignTraj(
    u,
    ref,
    select=protein_sel,
    in_memory=True
)
aligner.run()

print("Alignment completed.\n")

# ===============================
# LIGAND REFERENCE = FRAME 100 AFTER ALIGNMENT
# ===============================
u.trajectory[100]
ref_positions = ligand.positions.copy()

print("Ligand reference: positions at frame 100\n")

# ===============================
# COMPUTE LIGAND RMSD WITHOUT LOCAL FIT
# ===============================
print("Computing ligand RMSD (no fit, reference = frame 100)...")

rmsd_values = []

for ts in u.trajectory:
    diff = ligand.positions - ref_positions
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    rmsd_values.append(rmsd)

frames = np.arange(len(rmsd_values))
print(f"Frames analyzed: {len(frames)}\n")

# ===============================
# SAVE OUTPUT
# ===============================
df = pd.DataFrame({"Frame": frames, "RMSD": rmsd_values})
df.to_csv(output_csv, index=False)

print(f"✅ Data saved to: {output_csv}")

plt.figure(figsize=(7,4))
plt.plot(frames, rmsd_values, lw=1.2)
plt.xlabel("Frame")
plt.ylabel("Ligand RMSD (Å)")
plt.title("Ligand RMSD (no fit, reference = frame 100)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(output_plot, dpi=300)

print(f"✅ Plot saved as: {output_plot}")
print("\n=== Analysis complete ===")
