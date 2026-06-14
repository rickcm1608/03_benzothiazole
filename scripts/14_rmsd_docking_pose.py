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
ref_structure = "system.pdb"   # rigid reference (docking pose)

protein_sel = "protein and backbone"
ligand_sel = "resname UNK and not name H*"

output_csv = "rmsd_docking_pose.csv"
output_plot = "rmsd_docking_pose.png"

print("=== Ligand RMSD vs docking pose (reference = system.pdb) ===\n")

# ===============================
# LOAD UNIVERSES
# ===============================
u = mda.Universe(topology, trajectory)      # trajectory
ref = mda.Universe(topology, ref_structure) # rigid reference

ligand = u.select_atoms(ligand_sel)
ligand_ref = ref.select_atoms(ligand_sel)

print(f"Ligand selection: {ligand_sel}")
print(f"Ligand atoms: {len(ligand)}\n")

# ===============================
# ALIGN PROTEIN TO REFERENCE
# ===============================
print("Aligning trajectory using protein backbone only...")

aligner = align.AlignTraj(
    u,
    ref,
    select=protein_sel,
    weights="mass",
    in_memory=True
)
aligner.run()

print("Alignment completed.\n")

# ===============================
# REFERENCE = LIGAND IN system.pdb (DOCKING POSE)
# ===============================
ref_positions = ligand_ref.positions.copy()

print("Ligand reference taken from system.pdb (docking pose)\n")

# ===============================
# COMPUTE LIGAND RMSD VS DOCKING POSE
# ===============================
print("Computing ligand RMSD (reference = docking pose)...")

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
df = pd.DataFrame({"Frame": frames, "RMSD_docking_pose": rmsd_values})
df.to_csv(output_csv, index=False)

print(f"✓ Data saved to: {output_csv}")

plt.figure(figsize=(7,4))
plt.plot(frames, rmsd_values, lw=1.2)
plt.xlabel("Frame")
plt.ylabel("Ligand RMSD (Å)")
plt.title("Ligand RMSD vs Docking Pose (reference = system.pdb)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(output_plot, dpi=300)

print(f"✓ Plot saved as: {output_plot}")
print("\n=== Analysis complete ===")
