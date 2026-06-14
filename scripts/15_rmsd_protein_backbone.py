import pytraj as pt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ============================
# Configuration
# ============================

traj_file = "prot_lig_prod_1-21_whole.dcd"
top_file = "system_amber.prmtop"

Output_name = "rmsd_protein_backbone"

workDir = os.getcwd()

# Time per frame (1050 frames = 105 ns -> 0.1 ns/frame)
time_per_frame = 0.1  # ns

# ============================
# Load data
# ============================

print(f"Loaded trajectory: {traj_file}")
print(f"Loaded topology: {top_file}")

traj = pt.load(traj_file, top=top_file)
n_frames = traj.n_frames

print(f"Frames: {n_frames}")

# ============================
# RMSD calculation (Cα atoms)
# ============================

rmsd = pt.rmsd(traj, ref=0, mask="@CA")

# Time array
time_array = np.arange(0, n_frames * time_per_frame, time_per_frame)

# ============================
# Plot
# ============================

plt.figure(figsize=(7,4))
plt.plot(time_array, rmsd, alpha=0.6, linewidth=1.0)
plt.xlabel("Time (ns)", fontsize=14, fontweight='bold')
plt.ylabel("RMSD (Å)", fontsize=14, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig(os.path.join(workDir, Output_name + ".png"),
            dpi=600, bbox_inches='tight')

# ============================
# Save raw data
# ============================

raw_data = pd.DataFrame({"time_ns": time_array, "rmsd_ca": rmsd})
raw_data.to_csv(os.path.join(workDir, Output_name + ".csv"),
                index=False)

print("Output written to:")
print(" -", os.path.join(workDir, Output_name + ".png"))
print(" -", os.path.join(workDir, Output_name + ".csv"))
