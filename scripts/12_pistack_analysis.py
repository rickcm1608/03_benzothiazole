import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. Load trajectory and topology
traj = md.load('prot_lig_prod_1-21_whole.dcd', top='system_amber.prmtop')

# 2. ATOM SELECTION BY ABSOLUTE INDEX (based on your PDB)
# Note: MDTraj uses 0-based indexing, so subtract 1 from PDB numbers
# TYR 136: CG(2128), CD1(2129), CD2(2131), CE1(2133), CE2(2135), CZ(2137)
idx_136 = [2127, 2128, 2130, 2132, 2134, 2136]

# TYR 159: CG(2547), CD1(2548), CD2(2550), CE1(2552), CE2(2554), CZ(2556)
idx_159 = [2546, 2547, 2549, 2551, 2553, 2555]

# LIGAND UNK: aromatic ring atoms
nombres_lig = ['C3x', 'C4x', 'C5x', 'C6x', 'C7x', 'C8x']
idx_lig = [a.index for a in traj.topology.atoms if a.name in nombres_lig]

# --- VERIFICATION ---
print(f"TYR136 atoms selected: {[traj.topology.atom(i).name for i in idx_136]}")
print(f"TYR159 atoms selected: {[traj.topology.atom(i).name for i in idx_159]}")
print(f"Ligand atoms selected: {[traj.topology.atom(i).name for i in idx_lig]}")

# 3. Centroid calculation
def get_manual_centroid(traj, indices):
    return np.mean(traj.xyz[:, indices, :], axis=1)

c_136 = get_manual_centroid(traj, idx_136)
c_159 = get_manual_centroid(traj, idx_159)
c_lig = get_manual_centroid(traj, idx_lig)

# 4. Distance calculation (nm -> Å)
def calc_dist(p1, p2):
    return np.sqrt(np.sum((p1 - p2)**2, axis=1)) * 10

dist_136 = calc_dist(c_136, c_lig)
dist_159 = calc_dist(c_159, c_lig)

# 5. Save CSV
df = pd.DataFrame({
    'Time_ns': np.linspace(0, 105, len(traj)),
    'Dist_TYR136_A': dist_136,
    'Dist_TYR159_A': dist_159
})
df.to_csv('pistack_centroid_distances.csv', index=False)

# 6. Plot
plt.figure(figsize=(10, 5))
plt.plot(df['Time_ns'], dist_136, label='Dist. TYR136', color='#2c3e50')
plt.plot(df['Time_ns'], dist_159, label='Dist. TYR159', color='#e74c3c')
plt.axhline(y=4.0, color='green', linestyle='--', label='Target ~4.0 Å')
plt.ylabel('Distance (Å)')
plt.xlabel('Time (ns)')
plt.title('Pi-Pi Stacking: Centroid Distances (TYR136/TYR159 vs Ligand)')
plt.legend()
plt.savefig('pistack_centroid_distances.png')
plt.show()

print(f"\n--- Final Statistics ---")
print(f"Mean TYR136: {np.mean(dist_136):.2f} ± {np.std(dist_136):.2f} Å")
print(f"Mean TYR159: {np.mean(dist_159):.2f} ± {np.std(dist_159):.2f} Å")
