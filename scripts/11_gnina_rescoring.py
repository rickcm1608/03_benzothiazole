import MDAnalysis as mda
from MDAnalysis.analysis import align
import subprocess, os, sys
import pandas as pd
from rdkit import Chem

# ====== CONFIGURATION ======
workDir = os.getcwd()
pdb_ref = "prot_lig_prod_nw.pdb"
traj = "prot_lig_prod_1-21_nw.dcd"
minimization = "Yes"
Skip = 1

# ====== SETUP ======
pdb_dir = os.path.join(workDir, "PDBs")
os.makedirs(pdb_dir, exist_ok=True)
output = workDir

u1 = mda.Universe(pdb_ref, traj)
u2 = mda.Universe(pdb_ref)
align.AlignTraj(u1, u2, select="name CA", in_memory=True).run()

protein_ligand = u1.select_atoms("not resname HOH NA CL K")
with mda.Writer(os.path.join(output, "protein_ligand_GNINA.dcd"), protein_ligand.n_atoms) as W:
    for ts in u1.trajectory[::Skip]:
        W.write(protein_ligand)

minimize = "--minimize" if minimization == "Yes" else "--score_only"

# ====== GNINA RESCORING ======
affinity, CNNscore, CNNaffinity = [], [], []
prot = u1.select_atoms("protein")
lig = u1.select_atoms("resname UNK")

for i, ts in enumerate(u1.trajectory[::Skip]):
    prot.write(os.path.join(pdb_dir, "protein.pdb"))
    lig.write(os.path.join(pdb_dir, "ligand.pdb"))

    cmd = f"/home/jvaldiviezo/bin/gnina -r {pdb_dir}/protein.pdb -l {pdb_dir}/ligand.pdb {minimize} -o gnina_scored_pose.sdf"
    with open("gnina_score.sh", "w") as f:
        f.write(cmd + "\n")

    subprocess.run("chmod 700 gnina_score.sh && ./gnina_score.sh", shell=True)

    res = Chem.SDMolSupplier("gnina_scored_pose.sdf", False)[0]
    affinity.append(float(res.GetProp("minimizedAffinity")))
    CNNscore.append(float(res.GetProp("CNNscore")))
    CNNaffinity.append(float(res.GetProp("CNNaffinity")))

    print(f"Frame {i} done.")

# ====== EXPORT ======
pd.DataFrame(affinity).to_csv(os.path.join(output, "gnina_affinity.csv"))
pd.DataFrame(CNNscore).to_csv(os.path.join(output, "gnina_cnn_score.csv"))
pd.DataFrame(CNNaffinity).to_csv(os.path.join(output, "gnina_cnn_affinity.csv"))
