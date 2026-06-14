import os
import sys
import csv
import subprocess

# === Adjust local paths ===
workDir = os.getcwd()
receptor = os.path.join(workDir, "receptor.pdb")

# Path to P2Rank binary
prank_bin = "/home/jvaldiviezo/bin/p2rank_2.4/prank"

# Output folder
output_p2rank = os.path.join(workDir, "output_p2rank")

# Build P2Rank command
p2rank_cmd = f"{prank_bin} predict -f {receptor} -o {output_p2rank}"

# Save as executable script
script_path = os.path.join(workDir, "p2rank.sh")
with open(script_path, "w") as f:
    f.write(p2rank_cmd + "\n")

# Set permissions
subprocess.run(["chmod", "700", script_path])

# Run P2Rank via bash
subprocess.run(["bash", script_path])

# === Parse the CSV output ===
csv_file = os.path.join(output_p2rank, "receptor.pdb_predictions.csv")

with open(csv_file, "r") as f:
    csvreader = csv.reader(f)
    rows = list(csvreader)

# Columns
residue = []
score = []
cx = []
cy = []
cz = []

for row in rows:
    residue.append(row[9:10])
    score.append(row[2:3])
    cx.append(row[6:7])
    cy.append(row[7:8])
    cz.append(row[8:9])

# Print results for each predicted pocket
for i in range(1, len(residue)):
    file = str((residue[i])[0]).split()
    score_end = str((score[i])[0]).split()
    center_x_end = str((cx[i])[0]).split()
    center_y_end = str((cy[i])[0]).split()
    center_z_end = str((cz[i])[0]).split()

    print(f"Pocket {i}")
    print("Score =", score_end[0])

    final_res = []
    for r in file:
        final_res.append(int(r[2:]))  # Strip chain prefix (e.g. "A:" or "B:")

    print("Selected Residues =", final_res)
    print(f"Center x = {center_x_end[0]}  y = {center_y_end[0]}  z = {center_z_end[0]}\n")
