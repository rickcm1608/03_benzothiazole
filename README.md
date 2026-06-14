# HSA – Benzothiazole Docking and Molecular Dynamics

Computational pipeline for the study of benzothiazole derivatives as potential binders of Human Serum Albumin (HSA). This repository contains all analysis scripts used in the associated publication.

## Workflow Overview

```
01_prepare_receptor       → Clean PDB, fix missing atoms/H, convert to PDBQT
02_binding_site_prediction → P2Rank binding pocket prediction
03_residue_selection       → Extract pocket residues as individual PDB files
04_ligand_optimization     → ANI-2x geometry optimization (Neural Network Potential)
05_docking                 → GPU-accelerated docking with UniDock (Vina scoring)
06_build_topology          → AMBER19SB + GAFF2 + TIP3P topology (OpenMM)
07_interaction_fingerprint → Protein-ligand interaction fingerprint (ProLIF)
08_equilibration           → NVT + NPT equilibration (0.5 ns + 0.5 ns, OpenMM)
09_production_md           → NPT production MD — 21 x 5 ns = 105 ns (CUDA)
10_concatenate_trajectories→ PBC unwrap, complex centering, Ca alignment
11_gnina_rescoring         → CNN-based rescoring of MD frames with GNINA
12_pistack_analysis        → Pi-stacking centroid distances (TYR136, TYR159)
13_rmsd_ligand             → Ligand RMSD (no local fit, reference = frame 100)
14_rmsd_docking_pose       → Ligand RMSD vs initial docking pose
15_rmsd_protein_backbone   → Protein Ca RMSD over production trajectory
```

## Software & Force Fields

| Tool | Purpose |
|------|---------|
| OpenMM | MD engine |
| AMBER ff19SB | Protein force field |
| GAFF2 | Ligand force field |
| TIP3P | Water model |
| UniDock | GPU docking |
| P2Rank | Binding site prediction |
| GNINA | CNN rescoring |
| ANI-2x | Ligand geometry optimization |
| MDAnalysis | Trajectory analysis |
| ProLIF | Interaction fingerprints |
| MDTraj | Pi-stacking analysis |

## Data Files

| File | Description |
|------|-------------|
| `data/gnina_affinity.csv` | Minimized affinity per MD frame (kcal/mol) |
| `data/gnina_cnn_score.csv` | GNINA CNN score per frame |
| `data/gnina_cnn_affinity.csv` | GNINA CNN affinity per frame |
| `data/pistack_centroid_distances.csv` | Centroid distances TYR136/TYR159 vs ligand (A) |
| `data/rmsd_ligand.csv` | Ligand RMSD (ref = frame 100) |
| `data/rmsd_docking_pose.csv` | Ligand RMSD vs docking pose |
| `data/rmsd_protein_backbone.csv` | Protein Ca RMSD |

## Boltz-2 Predicted Poses

In addition to classical docking, Boltz-2 (protein-ligand co-folding, via Rowan Scientific) was used to independently predict the binding pose and estimate binding affinity for the benzothiazole derivative (LB1) at site IB of HSA. Three replicate runs were performed to assess reproducibility.

| Replicate | pTM | ipTM | Affinity Prob. | Predicted pIC50 | Predicted IC50 (M) |
|-----------|-----|------|----------------|-----------------|---------------------|
| 1 | 0.927 | 0.759 | 0.417 | 6.50 | 3.17e-7 |
| 2 | 0.922 | 0.752 | 0.347 | 6.33 | 4.67e-7 |
| 3 | 0.927 | 0.725 | 0.408 | 6.52 | 3.01e-7 |

All three replicates converge at site IB with consistent pIC50 (~6.4) and pTM > 0.92, confirming the binding prediction is reproducible and consistent with docking and MD results.

| File | Description |
|------|-------------|
| `boltz2/ligand_01_boltz2.pdb` | Boltz-2 predicted pose — Replicate 1 |
| `boltz2/ligand_02_boltz2.pdb` | Boltz-2 predicted pose — Replicate 2 |
| `boltz2/ligand_03_boltz2.pdb` | Boltz-2 predicted pose — Replicate 3 |
| `boltz2/boltz2_results.csv` | Full Boltz-2 scoring results |

## Key Results

**GNINA rescoring** (minimized affinity, kcal/mol, over 21 x 5 ns MD frames):
- Mean: -11.67 kcal/mol | Best frame: -12.37 kcal/mol
- Consistently strong binding predicted across all frames, indicating stable occupancy at the HSA binding site.

**Pi-stacking with TYR136 / TYR159** (centroid-centroid distance, A):
- Mean: 4.27 A | Min: 3.26 A | Max: 6.19 A
- The benzothiazole ring maintains persistent pi-pi stacking contacts with both tyrosine residues throughout the trajectory.

**Ligand RMSD** (no local fit, ref = frame 100):
- Mean: 3.72 A | Max: 6.44 A

**Protein backbone RMSD** (Ca):
- Mean: 3.19 A | Max: 4.25 A

**Boltz-2 predicted binding affinity** (3 replicates, ligand LB1):
- Mean pIC50: 6.45 | Mean ipTM: 0.745 | Mean affinity probability: 0.39
- Boltz-2 confirms binding at site IB, consistent with docking and MD results.

## Structure Inputs

- `receptor.pdb` — HSA crystal structure (cleaned, H-added at pH 7.4)
- `ligand.xyz` — Benzothiazole derivative initial geometry

## Notes

- Large trajectory files (`.dcd`, `.rst`, `.chk`, `.prmtop`, `.inpcrd`) are excluded via `.gitignore`. Contact the authors for trajectory data.
- Binary paths (P2Rank, UniDock, GNINA) are hardcoded for the HPC environment used. Update these before running locally.

## Citation

> *To be added upon publication.*
