# HSA – Benzothiazole Docking and Molecular Dynamics

Computational pipeline for the study of benzothiazole derivatives as potential binders of Human Serum Albumin (HSA). This repository contains all analysis scripts used in the associated publication.

## Workflow Overview

```
01_prepare_receptor       → Clean PDB, fix missing atoms/H, convert to PDBQT
02_binding_site_prediction → P2Rank binding pocket prediction
03_residue_selection       → Extract pocket residues as individual PDB files
04_ligand_optimization     → ANI-2x geometry optimization (Neural Network Potential)
05_docking                 → GPU-accelerated docking with UniDock (Vina scoring)
06_build_topology          → AMBER19SB + GAFF2.11 + TIP3P topology (OpenMM)
07_interaction_fingerprint → Protein-ligand interaction fingerprint (ProLIF)
08_equilibration           → NVT + NPT equilibration (0.5 ns + 0.5 ns, OpenMM)
09_production_md           → NPT production MD — 21 × 5 ns = 105 ns (CUDA)
10_concatenate_trajectories→ PBC unwrap, complex centering, Cα alignment
11_gnina_rescoring         → CNN-based rescoring of MD frames with GNINA
12_pistack_analysis        → Pi-stacking centroid distances (TYR136, TYR159)
13_rmsd_ligand             → Ligand RMSD (no local fit, reference = frame 100)
14_rmsd_docking_pose       → Ligand RMSD vs initial docking pose
15_rmsd_protein_backbone   → Protein Cα RMSD over production trajectory
```

## Software & Force Fields

| Tool | Version | Purpose |
|------|---------|---------|
| OpenMM | 8.x | MD engine |
| AMBER ff19SB | — | Protein force field |
| GAFF2 (v2.11) | — | Ligand force field |
| TIP3P | — | Water model |
| UniDock | — | GPU docking |
| P2Rank | 2.4 | Binding site prediction |
| GNINA | — | CNN rescoring |
| ANI-2x | — | Ligand geometry optimization |
| MDAnalysis | — | Trajectory analysis |
| ProLIF | — | Interaction fingerprints |
| MDTraj | — | Pi-stacking analysis |

## Data Files

| File | Description |
|------|-------------|
| `data/gnina_affinity.csv` | Minimized affinity per MD frame (kcal/mol) |
| `data/gnina_cnn_score.csv` | GNINA CNN score per frame |
| `data/gnina_cnn_affinity.csv` | GNINA CNN affinity per frame |
| `data/pistack_centroid_distances.csv` | Centroid distances TYR136/TYR159 vs ligand (Å) |
| `data/rmsd_ligand.csv` | Ligand RMSD (ref = frame 100) |
| `data/rmsd_docking_pose.csv` | Ligand RMSD vs docking pose |
| `data/rmsd_protein_backbone.csv` | Protein Cα RMSD |

## Boltz-2 Predicted Poses

In addition to classical docking, Boltz-2 (a deep learning structure prediction model) was used to predict binding poses for three benzothiazole derivatives at site IB of HSA. The predicted structures are consistent with the docking results obtained at site IB.

| File | Description |
|------|-------------|
| `boltz2/ligand_01_boltz2.pdb` | Boltz-2 predicted pose — Ligand 01 |
| `boltz2/ligand_02_boltz2.pdb` | Boltz-2 predicted pose — Lig