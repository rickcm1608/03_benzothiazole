#!/usr/bin/env python3
import os
import numpy as np
import parmed as pmd
from rdkit import Chem

# OpenMM imports
try:
    import openmm
    from openmm import app, unit, Vec3
except ImportError:
    from simtk import openmm, app, unit
    from simtk.openmm import Vec3

from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule
from openff.units.openmm import to_openmm

# -----------------------------
# FILE CONFIGURATION
# -----------------------------
workDir = os.getcwd()
ligand_out_sdf = os.path.join(workDir, "ligand_out.sdf")
receptor = os.path.join(workDir, "receptor.pdb")

# -----------------------------
# GENERAL PARAMETERS
# -----------------------------
Force_field = "AMBER19SB"
Water_type = "TIP3P"
Padding_distance = 4.0       # angstroms
Ions = "NaCl"
Concentration = 0.15         # M
pH = 7.4

# -----------------------------
# READ LIGAND
# -----------------------------
results = Chem.SDMolSupplier(ligand_out_sdf, removeHs=False)
ligand_mol = results[0]
Chem.MolToMolFile(ligand_mol, os.path.join(workDir, "ligand_prepared.sdf"))

ligand = Molecule.from_file(os.path.join(workDir, "ligand_prepared.sdf"))
ligand_positions = ligand.conformers[0]
ligand_topology = ligand.to_topology()

# -----------------------------
# LOAD FORCE FIELDS
# -----------------------------
ff_protein = "amber19/protein.ff19SB.xml"
ff_water = "amber19/tip3p.xml"
omm_forcefield = app.ForceField(ff_protein, ff_water)

# Ligand with GAFF 2.11
ligand_generator = GAFFTemplateGenerator(molecules=[ligand])
omm_forcefield.registerTemplateGenerator(ligand_generator.generator)

print("Force field XML files loaded.")

# -----------------------------
# LOAD RECEPTOR PDB
# -----------------------------
pdb = app.PDBFile(receptor)
modeller = app.Modeller(pdb.topology, pdb.positions)

# -----------------------------
# ADD LIGAND
# -----------------------------
modeller.add(ligand_topology.to_openmm(), to_openmm(ligand_positions))

# -----------------------------
# ADD HYDROGENS AT GIVEN PH
# -----------------------------
modeller.addHydrogens(omm_forcefield, pH=pH)

# -----------------------------
# SOLVATION AND IONS
# -----------------------------
positive_ion = 'Na+' if Ions == 'NaCl' else 'K+'

modeller.addSolvent(
    omm_forcefield,
    model='tip3p',
    padding=Padding_distance*unit.angstrom,
    ionicStrength=float(Concentration)*unit.molar,
    positiveIon=positive_ion,
    negativeIon='Cl-'
)

# -----------------------------
# GET TOPOLOGY AND POSITIONS
# -----------------------------
topology = modeller.getTopology()
positions = modeller.getPositions()

# -----------------------------
# DEFINE PERIODIC BOX (Vec3)
# -----------------------------
pos_nm = np.array(positions.value_in_unit(unit.nanometer))
min_pos = pos_nm.min(axis=0)
max_pos = pos_nm.max(axis=0)
vec = max_pos - min_pos + Padding_distance/10.0  # extra padding in nm

box_vectors = (Vec3(vec[0], 0, 0),
               Vec3(0, vec[1], 0),
               Vec3(0, 0, vec[2]))

modeller.topology.setPeriodicBoxVectors(box_vectors)
print(f"Periodic box (nm): {vec}")

# -----------------------------
# EXPORT PDB
# -----------------------------
system_pdb = os.path.join(workDir, "system.pdb")
app.PDBFile.writeFile(topology, positions, open(system_pdb, 'w'))
print(f"PDB file generated: {system_pdb}")

# -----------------------------
# CREATE SYSTEM FOR PARMED
# -----------------------------
system = omm_forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0*unit.nanometer,
                                     constraints=app.HBonds,
                                     flexibleConstraints=True)
pmd_struct = pmd.openmm.load_topology(topology, system, positions)

# -----------------------------
# EXPORT AMBER FILES
# -----------------------------
amber_prmtop = os.path.join(workDir, "system_amber.prmtop")
amber_inpcrd = os.path.join(workDir, "system_amber.inpcrd")
pmd_struct.save(amber_prmtop, overwrite=True)
pmd_struct.save(amber_inpcrd, overwrite=True)

print("AMBER files generated: ", amber_prmtop, amber_inpcrd)
print("Periodic topology generated successfully!")
