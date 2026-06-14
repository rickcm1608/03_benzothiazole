#!/usr/bin/env python3
import os
from sys import stdout
import parmed as pmd
import pytraj as pt
from openmm.app import (PDBFile, AmberPrmtopFile, AmberInpcrdFile,
                        Simulation, DCDReporter, StateDataReporter,
                        PME, HBonds)
from openmm import (LangevinIntegrator, MonteCarloBarostat,
                    CustomExternalForce, unit, XmlSerializer)

# -----------------------------
# FILE CONFIGURATION
# -----------------------------
workDir = os.getcwd()
prmtop_file = os.path.join(workDir, "system_amber.prmtop")
inpcrd_file = os.path.join(workDir, "system_amber.inpcrd")

pdb_file_out = os.path.join(workDir, "prot_lig_equil.pdb")
dcd_file = os.path.join(workDir, "prot_lig_equil.dcd")
log_file = os.path.join(workDir, "prot_lig_equil.log")
rst_file = os.path.join(workDir, "prot_lig_equil.rst")
min_log = os.path.join(workDir, "prot_lig_equil_min.log")

# -----------------------------
# SIMULATION PARAMETERS
# -----------------------------
temperature = 300.0 * unit.kelvin
pressure = 1.0 * unit.bar
friction = 1.0/unit.picosecond
timestep = 2.0*unit.femtoseconds

nsteps_min = 25000
nvt_ns = 0.5
npt_ns = 0.5

restraint_fc = 50.0  # kJ/mol/nm^2
save_stride_ps = 10.0

# -----------------------------
# LOAD SYSTEM
# -----------------------------
print("Loading Amber system...")
prmtop = AmberPrmtopFile(prmtop_file)
inpcrd = AmberInpcrdFile(inpcrd_file)

system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=HBonds,
    rigidWater=True
)

# -----------------------------
# ATOM SELECTION FOR RESTRAINTS
# -----------------------------
print("Selecting atoms for restraints using negation logic...")

traj = pt.load(inpcrd_file, prmtop_file)
# Exclude HOH, NA, CL, and H* atoms
selection = traj.top.select("!:HOH,NA,CL & !@H*")
print(f">> Atoms to restrain: {len(selection)}")

# Apply positional restraints
if restraint_fc > 0:
    restraint = CustomExternalForce("k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addPerParticleParameter("k")
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for idx in selection:
        pos = inpcrd.positions[idx]
        x, y, z = pos[0].value_in_unit(unit.nanometer), pos[1].value_in_unit(unit.nanometer), pos[2].value_in_unit(unit.nanometer)
        restraint.addParticle(idx, [restraint_fc, x, y, z])
    system.addForce(restraint)

# -----------------------------
# NVT INTEGRATOR AND SIMULATION
# -----------------------------
integrator_nvt = LangevinIntegrator(temperature, friction, timestep)
simulation_nvt = Simulation(prmtop.topology, system, integrator_nvt)
simulation_nvt.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation_nvt.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# -----------------------------
# ENERGY MINIMIZATION
# -----------------------------
print("Running energy minimization...")
simulation_nvt.minimizeEnergy(maxIterations=nsteps_min)
state = simulation_nvt.context.getState(getEnergy=True)
Emin = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
with open(min_log, "w") as f:
    f.write(f"Potential Energy after minimization (kJ/mol): {Emin}\n")
print(">> Minimization complete.")
print(f">> Energy = {Emin:.3f} kJ/mol")

# -----------------------------
# NVT PHASE (0.5 ns)
# -----------------------------
simulation_nvt.context.setVelocitiesToTemperature(temperature)
steps_nvt = int((nvt_ns * 1000 * unit.picoseconds) / timestep)
stride = int((save_stride_ps * unit.picoseconds) / timestep)

simulation_nvt.reporters.append(DCDReporter(dcd_file, stride))
simulation_nvt.reporters.append(StateDataReporter(
    log_file, stride,
    step=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, density=True,
    separator=','
))
simulation_nvt.reporters.append(StateDataReporter(
    stdout, stride,
    step=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, density=True,
    separator=','
))

print(f"Running NVT for {nvt_ns} ns = {steps_nvt} steps...")
simulation_nvt.step(steps_nvt)

# -----------------------------
# NPT SIMULATION SETUP
# -----------------------------
system.addForce(MonteCarloBarostat(pressure, temperature))
integrator_npt = LangevinIntegrator(temperature, friction, timestep)
simulation_npt = Simulation(prmtop.topology, system, integrator_npt)

# Transfer positions and velocities from NVT
state_nvt = simulation_nvt.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
simulation_npt.context.setPositions(state_nvt.getPositions())
simulation_npt.context.setVelocities(state_nvt.getVelocities())
if inpcrd.boxVectors is not None:
    simulation_npt.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# NPT reporters
simulation_npt.reporters.append(DCDReporter(dcd_file, stride))
simulation_npt.reporters.append(StateDataReporter(
    log_file, stride,
    step=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, density=True,
    separator=','
))
simulation_npt.reporters.append(StateDataReporter(
    stdout, stride,
    step=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, density=True,
    separator=','
))

# -----------------------------
# NPT PHASE (0.5 ns)
# -----------------------------
steps_npt = int((npt_ns * 1000 * unit.picoseconds) / timestep)
print(f"Running NPT for {npt_ns} ns = {steps_npt} steps...")
simulation_npt.step(steps_npt)

# -----------------------------
# SAVE FINAL STATE
# -----------------------------
positions = simulation_npt.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation_npt.topology, positions, open(pdb_file_out, "w"))

state_final = simulation_npt.context.getState(getPositions=True, getVelocities=True)
with open(rst_file, "w") as f:
    f.write(XmlSerializer.serialize(state_final))

print(">> Done.")
print(f"Generated:\n - {pdb_file_out}\n - {dcd_file}\n - {log_file}\n - {rst_file}\n - {min_log}")
