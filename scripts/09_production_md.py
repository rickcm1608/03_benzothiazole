#!/usr/bin/env python3
"""
NPT Production MD with OpenMM and GAFF 2.11
- Loads the equilibrated PDB and state file (.rst XML).
- Preserves periodic box.
- Saves DCD, PDB, state (.rst) and checkpoint at each stride.
- Recommended usage: launch with 'nohup python 09_production_md.py &'
"""

import os
from sys import stdout
import openmm
from openmm import unit, XmlSerializer, MonteCarloBarostat, LangevinIntegrator
from openmm.app import PDBFile, DCDReporter, StateDataReporter, CheckpointReporter, Simulation, PME, HBonds, ForceField
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator

# -----------------------------
# File configuration
# -----------------------------
workDir = os.getcwd()
Equilibrated_PDB = 'prot_lig_equil.pdb'
Equilibrated_RST = 'prot_lig_equil.rst'  # OpenMM XML state file
Ligand_file = 'ligand_prepared.sdf'
Jobname = 'prot_lig_prod'

# -----------------------------
# Simulation parameters
# -----------------------------
Stride_Time_ns = 5       # time per stride [ns]
Number_of_strides = 21     # number of strides
Integration_timestep_fs = 2  # fs
Temperature_K = 300       # K
Pressure_bar = 1          # bar
Frames_per_stride = 100   # frames saved per stride
Checkpoint_ps = 100       # checkpoint frequency in ps

# -----------------------------
# Unit conversion for OpenMM
# -----------------------------
stride_time_ps = Stride_Time_ns*1000
stride_time = stride_time_ps*unit.picoseconds
nstride = Number_of_strides
dt = Integration_timestep_fs*unit.femtoseconds
temperature = Temperature_K*unit.kelvin
pressure = Pressure_bar*unit.bar

# Steps per stride
nsteps = int(stride_time.value_in_unit(unit.picoseconds)/dt.value_in_unit(unit.picoseconds))

# Save exactly 100 frames per stride
nsavcrd = nsteps // Frames_per_stride
nchk = int(Checkpoint_ps*unit.picoseconds/dt)  # checkpoint interval in steps
total_steps = nsteps * nstride

# -----------------------------
# Load equilibrated system
# -----------------------------
pdbfile = os.path.join(workDir, Equilibrated_PDB)
rstfile = os.path.join(workDir, Equilibrated_RST)

pdb = PDBFile(pdbfile)
topology = pdb.topology
positions = pdb.positions

# -----------------------------
# Load ligand GAFF 2.11
# -----------------------------
ligand = Molecule.from_file(os.path.join(workDir, Ligand_file))
ligand_generator = GAFFTemplateGenerator(molecules=[ligand])

# -----------------------------
# Create force field
# -----------------------------
ff_protein = "amber19/protein.ff19SB.xml"
ff_water   = "amber19/tip3p.xml"
omm_forcefield = ForceField(ff_protein, ff_water)
omm_forcefield.registerTemplateGenerator(ligand_generator.generator)

# -----------------------------
# Create system
# -----------------------------
system = omm_forcefield.createSystem(
    topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)
system.addForce(MonteCarloBarostat(pressure, temperature))

# -----------------------------
# Integrator
# -----------------------------
friction = 1.0/unit.picoseconds
integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(1e-6)

platform = openmm.Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue("CudaDeviceIndex", "0")  # use GPU 0
platform.setPropertyDefaultValue("CudaPrecision", "mixed")

simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)

if pdb.topology.getPeriodicBoxVectors() is not None:
    simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# -----------------------------
# Production stride loop
# -----------------------------
for stride in range(1, nstride+1):
    print(f"\n>>> Stride #{stride} <<<")

    dcd_file = f"{Jobname}_{stride}.dcd"
    log_file = f"{Jobname}_{stride}.log"
    rst_file = f"{Jobname}_{stride}.rst"
    chk_file = f"{Jobname}_{stride}.chk"
    prev_rst = f"{Jobname}_{stride-1}.rst"

    # Load previous state
    if stride == 1:
        print(f"> Loading equilibration state: {rstfile}")
        with open(rstfile, 'r') as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))
    else:
        print(f"> Loading previous stride state: {prev_rst}")
        with open(prev_rst, 'r') as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))

    # Reporters
    simulation.reporters.append(DCDReporter(dcd_file, nsavcrd))
    simulation.reporters.append(CheckpointReporter(chk_file, nchk))
    simulation.reporters.append(StateDataReporter(
        stdout, nsavcrd,
        step=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=total_steps,
        separator='\t\t',
        kineticEnergy=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True
    ))
    simulation.reporters.append(StateDataReporter(
        log_file, nsavcrd,
        step=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=total_steps,
        separator='\t\t',
        kineticEnergy=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True
    ))

    print(f"> Simulating {nsteps} steps (Stride #{stride})...")
    simulation.step(nsteps)

    # Save final state
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(rst_file, 'w') as f:
        f.write(XmlSerializer.serialize(state))

    # Save final coordinates
    positions = state.getPositions()
    pdb_out = f"{Jobname}_{stride}.pdb"
    PDBFile.writeFile(simulation.topology, positions, open(pdb_out, 'w'))

    simulation.reporters.clear()  # Clear reporters for the next stride

print("\n>>> Simulation completed successfully <<<")
