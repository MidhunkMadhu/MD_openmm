# MD_openmm

## Overview

MD_openmm is a small OpenMM-based molecular dynamics workflow for a membrane protein system prepared using the CHARMM-GUI web server. OpenMM is a high-performance Python toolkit for molecular simulation. It can be used as a standalone simulation engine or as a flexible Python library for building custom molecular dynamics workflows.

This repository demonstrates how to continue a CHARMM-GUI-prepared membrane protein simulation in OpenMM. The CHARMM-GUI default equilibration protocol is followed up to step 6.6, which represents the final restrained equilibration stage. A new OpenMM script is then used for step 6.7, where a 20 ns unrestrained equilibration is performed. This step is equivalent to a production-level MD run and prepares the system for independent production-level molecular dynamics simulations.

The workflow is written for a membrane-containing system, but the production script is general enough to be adapted to other membrane protein systems prepared with CHARMM-GUI or converted into compatible OpenMM input formats.

## Installing OpenMM

A simple way to install OpenMM is through conda or mamba using the conda-forge channel.

```bash
conda create -n openmm -c conda-forge python=3.11 openmm parmed
conda activate openmm
```

Alternatively, using mamba:

```bash
mamba create -n openmm -c conda-forge python=3.11 openmm parmed
mamba activate openmm
```

ParmEd is included because the production script can write AMBER-style restart files during and after the simulation.

To check that OpenMM is available:

```bash
python -c "import openmm; print(openmm.version.version)"
```

To check available OpenMM platforms:

```bash
python -m openmm.testInstallation
```

For GPU runs, make sure the CUDA/OpenCL drivers are correctly installed on the machine or cluster where the simulation will be executed.

## System preparation

The membrane protein system in this repository was generated using the CHARMM-GUI web server. CHARMM-GUI was used to prepare the membrane protein system, solvate it, add ions, generate topology and parameter files, and create the default equilibration inputs.

The repository contains files generated from the CHARMM-GUI workflow, including topology, coordinate, parameter, equilibration, restart, and support files.

Example files in this repository include:

```text
b2r1.parm7
b2r1.rst7
step5_input.crd
step5_input.pdb
step5_input.psf
toppar.str
omm_barostat.py
omm_readinputs.py
omm_readparams.py
omm_restraints.py
omm_rewrap.py
omm_vfswitch.py
openmm_run.py
openmm_production.py
step6.1_equilibration.inp
step6.2_equilibration.inp
step6.3_equilibration.inp
step6.4_equilibration.inp
step6.5_equilibration.inp
step6.6_equilibration.inp
step6.7_equilibration.inp
```

The `step6.1` to `step6.6` files correspond to the default CHARMM-GUI equilibration protocol. These steps gradually equilibrate the system while applying restraints. The final restrained equilibration stage is step 6.6.

## Workflow

The intended workflow is:

```text
1. Build the membrane protein system using CHARMM-GUI
2. Run the default CHARMM-GUI equilibration protocol up to step 6.6
3. Use step6.6_equilibration.rst as the starting state
4. Run step 6.7 using openmm_production.py
5. Use the final step 6.7 restart/state files for independent production MD
```

In this repository, step 6.7 is a fresh unrestrained equilibration step added after the CHARMM-GUI restrained equilibration protocol. It runs for 20 ns without restraints and can be treated as a production-level equilibration before launching independent production simulations.

## Main production script

The main script is:

```text
openmm_production.py
```

This script is designed to run production-level or unrestrained equilibration MD using OpenMM. It reads all simulation settings from a single input file:

```text
step6.7_equilibration.inp
```

The script uses an AMBER-like starting-file logic. The input option:

```text
coordinate_file = ...
```

acts like the AMBER `-c` option. This means there is only one starting file. The same option can be used whether the starting file is an AMBER restart, an OpenMM XML state, an OpenMM checkpoint, or another supported coordinate format.

The script automatically detects the starting file type and loads it appropriately.

Supported starting formats include:

```text
AMBER restart or coordinate files
OpenMM XML state files
OpenMM checkpoint files
GROMACS GRO files
CHARMM coordinate files
```

This makes the script useful beyond this specific example system. It can be used for other membrane-containing systems as long as the required topology, coordinate/restart, and parameter files are available in a compatible format.

## Required input files

For the AMBER-style workflow demonstrated here, the main required files are:

```text
b2r1.parm7
step6.6_equilibration.rst
step6.7_equilibration.inp
openmm_production.py
omm_barostat.py
omm_restraints.py
omm_rewrap.py
omm_vfswitch.py
```

The topology file is:

```text
b2r1.parm7
```

The starting coordinate/state file for step 6.7 is:

```text
step6.6_equilibration.rst
```

Although the file extension is `.rst`, the script checks the file content and can detect whether it is an OpenMM XML state or an AMBER-style restart.

## Step 6.7 input file

The input file controls the topology file, starting file, output names, simulation length, timestep, temperature, pressure coupling, trajectory writing, checkpoint writing, and restart writing.

Example:

```text
force_field             = AMBER

topology_file           = b2r1.parm7
coordinate_file         = step6.6_equilibration.rst

output_xml              = step6.7_production.xml
output_dcd              = step6.7_production.dcd
output_report           = step6.7_production_report.txt
progress_file           = step6.7_production_progress.out

final_amber_restart     = step6.7_production.rst7
restart_snapshot_prefix = step6.7_production_snapshot.rst7

continuation            = yes
genvel                  = no
reset_step_and_time     = no
rewrap_coordinates      = yes
```

For true continuation from step 6.6:

```text
continuation = yes
genvel       = no
```

This keeps the velocities stored in the starting file.

For a fresh run from coordinates where new velocities should be generated:

```text
continuation = no
genvel       = yes
```

This generates velocities from the temperature specified in the input file.

By default, the script does not reset the step number or simulation time. If a reset is needed, it can be explicitly requested:

```text
reset_step_and_time = yes
```

## Running the simulation

Activate the OpenMM environment:

```bash
conda activate openmm
```

Run step 6.7:

```bash
python openmm_production.py step6.7_equilibration.inp
```

The default setup performs a 20 ns unrestrained equilibration.

## Output files

The script writes trajectory, progress, restart, checkpoint, and report files.

Typical outputs include:

```text
step6.7_production.dcd
step6.7_production_progress.out
step6.7_production.xml
step6.7_production.xml.chk
step6.7_production.rst7
step6.7_production_snapshot.rst7.<step>
step6.7_production_report.txt
```

The trajectory file is:

```text
step6.7_production.dcd
```

The progress file is:

```text
step6.7_production_progress.out
```

This file records quantities such as simulation step, time, potential energy, temperature, progress, remaining time, and simulation speed.

The final OpenMM XML state is:

```text
step6.7_production.xml
```

This can be used to continue a later OpenMM simulation while preserving the OpenMM state.

The OpenMM checkpoint is:

```text
step6.7_production.xml.chk
```

This can be useful for restarting on the same or compatible platform.

The final AMBER-style restart file is:

```text
step6.7_production.rst7
```

This can be used as a coordinate/restart file for later workflows that support AMBER restart formats.

Periodic AMBER restart snapshots are written using the prefix:

```text
step6.7_production_snapshot.rst7
```

These are useful for recovering intermediate simulation states or launching multiple downstream simulations.

The final report file is:

```text
step6.7_production_report.txt
```

This summarizes the run settings, input files, detected starting-file type, output files, timing, and performance.

## Using the output for further production MD

After step 6.7 is complete, the system can be used for independent production-level MD simulations.

Possible starting files for future runs include:

```text
step6.7_production.xml
step6.7_production.xml.chk
step6.7_production.rst7
```

For OpenMM-based continuation, the XML state or checkpoint is usually the most direct option. For workflows that need AMBER-style restart files, the `.rst7` file can be used.

A later production run can use the same script by changing the input file, for example:

```text
topology_file   = b2r1.parm7
coordinate_file = step6.7_production.xml

continuation    = yes
genvel          = no
```

or:

```text
topology_file   = b2r1.parm7
coordinate_file = step6.7_production.rst7

continuation    = yes
genvel          = no
```

This allows the same script to be reused for additional production runs, replicate simulations, or longer MD extensions.

