#!/usr/bin/env python3

from __future__ import print_function

import argparse
import datetime
import os
import sys
import time
import xml.etree.ElementTree as ET
from types import SimpleNamespace

from openmm import *
from openmm.app import *
from openmm.unit import *

from omm_vfswitch import *
from omm_barostat import *
from omm_restraints import *
from omm_rewrap import *


# ============================================================
# Input parsing
# ============================================================

def strip_comment(line):
    return line.split("#", 1)[0].strip()


def read_key_value_file(filename):
    data = {}

    with open(filename, "r") as f:
        for raw in f:
            line = strip_comment(raw)

            if not line:
                continue

            if "=" not in line:
                continue

            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            data[key] = value

    return data


def normalize_none(value):
    if value is None:
        return None

    v = str(value).strip()

    if v.lower() in ["none", "no", "null", "nil", "-", ""]:
        return None

    return v


def get_str(data, key, default=None, required=False):
    if key in data:
        return data[key]

    if required:
        sys.exit("Missing required key in input file: %s" % key)

    return default


def get_int(data, key, default=None, required=False):
    value = get_str(data, key, default, required)

    if value is None:
        return None

    return int(value)


def get_float(data, key, default=None, required=False):
    value = get_str(data, key, default, required)

    if value is None:
        return None

    return float(value)


def get_bool(data, key, default=False):
    if key not in data:
        return bool(default)

    value = str(data[key]).strip().lower()

    if value in ["yes", "true", "on", "1"]:
        return True

    if value in ["no", "false", "off", "0"]:
        return False

    sys.exit("Invalid yes/no value for %s: %s" % (key, data[key]))


def get_genvel(data):
    """
    Velocity policy.

    Default:
        genvel = yes

    genvel = yes:
        Generate new velocities from temp after coordinate_file is loaded.

    genvel = no:
        Do not generate new velocities.
        Use velocities from coordinate_file.
        Stop if coordinate_file has no velocities.

    Legacy compatibility:
        If old key generate_velocities_if_missing exists, use it as genvel.
    """

    if "genvel" in data:
        return get_bool(data, "genvel", True)

    if "generate_velocities_if_missing" in data:
        return get_bool(data, "generate_velocities_if_missing", True)

    return True


# ============================================================
# OpenMM option mapping
# ============================================================

def get_nonbonded_method(name):
    value = str(name).strip()

    mapping = {
        "NoCutoff": NoCutoff,
        "CutoffNonPeriodic": CutoffNonPeriodic,
        "CutoffPeriodic": CutoffPeriodic,
        "PME": PME,
        "Ewald": Ewald,
        "LJPME": LJPME,
    }

    if value not in mapping:
        sys.exit("Unknown coulomb/nonbonded method: %s" % value)

    return mapping[value]


def get_constraints(name):
    value = str(name).strip()

    mapping = {
        "None": None,
        "none": None,
        "HBonds": HBonds,
        "AllBonds": AllBonds,
        "HAngles": HAngles,
    }

    if value not in mapping:
        sys.exit("Unknown constraints option: %s" % value)

    return mapping[value]


def get_membrane_xy_mode(name):
    value = str(name).strip()

    mapping = {
        "XYIsotropic": MonteCarloMembraneBarostat.XYIsotropic,
        "XYAnisotropic": MonteCarloMembraneBarostat.XYAnisotropic,
    }

    if value not in mapping:
        sys.exit(
            "Unknown p_XYMode: %s\n"
            "Allowed values: %s"
            % (value, ", ".join(sorted(mapping.keys())))
        )

    return mapping[value]


def get_membrane_z_mode(name):
    value = str(name).strip()

    mapping = {
        "ZFree": MonteCarloMembraneBarostat.ZFree,
        "ZFixed": MonteCarloMembraneBarostat.ZFixed,
        "ConstantVolume": MonteCarloMembraneBarostat.ConstantVolume,
    }

    if value not in mapping:
        sys.exit(
            "Unknown p_ZMode: %s\n"
            "Allowed values: %s"
            % (value, ", ".join(sorted(mapping.keys())))
        )

    return mapping[value]


def build_inputs(data):
    inputs = SimpleNamespace()

    inputs.nstep = get_int(data, "nstep", required=True)
    inputs.dt = get_float(data, "dt", required=True)

    inputs.nstout = get_int(data, "nstout", 1000)
    inputs.nstdcd = get_int(data, "nstdcd", 1000)

    inputs.coulomb_name = get_str(data, "coulomb", "PME")
    inputs.coulomb = get_nonbonded_method(inputs.coulomb_name)

    inputs.ewald_Tol = get_float(data, "ewald_Tol", 0.0005)

    inputs.vdw = get_str(data, "vdw", "Force-switch")
    inputs.r_on = get_float(data, "r_on", 1.0)
    inputs.r_off = get_float(data, "r_off", 1.2)

    inputs.temp = get_float(data, "temp", 310.0)
    inputs.fric_coeff = get_float(data, "fric_coeff", 1.0)

    inputs.pcouple = get_str(data, "pcouple", "yes")
    inputs.p_ref = get_float(data, "p_ref", 1.0)
    inputs.p_type = get_str(data, "p_type", "membrane")

    inputs.p_XYMode_name = get_str(data, "p_XYMode", "XYIsotropic")
    inputs.p_ZMode_name = get_str(data, "p_ZMode", "ZFree")

    inputs.p_XYMode = get_membrane_xy_mode(inputs.p_XYMode_name)
    inputs.p_ZMode = get_membrane_z_mode(inputs.p_ZMode_name)

    inputs.p_tens = get_float(data, "p_tens", 0.0)
    inputs.p_tens_input = inputs.p_tens

    inputs.p_freq = get_int(data, "p_freq", 100)

    inputs.cons_name = get_str(data, "cons", "HBonds")
    inputs.cons = get_constraints(inputs.cons_name)

    inputs.rest = get_str(data, "rest", "no")
    inputs.lj_lrc = get_str(data, "lj_lrc", "no")
    inputs.e14scale = get_float(data, "e14scale", 1.0)

    return inputs


# ============================================================
# File type detection
# ============================================================

def read_file_head_bytes(filename, nbytes=8192):
    with open(filename, "rb") as f:
        return f.read(nbytes)


def head_as_text(head_bytes):
    try:
        return head_bytes.decode("utf-8", errors="ignore").lstrip()
    except Exception:
        return ""


def file_text_contains(filename, text):
    try:
        with open(filename, "r") as f:
            content = f.read()
    except Exception:
        return False

    return text in content


def looks_like_openmm_xml_state(filename):
    try:
        head = read_file_head_bytes(filename, 8192)
    except Exception:
        return False

    text = head_as_text(head)

    if text.startswith("<?xml") and "<State" in text[:2000]:
        return True

    if text.startswith("<State"):
        return True

    return False


def openmm_xml_state_has_velocities(filename):
    return file_text_contains(filename, "<Velocities>")


def detect_coordinate_file_type(filename):
    """
    coordinate_file is the only starting file, like Amber -c.

    Returns one of:
        openmm_checkpoint
        openmm_xml_state
        amber_restart
        gromacs_gro
        charmm_coordinate
        unknown
    """

    if not os.path.exists(filename):
        sys.exit("coordinate_file does not exist: %s" % filename)

    lower = filename.lower()

    if lower.endswith(".chk"):
        return "openmm_checkpoint"

    if looks_like_openmm_xml_state(filename):
        return "openmm_xml_state"

    if lower.endswith(".xml"):
        return "openmm_xml_state"

    if lower.endswith(".gro"):
        return "gromacs_gro"

    if lower.endswith((".rst7", ".rst", ".inpcrd", ".crd", ".ncrst", ".restrt")):
        return "amber_restart"

    if lower.endswith((".cor", ".pdb")):
        return "charmm_coordinate"

    return "unknown"


def extract_openmm_xml_step_and_time(filename):
    try:
        root = ET.parse(filename).getroot()
    except Exception:
        return None, None

    step_count = None
    time_ps = None

    if "stepCount" in root.attrib:
        try:
            step_count = int(root.attrib["stepCount"])
        except Exception:
            step_count = None

    if "time" in root.attrib:
        try:
            time_ps = float(root.attrib["time"])
        except Exception:
            time_ps = None

    return step_count, time_ps


# ============================================================
# Coordinate and velocity helpers
# ============================================================

def load_amber_coordinates(filename):
    """
    Load AMBER coordinates or restart.

    First tries to load positions, velocities, and box.
    If velocities are absent, falls back to positions and box.
    """

    try:
        crd = AmberInpcrdFile(
            filename,
            loadVelocities=True,
            loadBoxVectors=True,
        )
    except Exception:
        crd = AmberInpcrdFile(
            filename,
            loadBoxVectors=True,
        )

    return crd


def coordinate_has_velocities(crd):
    return hasattr(crd, "velocities") and crd.velocities is not None


def set_box_vectors_from_coordinate_object(simulation, crd):
    if hasattr(crd, "boxVectors") and crd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*crd.boxVectors)
        return True

    if hasattr(crd, "getPeriodicBoxVectors"):
        box_vectors = crd.getPeriodicBoxVectors()

        if box_vectors is not None:
            simulation.context.setPeriodicBoxVectors(*box_vectors)
            return True

    return False


def use_velocities_from_coordinate_object(simulation, crd, source_name, required=True):
    if coordinate_has_velocities(crd):
        simulation.context.setVelocities(crd.velocities)
        print("Using velocities from:", source_name)
        return True

    message = "No velocities found in: %s" % source_name

    if required:
        sys.exit(
            message
            + "\n"
            + "genvel = no was requested, so the run cannot start without velocities."
        )

    print(message)
    return False


# ============================================================
# Loading coordinate_file
# ============================================================

def load_openmm_xml_state_as_coordinate_file(simulation, filename, genvel):
    print("\nLoading coordinate_file as OpenMM XML state:", filename)

    has_vel = openmm_xml_state_has_velocities(filename)

    if (not genvel) and (not has_vel):
        sys.exit(
            "coordinate_file is an OpenMM XML state, but it has no <Velocities> block:\n"
            "%s\n"
            "genvel = no was requested, so the run cannot start without velocities."
            % filename
        )

    with open(filename, "r") as f:
        state_xml = f.read()

    state = XmlSerializer.deserialize(state_xml)
    simulation.context.setState(state)

    step_count, time_ps = extract_openmm_xml_step_and_time(filename)

    if step_count is not None:
        simulation.currentStep = step_count
        print("OpenMM XML stepCount:", step_count)

    if time_ps is not None:
        print("OpenMM XML time: %.6f ps" % time_ps)

    return {
        "kind": "openmm_xml_state",
        "coordinate_object": None,
        "velocity_source": filename,
        "has_velocities": has_vel,
        "step_count": step_count,
        "time_ps": time_ps,
    }


def load_openmm_checkpoint_as_coordinate_file(simulation, filename):
    print("\nLoading coordinate_file as OpenMM checkpoint:", filename)

    loaded_with_simulation = False

    try:
        simulation.loadCheckpoint(filename)
        loaded_with_simulation = True
    except Exception as exc:
        print(
            "WARNING: simulation.loadCheckpoint failed: %s" % exc,
            file=sys.stderr,
        )
        print(
            "Trying simulation.context.loadCheckpoint instead.",
            file=sys.stderr,
        )

        with open(filename, "rb") as f:
            simulation.context.loadCheckpoint(f.read())

    if loaded_with_simulation:
        print("Loaded checkpoint with simulation.loadCheckpoint.")
    else:
        print("Loaded checkpoint with context.loadCheckpoint.")

    return {
        "kind": "openmm_checkpoint",
        "coordinate_object": None,
        "velocity_source": filename,
        "has_velocities": True,
        "step_count": int(simulation.currentStep),
        "time_ps": None,
    }


def load_amber_restart_as_coordinate_file(simulation, filename, fftype):
    if fftype != "AMBER":
        sys.exit(
            "coordinate_file looks like an AMBER restart, but force_field is not AMBER:\n"
            "%s" % filename
        )

    print("\nLoading coordinate_file as AMBER restart or coordinate:", filename)

    crd = load_amber_coordinates(filename)

    simulation.context.setPositions(crd.positions)
    set_box_vectors_from_coordinate_object(simulation, crd)

    return {
        "kind": "amber_restart",
        "coordinate_object": crd,
        "velocity_source": filename,
        "has_velocities": coordinate_has_velocities(crd),
        "step_count": None,
        "time_ps": None,
    }


def load_gromacs_gro_as_coordinate_file(simulation, filename, fftype):
    if fftype != "GROMACS":
        sys.exit(
            "coordinate_file looks like a GROMACS .gro file, but force_field is not GROMACS:\n"
            "%s" % filename
        )

    print("\nLoading coordinate_file as GROMACS GRO:", filename)

    gro = GromacsGroFile(filename)

    simulation.context.setPositions(gro.positions)
    simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())

    return {
        "kind": "gromacs_gro",
        "coordinate_object": gro,
        "velocity_source": filename,
        "has_velocities": coordinate_has_velocities(gro),
        "step_count": None,
        "time_ps": None,
    }


def load_charmm_coordinate_as_coordinate_file(simulation, filename, fftype):
    if fftype != "CHARMM":
        sys.exit(
            "coordinate_file looks like a CHARMM or PDB coordinate file, but force_field is not CHARMM:\n"
            "%s" % filename
        )

    print("\nLoading coordinate_file as CHARMM coordinate:", filename)

    from omm_readparams import read_crd

    crd = read_crd(filename, "CHARMM")

    simulation.context.setPositions(crd.positions)

    return {
        "kind": "charmm_coordinate",
        "coordinate_object": crd,
        "velocity_source": filename,
        "has_velocities": coordinate_has_velocities(crd),
        "step_count": None,
        "time_ps": None,
    }


def load_coordinate_file_auto(simulation, filename, fftype, genvel):
    """
    Load coordinate_file automatically.

    This is the Amber -c equivalent.
    There is no input_restart.
    """

    detected = detect_coordinate_file_type(filename)

    print("\ncoordinate_file auto-detection")
    print("coordinate_file:", filename)
    print("Detected coordinate_file type:", detected)

    if detected == "openmm_checkpoint":
        return load_openmm_checkpoint_as_coordinate_file(simulation, filename)

    if detected == "openmm_xml_state":
        return load_openmm_xml_state_as_coordinate_file(simulation, filename, genvel)

    if detected == "amber_restart":
        return load_amber_restart_as_coordinate_file(simulation, filename, fftype)

    if detected == "gromacs_gro":
        return load_gromacs_gro_as_coordinate_file(simulation, filename, fftype)

    if detected == "charmm_coordinate":
        return load_charmm_coordinate_as_coordinate_file(simulation, filename, fftype)

    print("coordinate_file type was unknown. Trying fallback loading.")

    if looks_like_openmm_xml_state(filename):
        return load_openmm_xml_state_as_coordinate_file(simulation, filename, genvel)

    if fftype == "AMBER":
        try:
            return load_amber_restart_as_coordinate_file(simulation, filename, fftype)
        except Exception as exc:
            sys.exit(
                "Could not load coordinate_file as AMBER restart or coordinate:\n"
                "%s\n"
                "Original error: %s" % (filename, exc)
            )

    sys.exit(
        "Could not auto-detect coordinate_file type:\n"
        "%s\n"
        "Use AMBER rst7/inpcrd/ncrst, OpenMM XML state, OpenMM .chk, GROMACS .gro, or CHARMM coordinate."
        % filename
    )


# ============================================================
# AMBER restart writer
# ============================================================

def write_amber_restart(simulation, filename, netcdf=False, enforce_pbc=True):
    try:
        from parmed.openmm.reporters import RestartReporter
    except ImportError:
        print(
            "WARNING: ParmEd not found. Cannot write AMBER restart: %s" % filename,
            file=sys.stderr,
        )
        return None

    try:
        state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            enforcePeriodicBox=enforce_pbc,
        )

        reporter = RestartReporter(
            filename,
            1,
            write_multiple=False,
            netcdf=netcdf,
        )

        reporter.report(simulation, state)

        if hasattr(reporter, "finalize"):
            try:
                reporter.finalize()
            except Exception:
                pass

        return filename

    except Exception as exc:
        print(
            "WARNING: failed to write AMBER restart %s: %s" % (filename, exc),
            file=sys.stderr,
        )
        return None


# ============================================================
# Main
# ============================================================

parser = argparse.ArgumentParser(
    description="OpenMM production runner using one Amber-style coordinate_file."
)

parser.add_argument(
    "input_file",
    help="Single input file containing filenames and all run options.",
)

args = parser.parse_args()

cfg = read_key_value_file(args.input_file)
inputs = build_inputs(cfg)

fftype = get_str(cfg, "force_field", "AMBER").upper()

platform_request = get_str(cfg, "platform", "CUDA")
cuda_precision = get_str(cfg, "cuda_precision", "single")

topfile = get_str(cfg, "topology_file", required=True)
crdfile = get_str(cfg, "coordinate_file", required=True)

toppar = normalize_none(get_str(cfg, "toppar_file", None))
gmx_include = get_str(cfg, "gmx_include", "toppar")

continuation = get_bool(cfg, "continuation", False)
genvel = get_genvel(cfg)

output_xml = get_str(cfg, "output_xml", "step7_production.xml")
output_dcd = get_str(cfg, "output_dcd", "step7_production.dcd")
output_report = get_str(cfg, "output_report", "step7_production_report.txt")

progress_file = get_str(
    cfg,
    "progress_file",
    os.path.splitext(output_report)[0] + "_progress.out",
)

progress_stdout = get_bool(cfg, "progress_stdout", True)

checkpoint_steps = get_int(cfg, "checkpoint_steps", 0)

rst7_netcdf = get_bool(cfg, "rst7_netcdf", False)
rst7_steps = get_int(cfg, "rst7_steps", 0)

rst7_ext = ".ncrst" if rst7_netcdf else ".rst7"

final_amber_restart = normalize_none(
    get_str(cfg, "final_amber_restart", None)
)

if final_amber_restart is None:
    final_amber_restart = os.path.splitext(output_dcd)[0] + rst7_ext

restart_snapshot_prefix = normalize_none(
    get_str(cfg, "restart_snapshot_prefix", None)
)

if restart_snapshot_prefix is None:
    restart_snapshot_prefix = (
        os.path.splitext(final_amber_restart)[0]
        + "_snapshot"
        + rst7_ext
    )

rewrap_coordinates = get_bool(cfg, "rewrap_coordinates", True)

# Default is no. Reset only if explicitly requested.
reset_step_and_time = get_bool(cfg, "reset_step_and_time", False)

print("Loading parameters")
print("Input file:", args.input_file)
print("Force field:", fftype)
print("Topology file:", topfile)
print("Coordinate file:", crdfile)


# ============================================================
# Load topology only
# ============================================================

params = None
top = None
topology = None

if fftype == "CHARMM":
    if toppar is None:
        sys.exit("Error: CHARMM requires toppar_file")

    from omm_readparams import read_top, read_params, read_crd, gen_box

    top = read_top(topfile, "CHARMM")
    params = read_params(toppar)

    coord_type_for_charmm = detect_coordinate_file_type(crdfile)

    if coord_type_for_charmm == "charmm_coordinate":
        try:
            crd_for_box = read_crd(crdfile, "CHARMM")
            top = gen_box(top, crd_for_box)
        except Exception as exc:
            print(
                "WARNING: could not apply CHARMM box from coordinate_file: %s" % exc,
                file=sys.stderr,
            )

    topology = top.topology

elif fftype == "AMBER":
    top = AmberPrmtopFile(topfile)
    topology = top.topology

elif fftype == "GROMACS":
    coord_type_for_top = detect_coordinate_file_type(crdfile)

    if coord_type_for_top == "gromacs_gro":
        gro_for_box = GromacsGroFile(crdfile)

        top = GromacsTopFile(
            topfile,
            periodicBoxVectors=gro_for_box.getPeriodicBoxVectors(),
            includeDir=gmx_include,
        )
    else:
        top = GromacsTopFile(
            topfile,
            includeDir=gmx_include,
        )

    topology = top.topology

else:
    sys.exit("Error: force_field must be CHARMM, AMBER, or GROMACS")


# ============================================================
# Create system
# ============================================================

nboptions = dict(
    nonbondedMethod=inputs.coulomb,
    nonbondedCutoff=inputs.r_off * nanometers,
    constraints=inputs.cons,
    ewaldErrorTolerance=inputs.ewald_Tol,
)

if inputs.vdw == "Switch":
    nboptions["switchDistance"] = inputs.r_on * nanometers

if inputs.vdw == "LJPME":
    nboptions["nonbondedMethod"] = LJPME

if fftype == "CHARMM":
    system = top.createSystem(params, **nboptions)
else:
    system = top.createSystem(**nboptions)

if fftype == "CHARMM" and inputs.vdw == "Force-switch":
    system = vfswitch(system, top, inputs)

if inputs.lj_lrc == "yes":
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            force.setUseDispersionCorrection(True)

        if isinstance(force, CustomNonbondedForce) and force.getNumTabulatedFunctions() != 1:
            force.setUseLongRangeCorrection(True)

if inputs.e14scale != 1.0:
    nonbonded = None

    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nonbonded = force
            break

    if nonbonded is not None:
        for i in range(nonbonded.getNumExceptions()):
            atom1, atom2, chg, sig, eps = nonbonded.getExceptionParameters(i)

            nonbonded.setExceptionParameters(
                i,
                atom1,
                atom2,
                chg * inputs.e14scale,
                sig,
                eps,
            )

if inputs.pcouple == "yes":
    system = barostat(system, inputs)

if fftype == "CHARMM" and inputs.rest == "yes":
    coord_type = detect_coordinate_file_type(crdfile)

    if coord_type != "charmm_coordinate":
        sys.exit(
            "CHARMM restraints require coordinate_file to be a CHARMM coordinate or PDB file."
        )

    from omm_readparams import read_crd

    crd_for_rest = read_crd(crdfile, "CHARMM")
    system = restraints(system, crd_for_rest, inputs)


# ============================================================
# Integrator and platform
# ============================================================

integrator = LangevinIntegrator(
    inputs.temp * kelvin,
    inputs.fric_coeff / picosecond,
    inputs.dt * picoseconds,
)

enabled_platforms = [
    Platform.getPlatform(i).getName()
    for i in range(Platform.getNumPlatforms())
]

if platform_request in enabled_platforms:
    platform_name = platform_request
else:
    print("Requested platform %s not available." % platform_request, file=sys.stderr)
    print("Available platforms: %s" % ", ".join(enabled_platforms), file=sys.stderr)

    platform_name = None

    for candidate in ["CUDA", "OpenCL", "CPU"]:
        if candidate in enabled_platforms:
            platform_name = candidate
            break

    if platform_name is None:
        sys.exit("Error: no usable OpenMM platform found.")

    print("Falling back to platform %s." % platform_name, file=sys.stderr)

platform = Platform.getPlatformByName(platform_name)

print("Requested platform:", platform_request)
print("Using platform:", platform.getName())

if platform.getName() == "CUDA":
    prop = dict(CudaPrecision=cuda_precision)
else:
    prop = dict()


# ============================================================
# Create simulation
# ============================================================

simulation = Simulation(topology, system, integrator, platform, prop)


# ============================================================
# Load coordinate_file, then apply velocity policy
# ============================================================

print("\nStart mode")
print("Continuation:", "yes" if continuation else "no")
print("Generate velocities:", "yes" if genvel else "no")
print("Reset step and time:", "yes" if reset_step_and_time else "no")

start_info = load_coordinate_file_auto(
    simulation=simulation,
    filename=crdfile,
    fftype=fftype,
    genvel=genvel,
)

if genvel:
    print("\nGenerating velocities from temperature:", inputs.temp, "K")
    simulation.context.setVelocitiesToTemperature(inputs.temp * kelvin)
else:
    if start_info["kind"] in ["openmm_xml_state", "openmm_checkpoint"]:
        if start_info["has_velocities"]:
            print("\nKeeping velocities from:", crdfile)
        else:
            sys.exit(
                "genvel = no was requested, but coordinate_file has no velocities:\n"
                "%s" % crdfile
            )
    else:
        use_velocities_from_coordinate_object(
            simulation,
            start_info["coordinate_object"],
            start_info["velocity_source"],
            required=True,
        )


# ============================================================
# Optional rewrap and optional reset
# ============================================================

if rewrap_coordinates:
    print("\nRewrapping starting coordinates")
    simulation = rewrap(simulation)

if reset_step_and_time:
    simulation.currentStep = 0
    simulation.context.setTime(0.0 * picoseconds)
    print("Step and time reset to zero.")
else:
    print("Step and time were not reset.")


# ============================================================
# Initial energy
# ============================================================

print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())


# ============================================================
# Reporters and run
# ============================================================

sim_time_ps = inputs.nstep * inputs.dt
sim_time_ns = sim_time_ps / 1000.0

start_step_for_report = int(simulation.currentStep)
planned_final_step = start_step_for_report + inputs.nstep

print("\nMD run: %s steps" % inputs.nstep)
print("Starting step:", start_step_for_report)
print("Planned final step:", planned_final_step)
print("Planned simulation length: %.3f ns" % sim_time_ns)

if inputs.nstdcd > 0:
    simulation.reporters.append(
        DCDReporter(
            output_dcd,
            inputs.nstdcd,
            enforcePeriodicBox=True,
        )
    )

progress_handle = open(progress_file, "w", buffering=1)

simulation.reporters.append(
    StateDataReporter(
        progress_handle,
        inputs.nstout,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=planned_final_step,
        separator="\t",
    )
)

print("Progress file:", progress_file)

if progress_stdout:
    simulation.reporters.append(
        StateDataReporter(
            sys.stdout,
            inputs.nstout,
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=planned_final_step,
            separator="\t",
        )
    )

chk_interval = 0

if checkpoint_steps != 0:
    chk_interval = inputs.nstout if checkpoint_steps < 0 else checkpoint_steps
    chk_path = output_xml + ".chk"

    simulation.reporters.append(
        CheckpointReporter(chk_path, chk_interval)
    )

    print("Checkpoint every %d steps -> %s" % (chk_interval, chk_path))


rst7_interval = 0

if rst7_steps != 0:
    try:
        from parmed.openmm.reporters import RestartReporter

        rst7_interval = inputs.nstout if rst7_steps < 0 else rst7_steps

        simulation.reporters.append(
            RestartReporter(
                restart_snapshot_prefix,
                rst7_interval,
                write_multiple=True,
                netcdf=rst7_netcdf,
            )
        )

        print(
            "Periodic AMBER restart every %d steps -> %s.<step>"
            % (rst7_interval, restart_snapshot_prefix)
        )

    except ImportError:
        print(
            "WARNING: ParmEd not found. Periodic AMBER restarts disabled.",
            file=sys.stderr,
        )
        rst7_interval = 0


start_wall = time.time()
start_datetime = datetime.datetime.now()

try:
    simulation.step(inputs.nstep)
finally:
    progress_handle.flush()
    progress_handle.close()

end_wall = time.time()
end_datetime = datetime.datetime.now()

elapsed_seconds = end_wall - start_wall
elapsed_hours = elapsed_seconds / 3600.0
elapsed_days = elapsed_seconds / 86400.0

if elapsed_days > 0:
    final_speed_ns_per_day = sim_time_ns / elapsed_days
else:
    final_speed_ns_per_day = 0.0


# ============================================================
# Write final restarts
# ============================================================

if rewrap_coordinates:
    print("\nRewrapping final coordinates before writing restart")
    simulation = rewrap(simulation)

state = simulation.context.getState(
    getPositions=True,
    getVelocities=True,
)

with open(output_xml, "w") as f:
    f.write(XmlSerializer.serialize(state))

print("Final OpenMM XML restart written to:", output_xml)

print("\nWriting final AMBER restart")

written_rst7 = write_amber_restart(
    simulation,
    final_amber_restart,
    netcdf=rst7_netcdf,
    enforce_pbc=True,
)

if written_rst7:
    print("Final AMBER restart written to:", written_rst7)


# ============================================================
# Write report
# ============================================================

with open(output_report, "w") as f:
    f.write("OpenMM production report\n")
    f.write("========================\n\n")

    f.write("Working directory: %s\n" % os.getcwd())
    f.write("Input file: %s\n" % args.input_file)

    f.write("\nFiles\n")
    f.write("-----\n")
    f.write("Force field: %s\n" % fftype)
    f.write("Topology file: %s\n" % topfile)
    f.write("Coordinate file: %s\n" % crdfile)
    f.write("Detected coordinate_file type: %s\n" % start_info["kind"])
    f.write("Output XML restart: %s\n" % output_xml)
    f.write("Output trajectory: %s\n" % output_dcd)
    f.write("Progress file: %s\n" % progress_file)
    f.write("Final AMBER restart: %s\n" % (written_rst7 if written_rst7 else "FAILED"))
    f.write("Output report: %s\n" % output_report)

    f.write("\nStart mode\n")
    f.write("----------\n")
    f.write("Continuation: %s\n" % ("yes" if continuation else "no"))
    f.write("Generate velocities: %s\n" % ("yes" if genvel else "no"))
    f.write("Rewrap coordinates: %s\n" % ("yes" if rewrap_coordinates else "no"))
    f.write("Reset step and time: %s\n" % ("yes" if reset_step_and_time else "no"))
    f.write("Coordinate velocity source: %s\n" % start_info["velocity_source"])
    f.write("Coordinate file has velocities: %s\n" % ("yes" if start_info["has_velocities"] else "no"))

    if start_info["step_count"] is not None:
        f.write("Coordinate file step count: %s\n" % start_info["step_count"])

    if start_info["time_ps"] is not None:
        f.write("Coordinate file time: %.6f ps\n" % start_info["time_ps"])

    f.write("\nRun options\n")
    f.write("-----------\n")
    f.write("Platform requested: %s\n" % platform_request)
    f.write("Platform used: %s\n" % platform.getName())
    f.write("CUDA precision: %s\n" % cuda_precision)
    f.write("Available platforms: %s\n" % ", ".join(enabled_platforms))
    f.write("Progress also to stdout: %s\n" % ("yes" if progress_stdout else "no"))
    f.write("Checkpoint every N steps: %s\n" % (chk_interval if chk_interval > 0 else "off"))
    f.write("Periodic rst7 every N steps: %s\n" % (rst7_interval if rst7_interval > 0 else "off"))
    f.write("Periodic rst7 prefix: %s\n" % restart_snapshot_prefix)
    f.write("rst7 format: %s\n" % ("NetCDF" if rst7_netcdf else "ASCII"))

    f.write("\nMD settings\n")
    f.write("-----------\n")
    f.write("nstep: %d\n" % inputs.nstep)
    f.write("dt: %.6f ps\n" % inputs.dt)
    f.write("nstout: %d\n" % inputs.nstout)
    f.write("nstdcd: %d\n" % inputs.nstdcd)
    f.write("Simulation length: %.3f ps\n" % sim_time_ps)
    f.write("Simulation length: %.3f ns\n" % sim_time_ns)
    f.write("Temperature: %.3f K\n" % inputs.temp)
    f.write("Friction coefficient: %.3f /ps\n" % inputs.fric_coeff)
    f.write("Coulomb method: %s\n" % inputs.coulomb_name)
    f.write("vdW method: %s\n" % inputs.vdw)
    f.write("r_on: %.3f nm\n" % inputs.r_on)
    f.write("r_off: %.3f nm\n" % inputs.r_off)
    f.write("constraints: %s\n" % inputs.cons_name)
    f.write("pcouple: %s\n" % inputs.pcouple)
    f.write("p_type: %s\n" % inputs.p_type)
    f.write("p_ref: %.3f bar\n" % inputs.p_ref)
    f.write("p_XYMode: %s\n" % inputs.p_XYMode_name)
    f.write("p_ZMode: %s\n" % inputs.p_ZMode_name)
    f.write("p_tens: %.3f\n" % inputs.p_tens_input)
    f.write("p_freq: %d\n" % inputs.p_freq)

    f.write("\nStep accounting\n")
    f.write("---------------\n")
    f.write("Start step: %d\n" % start_step_for_report)
    f.write("Planned final step: %d\n" % planned_final_step)
    f.write("Actual final step: %d\n" % int(simulation.currentStep))

    f.write("\nTiming\n")
    f.write("------\n")
    f.write("Start time: %s\n" % start_datetime)
    f.write("End time: %s\n" % end_datetime)
    f.write("Wall time: %.2f seconds\n" % elapsed_seconds)
    f.write("Wall time: %.3f hours\n" % elapsed_hours)
    f.write("Average performance: %.3f ns/day\n" % final_speed_ns_per_day)


print("Production completed successfully.")
print("Report written to:", output_report)
print("Final average performance: %.3f ns/day" % final_speed_ns_per_day)
