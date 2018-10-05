import simtk.openmm as mm
import membraneCompatibility
import simtk.openmm.app as app
from simtk.openmm import unit
from membraneCompatibility import Modeller
from simtk.openmm.app import Topology, PDBFile, ForceField, PDBxFile
import atomsmm
import logging
import yank
import mdtraj as md

import porinMembraneSystem
from porinMembraneSystem import PorinMembraneSystem

logger = logging.getLogger(__name__)
# Setup general logging
logging.root.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.DEBUG)
yank.utils.config_root_logger(verbose=True, log_file_path=None)

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

pdbx = PDBxFile('atomistic-system-with-hydrogens.pdbx')

system = forcefield.createSystem(pdbx.topology,
                                 nonbondedMethod=app.PME,
                                 rigidWater=False,
                                 nonbondedCutoff=1*unit.nanometer)

integrator = mm.VerletIntegrator(0.5*unit.femtosecond)
platform = mm.Platform.getPlatformByName('CUDA')
simulation = app.Simulation(pdbx.topology, system, integrator, platform)
simulation.loadCheckpoint('state.chk')
positions = simulation.context.getState(getPositions=True).getPositions()

#rigidWater = False is necessary because water parameters are required by Parmed
#when loading the ligand
ligand_atomistic_system = PorinMembraneSystem('comp7', system, pdbx.topology, positions, platform, membrane = 'DPPC')
