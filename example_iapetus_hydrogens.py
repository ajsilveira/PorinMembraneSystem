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

# File from MemProt Data Base
pdb = PDBFile('atomistic-system.pdb')

modeller = Modeller(pdb.topology,pdb.positions)
modeller.changeAtomNames()
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller.addHydrogens(forcefield=forcefield)
PDBxFile.writeFile(modeller.topology,modeller.positions,open("atomistic-system-with-hydrogens.pdbx", 'w'))

system = forcefield.createSystem(modeller.topology,
                                 nonbondedMethod=app.PME,
                                 rigidWater=True,
                                 nonbondedCutoff=1*unit.nanometer)
integrator = mm.VerletIntegrator(0.5*unit.femtoseconds)
platform = mm.Platform.getPlatformByName('CUDA')
simulation = app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
# Minimize the system after adding hydrogens to the membrane
tolerance = 0.1*unit.kilojoules_per_mole/unit.angstroms
simulation.minimizeEnergy(tolerance=tolerance,maxIterations=0)
simulation.reporters.append(app.StateDataReporter('relax-hydrogens.log', 1000, step=True,
                                                           temperature=True,
                                                           potentialEnergy=True,
                                                           totalEnergy=True,
                                                           speed=True))
simulation.step(10000)
positions = simulation.context.getState(getPositions=True).getPositions()
min_potential = atomsmm.splitPotentialEnergy(system, modeller.topology, positions)
for key, value in min_potential.items():
    logger.debug(key, value)
logger.debug(modeller.topology)
simulation.saveCheckpoint('state.chk')
#del simulation.context, integrator
# rigidWater = False is necessary because water parameters are required by Parmed
# when loading the ligand
#system = forcefield.createSystem(modeller.topology,
#                                 nonbondedMethod=app.PME,
#                                 rigidWater=False,
#                                 nonbondedCutoff=1*unit.nanometer)
#ligand_atomistic_system = PorinMembraneSystem('comp7', system, modeller.topology, positions, platform, membrane = 'DPPC')
