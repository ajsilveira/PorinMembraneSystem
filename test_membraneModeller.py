"""
Test for the membraneModeller module.
"""

# Import packages
import simtk.openmm as mm
from simtk.openmm import unit
import simtk.openmm.app as app
from membraneModeller import Modeller
from simtk.openmm.app import PDBFile, ForceField

def test_membraneModeller():
    # pdb file containing a solvated single lipid molecule
    pdb = PDBFile('solvated-lipid.pdb')
    modeller = Modeller(pdb.topology,pdb.positions)
    modeller.modifyTopology()
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller.addHydrogens(forcefield=forcefield)

    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME,
                                     rigidWater=True,
                                     nonbondedCutoff=1*unit.nanometer)
    integrator = mm.VerletIntegrator(0.5*unit.femtoseconds)
    platform = mm.Platform.getPlatformByName('Reference')
    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    # Minimize the system after adding hydrogens
    tolerance = 0.1*unit.kilojoules_per_mole/unit.angstroms
    simulation.minimizeEnergy(tolerance=tolerance,maxIterations=200)
    simulation.reporters.append(app.StateDataReporter('relax-hydrogens.log',
                                                       1000,
                                                       step=True,
                                                       potentialEnergy=True))
    simulation.step(1000)
