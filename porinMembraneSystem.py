import atomsmm
import logging
import openmmtools
import mdtraj as md
import parmed as pmd
import simtk.unit as unit
from copy import deepcopy
import simtk.openmm as mm
import simtk.openmm.app as app
from yank.experiment import ExperimentBuilder
from mdtraj.reporters import NetCDFReporter
from simtk.openmm.app import Topology, PDBFile, ForceField, PDBxFile
class PorinMembraneSystem(object):

    def __init__(self, ligand_name, system, topology, positions, platform, membrane=None):
        self.membrane = membrane
        self.logger = logging.getLogger(__name__)
        with open(ligand_name + '.prmtop') as f:
            res_label = False
            for line in f:
                if (res_label):
                    self.lig_name = str(next(f).split()[0])
                    break
                if line.startswith('%'+'FLAG RESIDUE_LABEL'):
                    res_label = True

        self.structure =  pmd.openmm.load_topology(topology,
                                                   system=system,
                                                   xyz=positions)
        self.ligand = pmd.load_file(ligand_name + '.prmtop', xyz=ligand_name + '.inpcrd')
        topo = md.Topology.from_openmm(topology)
        if membrane is not None:
            porin_indices = topo.select('(protein and not resname ' + membrane + ')')
        else:
            porin_indices = topo.select('protein')

        self.structure = self._PlaceLigand(self.structure, porin_indices)
        #indices_test = []
        #for res in self.structure.residues:
         #   if res.name != 'CL' and res.name !='HOH' and res.name !='NA':
         #       atoms = res.atoms
         #       for item in atoms:
         #           indices_test.append(item.idx)
       
        new_topo = md.Topology.from_openmm(self.structure.topology)
        self.atom_to_freeze = new_topo.select('protein or resname  ' + membrane + ' or resname ' + self.lig_name)

        self._minimizeOverlaps(platform, membrane)
    def _PlaceLigand(self, structure, porin_indices):

        coords = structure.get_coordinates(frame=0)
        coords_ligand = self.ligand.get_coordinates(frame=0)
        coords_ligand += coords[porin_indices,:].mean(0) - coords_ligand[:,:].mean(0)
        coords_ligand[:,2] += coords[porin_indices,:].min(0)[2] - coords[porin_indices,:].mean(0)[2]
        self.ligand.coordinates = coords_ligand[:,:]
        new_structure = structure + self.ligand

        return new_structure

    def _minimizeOverlaps(self, platform, membrane):

        hydrogen_mass = 3*unit.amu
        system = self.structure.createSystem(nonbondedMethod=app.PME,
                                 nonbondedCutoff=1*unit.nanometer,
                                 rigidWater=True,
                                 flexibleConstraints=False,
                                 constraints=app.HBonds,
                                 hydrogenMass=hydrogen_mass,
                                 removeCMMotion=False)
        backup_system = deepcopy(system)

        for index in self.atom_to_freeze:
            system.setParticleMass(int(index), 0.0*unit.dalton)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2*unit.femtosecond)
        simulation = app.Simulation(self.structure.topology, system, integrator, platform)
        simulation.context.setPositions(self.structure.positions)
        simulation.minimizeEnergy(tolerance=0.1*unit.kilojoule/unit.mole, maxIterations=0)
        positions = simulation.context.getState(getPositions=True).getPositions()
        min_potential = atomsmm.splitPotentialEnergy(system, self.structure.topology, positions)
        for key, value in min_potential.items():
            self.logger.debug(key, value)
        simulation.saveCheckpoint('porinMembraneSystem.chk')
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2*unit.femtosecond)
        simulation = app.Simulation(self.structure.topology, backup_system, integrator, platform)
        simulation.context.setPositions(positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.reporters.append(app.StateDataReporter('relax-ligand-system.log', 100, step=True,
                                                           temperature=True,
                                                           potentialEnergy=True,
                                                           totalEnergy=True,
                                                           speed=True,
                                                          ))
        topology = md.Topology.from_openmm(self.structure.topology)
        atom_indices = topology.select('(protein and not resname  ' + membrane + ')' + ' or resname ' + self.lig_name)
        simulation.reporters.append(NetCDFReporter('relax-ligand-system.nc', 100, atomSubset=atom_indices))

        simulation.step(100000)
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBxFile.writeFile(self.structure.topology,positions,open("atomistic-system-combined.pdbx", 'w'))
