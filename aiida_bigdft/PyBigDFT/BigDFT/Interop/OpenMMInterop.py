"""
This module contains some wrappers for using OpenMM to perform
various operations on BigDFT molecules.

https://simtk.org/api_docs/openmm/api6_0/python/index.html
"""

from BigDFT.Systems import System


def get_available_ff_names():
    import os
    import sys
    names = []
    for path in sys.path:
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(".xml"):
                    filename = os.path.join(root, file)
                    if 'data' in filename and 'openmm' in filename:
                        iname = filename.index('data')+5
                        names.append(filename[iname:])
    return names


class OMMSystem(System):
    """
    Class of OpenMM binder system which enables functionalities
    of openMM on a BigDFT system

    Args:
        system (System): a instance of a system class
        filename (str): name of the PDB file to instantiate the system
    """
    def __init__(self, system=None, filename=None):
        from simtk.openmm.app.pdbfile import PDBFile
        from simtk.openmm import app
        from BigDFT.IO import read_pdb, write_pdb
        from tempfile import NamedTemporaryFile as tmp
        if filename is not None:
            pdb = PDBFile(open(filename))
            sys = read_pdb(open(filename))
        elif system is not None:
            sys = system
            ofile = tmp('w+')
            write_pdb(system=system, ofile=ofile)
            ofilename = ofile.name
            pdb = PDBFile(open(ofilename))
            ofile.close()
        System.__init__(self, **sys.dict())
        self.pdb = pdb
        self.modeller = app.Modeller(pdb.topology, pdb.positions)

    def set_forcefields(self, *ff_list):
        """
        Define a set of force fields that will be used in the

        Args:
            *ff_list: list of the force fields to be included,
                in priority order.
                Use the :func:py:`get_available_ff_names` to identify the
                available force fields
        """
        from simtk.openmm import app
        self.forcefield = app.ForceField(*ff_list)

    @property
    def OMMsystem(self):
        from simtk.openmm import app
        if not hasattr(self, '_system'):
            self._system = self.forcefield.createSystem(
                self.modeller.topology, nonbondedMethod=app.NoCutoff,
                constraints=None)
        return self._system

    def set_integrator(self, T=298.15):
        from simtk import unit as u
        import simtk.openmm as mm
        temperature = T * u.kelvin
        self.integrator = mm.LangevinIntegrator(
            temperature, 1 / u.picosecond,  0.0005 * u.picoseconds)

    @property
    def OMMsimulation(self):
        from simtk.openmm import app
        if not hasattr(self, '_simulation'):
            self._simulation = app.Simulation(self.modeller.topology,
                                              self.OMMsystem, self.integrator)
            self._simulation.context.setPositions(self.modeller.positions)
        return self._simulation

    def OMMenergy(self, units='kcal/mol'):
        from simtk.openmm import KcalPerKJ
        energy = self.OMMsimulation.context.getState(
             getEnergy=True).getPotentialEnergy()
        return energy._value * KcalPerKJ

    @property
    def OMMposition(self):
        return self.OMMsimulation.context.getState(
            getPositions=True).getPositions()

    def write(self, ofile):
        from simtk.openmm import app
        app.PDBFile.writeFile(self.OMMsimulation.topology, self.OMMposition,
                              open(ofile, 'w'))

    def optimize(self, iters):
        self.OMMsimulation.minimizeEnergy(maxIterations=iters)

    @property
    def forcegroups(self):
        forcegroups = {}
        for i in range(self.OMMsystem.getNumForces()):
            force = self.OMMsystem.getForce(i)
            force.setForceGroup(i)
            forcegroups[force] = i
        return forcegroups

    def get_energies(self):
        from simtk.openmm import KcalPerKJ

        def component_name(k):
            return str(k).split('.')[-1].split(';')[0].lstrip('"')

        energies = {}
        for f, i in self.forcegroups.items():
            en = self.OMMsimulation.context.getState(
                getEnergy=True, groups=2**i).getPotentialEnergy()
            name = component_name(f)
            energies[name] = en._value * KcalPerKJ
        return energies
