"""
This module contains some wrappers for using ASE to perform
various operations on BigDFT molecules.

https://wiki.fysik.dtu.dk/ase/
"""
from ase.calculators.calculator import Calculator


def bigdft_to_ase(sys):
    """
    Given a BigDFT system, this transforms that system into an ASE
    collection of atoms.

    Args:
      frag (BigDFT.Systems.System): the system to convert.

    Returns:
      (ase.Atoms): a collection of atoms for ASE.
    """
    from ase import Atom, Atoms
    from ase.units import Bohr

    atlist = []
    for frag in sys.values():
        for at in frag:
            pos = at.get_position("angstroem")
            atlist += [Atom(at.sym, (pos[0], pos[1], pos[2]))]

    cellargs = {}
    if sys.cell is not None and sys.cell.get_boundary_condition() != "free":
        cellargs["cell"] = []
        cellargs["pbc"] = []

        for i in range(3):
            if sys.cell[i, i] == float("inf"):
                cellargs["cell"].append([0, 0, 0])
                cellargs["pbc"].append(False)
            else:
                cellargs["cell"].append([x*Bohr for x in sys.cell[i, :]])
                cellargs["pbc"].append(True)

    return Atoms(atlist, **cellargs)


def ase_to_bigdft(mol):
    """
    Given an ASE collection of atoms, this transforms that collection into
    a BigDFT fragment.

    Args:
      mol (ase.Atoms): a collection of atoms used by ASE.

    Returns:
      (BigDFT.Fragments.Fragment): a BigDFT fragment.
    """
    from BigDFT.Fragments import Fragment
    from BigDFT.Atoms import Atom

    frag = Fragment()
    for at in mol:
        frag += [Atom({"sym": at.symbol, "r": at.position,
                       "units": "angstroem"})]

    return frag


def ase_potential_energy(sys, ase_calculator):
    """
    Given a ASE calculator, calculates the potential energy
    of a BigDFT System.

    Args:
        sys (BigDFT.Systems.System): a BigDFT system
        ase_calculator (ase.calculators.calculator.Calculator): ASE calculator

    Returns:
        float: the potential energy, in Hartree
    """
    from ase.units import Hartree
    asys = bigdft_to_ase(sys)
    asys.set_calculator(ase_calculator)
    return asys.get_potential_energy() / Hartree


class BigASECalculator(Calculator):
    """
    This calculator can be used by ASE to run BigDFT calculations.

    Args:
      parameters (BigDFT.Inputfiles.Inputfile): the input parameters to use.
      calc (BigDFT.Calculators.SystemCalculator): the calculator to use.
      directory (str): the directory to run the calculation in.
      label (str): a label for this particular calculation.
      atoms (ase.Atoms): an ase collection of atoms to compute.
    """
    def __init__(self, parameters, calc, restart=None,
                 ignore_bad_restart_file=False, directory=".", label=None,
                 atoms=None, **kwargs):

        # Call the Super Method
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

        # Internal data.
        self.log = None
        self.code = calc
        self.parameters = parameters
        self.rundir = directory

    def calculation_required(self, atoms, quantities):
        """
        https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/interface.html
        """
        from copy import deepcopy

        if atoms != self.atoms:
            self.atoms = deepcopy(atoms)
            return True
        elif self.log is None:
            return True
        else:
            return False

    def get_forces(self, atoms):
        """
        https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/interface.html
        """
        from numpy import array
        from ase.units import Bohr, Hartree

        if self.calculation_required(atoms, quantities=["forces"]):
            self._run_calculation()
        sys = self._create_system()
        sys.set_atom_forces(self.log)

        forces = []
        for at in sys["FRAG:0"]:
            forces.append(at.get_force())

        self.results["forces"] = array(forces) * (Hartree / Bohr)
        return self.results["forces"]

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """
        https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/interface.html
        """
        from ase.units import Hartree
        if self.calculation_required(atoms, quantities=["energy"]):
            self._run_calculation()
        self.results["energy"] = self.log.energy * Hartree
        return self.results["energy"]

    def _create_system(self):
        """
        Create a BigDFT system from the atoms of this calculator.
        """
        from BigDFT.Systems import System
        from BigDFT.UnitCells import UnitCell
        from ase.units import Bohr

        sys = System()
        sys["FRAG:0"] = ase_to_bigdft(self.atoms)

        if self.atoms.cell is not None:
            if any(self.atoms.get_pbc()):
                if not self.atoms.cell.orthorhombic:
                    raise ValueError("Only orthorhombic cells are supported.")
                vec = [float(self.atoms.cell[i, i]) if self.atoms.get_pbc()[i]
                       else float("inf") for i in range(3)]
                sys.cell = UnitCell(vec, units="angstroem")

        return sys

    def _get_posinp(self):
        """
        Get the posinp associated with the atoms of this calculator.
        """
        return self._create_system().get_posinp()

    def _run_calculation(self):
        self.log = self.code.run(input=self.parameters,
                                 posinp=self._get_posinp(),
                                 name=self.label, run_dir=self.rundir)


def _example():
    """Example of using ASE interoperability"""
    from BigDFT.IO import XYZReader
    from BigDFT.Inputfiles import Inputfile
    from BigDFT.Calculators import SystemCalculator
    from BigDFT.Fragments import Fragment
    from BigDFT.Systems import System
    from ase.calculators.lj import LennardJones
    from BigDFT.UnitCells import UnitCell
    from ase.units import Hartree
    from ase.optimize import BFGS
    from os.path import join
    from copy import deepcopy

    # Create a system.
    sys = System()
    reader = XYZReader("Ar")
    sys["FRA:1"] = Fragment(xyzfile=reader)
    sys["FRA:2"] = deepcopy(sys["FRA:1"])
    sys["FRA:2"].translate([-2, 0, 0])

    # Skip straight to the potential energy.
    print(ase_potential_energy(sys, LennardJones()))

    # Advanced used.
    asys = bigdft_to_ase(sys)
    asys.set_calculator(LennardJones())
    dyn = BFGS(asys)
    dyn.run(fmax=0.05)
    print(asys.get_potential_energy() / Hartree)

    # Unit cells
    sys.cell = UnitCell([5, 5, 5], units="bohr")
    print(ase_potential_energy(sys, LennardJones()))
    sys.cell = UnitCell([5, float("inf"), 5], units="bohr")
    print(ase_potential_energy(sys, LennardJones()))
    sys.cell = UnitCell([float("inf"), float("inf"), 5], units="bohr")
    print(ase_potential_energy(sys, LennardJones()))

    # We can also use BigDFT with ase.
    inp = Inputfile()
    inp.set_xc("PBE")
    inp.set_hgrid(0.4)
    code = SystemCalculator(verbose=False)
    sys.cell = UnitCell()
    print(ase_potential_energy(sys, BigASECalculator(inp, code,
                                                     directory="work",
                                                     label="ase-free")))
    sys.cell = UnitCell([5, 5, 5], units="bohr")
    print(ase_potential_energy(sys, BigASECalculator(inp, code,
                                                     directory="work",
                                                     label="ase-periodic")))
    sys.cell = UnitCell([5, float("inf"), 5], units="bohr")
    print(ase_potential_energy(sys, BigASECalculator(inp, code,
                                                     directory="work",
                                                     label="ase-surface")))
    sys.cell = UnitCell([float("inf"), float("inf"), 5], units="bohr")
    print(ase_potential_energy(sys, BigASECalculator(inp, code,
                                                     directory="work",
                                                     label="ase-wire")))


if __name__ == "__main__":
    _example()
