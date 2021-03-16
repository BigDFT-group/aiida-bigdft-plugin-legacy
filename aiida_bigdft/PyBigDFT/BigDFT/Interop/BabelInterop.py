"""
This module contains some wrappers for using OpenBabel to perform
various operations on BigDFT molecules.

https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python.html
"""

# openbabel's forcefields can have different energy units.
_energy_conversion = {"kJ/mol": 0.00038, "kcal/mol": 0.0016}


def convert_system_to_babel(sys):
    """
    Convert a BigDFT system to an open babel molecule.

    Args:
      sys (BigDFT.Systems.System): the system to convert.

    Returns:
      (openbabel.OBMol): an open babel type molecule.
    """
    from BigDFT.IO import write_pdb
    from openbabel.openbabel import OBMol, OBConversion
    # py2 workaround
    from sys import version_info
    if version_info[0] < 3:
        from io import BytesIO as StringIO
    else:
        try:
            from io import StringIO
        except ImportError:
            from StringIO import StringIO

    # We convert by way of pdb file.
    conv = OBConversion()
    conv.SetInFormat("pdb")

    sval = StringIO()
    write_pdb(sys, sval)

    mol = OBMol()
    conv.ReadString(mol, sval.getvalue())

    return mol


def convert_babel_to_system(mol):
    """
    Convert a BigDFT fragment to an open babel molecule.

    Args:
      mol (openbabel.OBMol): the molecule to convert.

    Returns:
      (BigDFT.Systems.System): bigdft system.
    """
    from BigDFT.IO import read_mol2
    from openbabel.openbabel import OBConversion
    # py2 workaround
    from sys import version_info
    if version_info[0] < 3:
        from io import BytesIO as StringIO
    else:
        try:
            from io import StringIO
        except ImportError:
            from StringIO import StringIO

    conv = OBConversion()
    conv.SetOutFormat("mol2")

    sval = StringIO(conv.WriteString(mol))
    return read_mol2(sval)


def compute_smiles(sys):
    """
    Computes the SMILES representation of a given system.

    Args:
      sys (BigDFT.System.Systems): the system to compute the
        representation of.

    Return:
      (str): the smiles representation of this molecule.
    """
    from openbabel.openbabel import OBConversion

    conv = OBConversion()
    mol = convert_system_to_babel(sys)
    conv.SetOutFormat("SMI")

    retstr = conv.WriteString(mol)
    retstr = retstr.replace("\n", "")
    retstr = retstr.replace("\t", "")

    return retstr


def compute_fingerprint(sys, fps="fp2"):
    """
    Computes the fingerprint for a particular fragment.

    Args:
      sys (BigDFT.Systems.System): the fragment to compute the
        representation of.
      fps (str): the type of finger print to compute.

    Return:
      (openbabel.OBFingerprint): a fingerprint for this fragment.
    """
    from pybel import Molecule

    mol = convert_system_to_babel(sys)
    pybelmol = Molecule(mol)

    return pybelmol.calcfp(fps)


def system_energy(sys, forcefield="MMFF94", verbose=False):
    """
    Compute the energy of a system using an openbabel forcefield.

    Args:
      sys (BigDFT.Systems.System): the system to compute.
      forcefield (str): the type of forcefield to use.
      verbose (bool): whether to have openbabel run in verbose mode.

    Returns:
      (float): the energy value computed in Hartree.
    """
    from openbabel.openbabel import OBForceField, OBFF_LOGLVL_LOW

    # Setup the forcefield
    ff = OBForceField.FindForceField(forcefield)
    mol = convert_system_to_babel(sys)
    ff.Setup(mol)

    if verbose:
        ff.SetLogToStdOut()
        ff.SetLogLevel(OBFF_LOGLVL_LOW)

    # Call the energy routine.
    return ff.Energy() * _energy_conversion[ff.GetUnit()]


def optimize_system(sys, forcefield="MMFF94", method="SteepestDescent",
                    steps=1000, econv=1e-6, verbose=False):
    """
    Optimize the geometry of a given fragment.

    Args:
      sys (BigDFT.Systems.System): the fragment to optimize.
      forcefield (str): the type of forcefield to use.
      verbose (bool): if True, the openbabel output will be printed.

    Returns:
      (BigDFT.Fragments.Fragment): a new fragment with the optimized positions.
    """
    from openbabel.openbabel import OBForceField, OBFF_LOGLVL_LOW
    from copy import deepcopy

    # Setup the forcefield
    ff = OBForceField.FindForceField(forcefield)
    mol = convert_system_to_babel(sys)
    ff.Setup(mol)
    if verbose:
        ff.SetLogToStdOut()
        ff.SetLogLevel(OBFF_LOGLVL_LOW)

    # Call the optimization routine.
    if method == "SteepestDescent":
        ff.SteepestDescent(steps, econv)
    elif method == "ConjugateGradients":
        ff.ConjugateGradients(steps, econv)
    else:
        raise ValueError("Invalid minimization method.")

    # Extract out the the positions.
    ff.GetCoordinates(mol)
    newsys = convert_babel_to_system(mol)

    return newsys


def molecular_dynamics(sys, steps, temperature, forcefield="MMFF94",
                       timestep=0.001, verbose=False):
    """
    Run molecular dynamics on a given fragment..

    Args:
      sys (BigDFT.Systemtems.System): the system to run.
      steps (int): the number of MD steps to take.
      temperature (float): temperature in K.
      timestep (float): time step in picoseconds.
      forcefield (str): the type of forcefield to use.
      constraints (list): for each atom, list whether it if frozen or not.
      verbose (bool): if True, the openbabel output will be printed.

    Returns:
      (BigDFT.Systems.System): a new system with the optimized positions.
    """
    from openbabel.openbabel import OBForceField, OBFF_LOGLVL_LOW

    # Setup the calculation
    ff = OBForceField.FindForceField(forcefield)
    mol = convert_system_to_babel(sys)
    ff.Setup(mol)
    if verbose:
        ff.SetLogToStdOut()
        ff.SetLogLevel(OBFF_LOGLVL_LOW)

    ff.MolecularDynamicsTakeNSteps(steps, temperature, timestep)

    ff.GetCoordinates(mol)
    newsys = convert_babel_to_system(mol)

    return newsys


def compute_system_forces(sys, forcefield="MMFF94", verbose=False):
    """
    Assign the forces of a system using an openbabel forcefield.

    Args:
      sys (BigDFT.Systems.System): the system to compute.
      forcefield (str): the type of forcefield to use.
      verbose (bool): whether to have openbabel run in verbose mode.

    Returns:
      (float): the energy of the system.
    """
    from openbabel.openbabel import OBForceField, OBFF_LOGLVL_LOW
    from BigDFT.Atoms import Atom, number_to_symbol, AU_to_A

    # Handle verbository.
    if verbose:
        ff.SetLogToStdOut()
        ff.SetLogLevel(OBFF_LOGLVL_LOW)

    # Setup the forcefield
    ff = OBForceField.FindForceField(forcefield)
    mol = convert_system_to_babel(sys)
    ff.Setup(mol)

    # Compute the forces.
    energy_out = ff.Energy() * _energy_conversion[ff.GetUnit()]
    gradients = []
    for idx in range(1, mol.NumAtoms()+1):
        at = mol.GetAtom(idx)
        grad = ff.GetGradient(at)
        convgrad = [grad.GetX(), grad.GetY(), grad.GetZ()]
        convgrad = [x * _energy_conversion[ff.GetUnit()]/AU_to_A
                    for x in convgrad]
        gradients.append(convgrad)

    # Create the atom list for the compute matching procedure.
    atom_list = []
    for idx in range(1, mol.NumAtoms()+1):
        obat = mol.GetAtom(idx)
        atnum = obat.GetAtomicNum()
        pos = [obat.GetX(), obat.GetY(), obat.GetZ()]
        atom_list.append(Atom({number_to_symbol(atnum): pos,
                               "units": "angstroem"}))
    lookup = sys.compute_matching(atom_list)

    # Assign
    for fragid, frag in sys.items():
        for i, at in enumerate(frag):
            idx = lookup[fragid][i]
            at.set_force(gradients[idx])

    return energy_out


def _setup_constraints(forcefield, constraints):
    """
    This helper routine takes a list of constraints and updates
    as forcefield with those values.
    """
    constr = forcefield.GetConstraints()
    for i in range(0, len(constraints)):
        if constraints[i] is None:
            continue
        elif constraints[i] == "fxyz" or constraints[i] == "f":
            constr.AddAtomConstraint(i+1)
        elif constraints[i] == "fx":
            constr.AddAtomXConstraint(i+1)
        elif constraints[i] == "fy":
            constr.AddAtomYConstraint(i+1)
        elif constraints[i] == "fz":
            constr.AddAtomZConstraint(i+1)
    forcefield.SetConstraints(constr)


def _example():
    """Example of using OpenBabel interoperability"""
    from BigDFT.Systems import System
    from BigDFT.Fragments import Fragment
    from BigDFT.IO import XYZReader
    from os.path import join

    # Read in a system.
    sys = System()
    sys["FRA:1"] = Fragment()
    with XYZReader("CH4") as ifile:
        for at in ifile:
            sys["FRA:1"] += Fragment([at])

    # We can compute the smiles representation.
    print(compute_smiles(sys))

    # The energy.
    print(system_energy(sys, forcefield="UFF"))

    # Extract the forces.
    compute_system_forces(sys, forcefield="UFF")
    for frag in sys.values():
        for at in frag:
            print(at["force"])

    # Optimize the geometry.
    sys2 = optimize_system(sys, forcefield="UFF")
    print(system_energy(sys2, forcefield="UFF"))


if __name__ == "__main__":
    _example()
