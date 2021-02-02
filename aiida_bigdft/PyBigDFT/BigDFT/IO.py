"""
This module defines some tools used for reading/writing fragment files.

"""
from futile.Utils import write as safe_print

# Standard residue names
# https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf
_standard_residues = ["ALA", "ACD", "ALA", "ALB", "ALI", "ARG", "ARO", "ASN",
                      "ASP", "ASX", "BAS", "CYS", "GLN", "GLU", "GLX", "GLY",
                      "HIS", "HYP", "ILE", "LEU", "LYS", "MET", "PCA", "PHE",
                      "PRO", "SAR", "SER", "THR", "TRP", "TYR", "VAL"]


def _split_line(line, limit):
    split = line.split()
    # print line, split
    if len(split) < 3:  # there should be at least two atoms
        tmp = line.lstrip('CONECT')
        split = ['CONECT'] + _split_malformed_connection(tmp)
    resplit = [split[0]]
    for at2 in split[1:]:
        atom = int(at2)
        if atom > limit:
            atoms = _split_malformed_connection(at2)
        else:
            atoms = [atom]
        resplit += atoms
    return resplit


def _split_malformed_connection(conect):
    maxl = len(conect)
    if maxl > 10:  # we have three atoms at least
        ll = 5 if conect[0] == '9' else 6
        return [conect[:ll]] + _split_malformed_connection(conect[ll:])
    else:
        ll = int(maxl/2)
        tmp1 = conect[:ll]
        tmp2 = conect[ll:]
        return [tmp1, tmp2]


def read_pdb(ifile, include_chain=False, disable_warnings=False):
    """
    Read a system from a PDB file.

    Args:
      ifile (TextIOBase): the file to read from.
      disable_warnings (bool): whether to print warnings about possible file
        issues.
      include_chain (bool): include the chain id if True

    Warning:
       This will read in the connectivity information from the pdb as well.
       However, a pdb file does not provide any information about the bond
       order. Thus, the bond order of each bond will be set to one.

    Returns:
      (BigDFT.Systems.System): the system file.
    """
    from BigDFT.Fragments import Fragment
    from BigDFT.Systems import System
    from warnings import warn
    from BigDFT.UnitCells import UnitCell

    # First pass read in the atoms.
    sys = System()
    lookup = {}
    sys.conmat = {}
    found = False

    for line in ifile:
        try:  # See if there is an atom on this line
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                at, atid, fragid = _process_atom(line,
                                                 include_chain=include_chain)

                # We can ignore lone pairs
                if at.sym == "Lp":
                    continue

                # Add to the system
                if fragid not in sys:
                    sys[fragid] = Fragment()
                sys[fragid] += [at]

                # Build the lookup table
                lookup[atid] = (fragid, len(sys[fragid]) - 1)

            elif line[:6] == "CONECT":
                found = True
                split = _split_line(line, len(lookup))
                (fragid, atnum) = lookup[int(split[1])]

                for at2 in split[2:]:
                    fragid2, atnum2 = lookup[int(at2)]

                    if fragid not in sys.conmat:
                        sys.conmat[fragid] = []
                        for i in range(0, len(sys[fragid])):
                            sys.conmat[fragid].append({})

                    sys.conmat[fragid][atnum][(fragid2, atnum2)] = 1.0

            elif line[:6] == "CRYST1":
                a = float(line[7:15])
                b = float(line[16:24])
                c = float(line[25:33])
                alpha = float(line[34:40])
                beta = float(line[41:47])
                gamma = float(line[48:54])

                sys.cell = UnitCell([a, b, c], units="angstroem")

                if not disable_warnings:
                    if (alpha != 90 or beta != 90 or gamma != 90):
                        warn("Cell angles must be 90 degrees", UserWarning)

        except IndexError:  # For shorter lines
            continue

    if not found:
        sys.conmat = None
    else:
        # for any connectivity not specified we give default values.
        for fragid in sys:
            if fragid not in sys.conmat:
                sys.conmat[fragid] = []
                for i in range(0, len(sys[fragid])):
                    sys.conmat[fragid].append({})

    if not disable_warnings:
        if sum([len(x) for x in sys.values()]) == 0:
            warn("Warning: zero atoms found", UserWarning)

    return sys


def reorder_fragments(frags):
    chs = {}
    for frag in frags:
        chres, num = frag.split(':')
        if "-" in chres:
            ch, res = chres.split('-')
        else:
            ch = 'A'
        subch = chs.setdefault(ch, {})
        if int(num) in subch:
            return frags  # reordering not possible - incoherent fragments
        else:
            subch[int(num)] = frag
    lis = []
    for ch in sorted(chs.keys()):
        for fr in sorted(chs[ch].keys()):
            lis.append(chs[ch][fr])
    return lis


def write_pdb(system, ofile):
    """
    Write out a system to a string in the pdb format.

    Args:
      system (BigDFT.Systems.System): the system to write.
      ofile (TextIOBase): the output stream to write to.
    """

    from BigDFT.Systems import GetFragTuple

    outstr = ""

    # Write Cell
    if system.cell.get_boundary_condition() != "free":
        line = list(" " * 80)
        line[0:6] = "CRYST1".ljust(6)  # CRYST1
        pos = ["{:.3f}".format(x) for x in
               system.cell.get_posinp(units="angstroem")]
        line[7:15] = pos[0].rjust(8)  # a (Angstrom)
        line[16:24] = pos[1].rjust(8)  # b (Angstrom)
        line[25:33] = pos[2].rjust(8)  # c (Angstrom)
        line[34:40] = "90".rjust(6)  # alpha (degrees)
        line[41:47] = "90".rjust(6)  # beta (degrees)
        line[48:54] = "90".rjust(6)  # gamma (degrees)
        outstr += "".join(line) + "\n"

    # Write Postions
    idx = 1
    lookup = {}
    reorder = reorder_fragments(system)
    for fragid in reorder:  # system.items():
        frag = system[fragid]
        if any([x in fragid for x in _standard_residues]):
            atom_type = "ATOM".ljust(6)
        else:
            atom_type = "HETATM"
        for i, at in enumerate(frag):
            pos = [str("{:.3f}".format(x))
                   for x in at.get_position("angstroem", cell=system.cell)]
            fragtuple = GetFragTuple(fragid)
            resname = fragtuple[0]
            chres = resname.split('-')
            chain = 'A' if '-' not in resname else chres[0]
            resname = chres[-1][:3].ljust(3)
            line = list(" " * 80)
            line[0:6] = atom_type  # HETATM
            line[6:11] = str(idx).rjust(5)  # SERIAL NUMBER
            if "name" in at:  # ATOM NAME
                line[12:16] = at["name"].ljust(4)
            else:
                line[12:16] = at.sym.ljust(4)
            line[16:17] = " "  # ALTERNATIVE LOCATION INDICATOR
            line[17:20] = resname  # RESIDUE NAME
            line[21:22] = chain  # CHAIN IDENTIFIER
            line[22:26] = fragtuple[1].rjust(4)  # RESIDUE SEQUENCE NUMBER
            line[26:27] = " "  # CODE FOR INSERTION OF RESIDUES
            line[30:38] = pos[0].rjust(8)  # X COORDINATE
            line[38:46] = pos[1].rjust(8)  # Y COORDINATE
            line[46:54] = pos[2].rjust(8)  # Z COORDINATE
            line[55:61] = " 1.00 "  # OCCUPANCY
            line[61:67] = " 0.00 "  # TEMPERATURE
            line[72:76] = " B  "  # SEGMENT IDENTIFIER
            line[76:78] = at.sym.rjust(2)  # ELEMENT SYMBOL
            line[78:80] = "  "  # CHARGE
            outstr += "".join(line) + "\n"

            # Keep track of the indexes
            lookup[(fragid, i)] = idx
            idx = idx + 1

    # Write the connectivity information
    if system.conmat is not None:
        for fragid, frag in system.items():
            for i, at in enumerate(frag):
                connections = system.conmat[fragid][i]
                connections = [lookup[x] for x in connections.keys()]

                line = list(" " * 80)
                line[0:6] = "CONECT"

                # SERIAL NUMBER
                line[7:11] = str(lookup[(fragid, i)]).rjust(4)
                if len(connections) > 0:  # BOND SERIAL NUMBERS
                    line[12:16] = str(connections[0]).rjust(4)
                if len(connections) > 1:
                    line[17:21] = str(connections[1]).rjust(4)
                if len(connections) > 2:
                    line[22:26] = str(connections[2]).rjust(4)
                if len(connections) > 3:
                    line[27:31] = str(connections[3]).rjust(4)

                outstr += "".join(line) + "\n"

    ofile.write(outstr)


def read_mol2(ifile, disable_warnings=False):
    """
    Read a system from a mol2 file.

    Args:
      ifile (TextIOBase): the file to read from.
      disable_warnings (bool): whether to print warnings about possible file
        issues.

    Returns:
      (BigDFT.Systems.System): the system file.
    """
    from BigDFT.Systems import System
    from BigDFT.Fragments import Fragment
    from BigDFT.Atoms import Atom
    from BigDFT.UnitCells import UnitCell
    from warnings import warn

    sys = System()

    # Just go ahead and read the whole file into a string.
    lines = [x for x in ifile]

    # First pass, read in the number of atoms.
    for start, line in enumerate(lines):
        if ("@<TRIPOS>MOLECULE" in line):
            break
    start += 1

    split = lines[start+1].split()
    natoms = int(split[0])
    nbonds = int(split[1])

    # Second pass read in the atoms.
    for start, line in enumerate(lines):
        if ("@<TRIPOS>ATOM" in line):
            break
    start += 1

    lookup = []
    for i in range(0, natoms):
        split = lines[start + i].split()
        pos = [float(x) for x in split[2:5]]
        name = split[5]
        sym = name.split(".")[0]
        fragid = split[7] + ":" + split[6]
        charge = [float(split[8])]

        # Add fragment
        if fragid not in sys:
            sys[fragid] = Fragment()
        at = Atom({sym: pos, "units": "angstroem", "q0":
                   charge, "name": name})
        sys[fragid] += [at]

        # Lookup table for connectivity
        lookup.append((fragid, len(sys[fragid]) - 1))

    # Third pass reads the connectivity.
    for start, line in enumerate(lines):
        if ("@<TRIPOS>BOND" in line):
            break
    start += 1

    if start < len(lines):
        sys.conmat = {}
        for fragid, frag in sys.items():
            sys.conmat[fragid] = []
            for i in range(0, len(frag)):
                sys.conmat[fragid].append({})

    bowarn = False
    for i in range(0, nbonds):
        split = lines[start + i].split()
        frag1, at1 = lookup[int(split[1])-1]
        frag2, at2 = lookup[int(split[2])-1]
        bo = split[3]
        try:
            bo = float(split[3])
        except ValueError:
            bowarn = True
            bo = 1
        sys.conmat[frag1][at1][(frag2, at2)] = bo

        # Since mol2 doesn't include the symmetric bonds.
        if frag1 != frag2 or at1 != at2:
            sys.conmat[frag2][at2][(frag1, at1)] = bo

    # Fourth path reads the unit cell.
    for start, line in enumerate(lines):
        if ("@<TRIPOS>CRYSIN" in line):
            break
    start += 1

    if start < len(lines):
        split = lines[start].split()
        a = float(split[0])
        b = float(split[1])
        c = float(split[2])
        alpha = float(split[3])
        beta = float(split[4])
        gamma = float(split[5])

        sys.cell = UnitCell([a, b, c], units="angstroem")

        if not disable_warnings:
            if (alpha != 90 or beta != 90 or gamma != 90):
                warn("Cell angles must be 90 degrees", UserWarning)

    if not disable_warnings:
        if sum([len(x) for x in sys.values()]) == 0:
            warn("Warning: zero atoms found", UserWarning)
        if bowarn:
            warn("Unsupported bond type had to be set to 1 (i.e. aromatic)",
                 UserWarning)

    return sys


def write_mol2(system, ofile):
    """
    Write out a system to a string in the mol2 format.

    Args:
      system (BigDFT.Systems.System): the system to write.
      ofile (TextIOBase): the output stream to write to.
    Returns:
      (str): a string representation of the file.
    """
    from BigDFT.Systems import GetFragTuple

    charges = True

    idx = 1
    lookup = {}
    atomstr = "@<TRIPOS>ATOM\n"
    for fragid, frag in system.items():
        for i, at in enumerate(frag):
            pos = [str("{:.3f}".format(x))
                   for x in at.get_position("angstroem", cell=system.cell)]
            fragtuple = GetFragTuple(fragid)

            line = str(idx) + " " + at.sym+str(idx) + " "
            line += pos[0] + " " + pos[1] + " " + pos[2] + " "
            line += at.sym + " "
            line += fragtuple[1] + " "
            line += fragtuple[0] + " "

            if at.q0:
                line += str(at.q0)
            else:
                line += str(0.0)

            atomstr += line + "\n"

            # Keep track of the indexes
            lookup[(fragid, i)] = idx
            idx = idx + 1
    numatoms = idx - 1

    # Write the connectivity information
    bondstr = "@<TRIPOS>BOND\n"
    numbonds = 0
    idx = 1
    if system.conmat is not None:
        for fragid, frag in system.items():
            for i, at in enumerate(frag):
                connections = system.conmat[fragid][i]
                bonds = connections.values()
                connections = [lookup[x] for x in connections.keys()]

                for con, bo in zip(connections, bonds):
                    if con <= idx:
                        continue
                    bondstr += str(numbonds + 1) + " "
                    bondstr += str(idx) + " " + str(con) + " "
                    bondstr += str(int(bo)) + "\n"
                    numbonds = numbonds + 1
                idx = idx + 1

    # Write out the molecule part
    molstr = "@<TRIPOS>MOLECULE\n"
    molstr += "generated by BigDFT\n"
    molstr += str(numatoms) + " " + str(numbonds) + " "
    molstr += str(len(system.keys())) + " 0 0\n"
    molstr += "****\n"
    if charges:
        molstr += "MULLIKEN_CHARGES\n"
    else:
        molstr += "NO_CHARGES\n"
    molstr += "****\n"
    molstr += "****\n"

    # Write Cell
    cellstr = ""
    if system.cell.get_boundary_condition() != "free":
        cellstr += "@<TRIPOS>CRYSIN\n"
        pos = ["{:.3f}".format(x) for x in
               system.cell.get_posinp(units="angstroem")]
        cellstr += " ".join(pos)
        cellstr += " 90 90 90"
        cellstr += " 1 1\n"

    ofile.write(molstr + "\n" + atomstr + "\n" + bondstr + "\n" + cellstr)


def write_xyz(system, ofile):
    """
    Write out a system to a file in the xyz format.

    Args:
      system (BigDFT.Systems.System): the system to write.
      ofile (TextIOBase): the output stream to write to.
      cell (list): the unit cell.
    """
    outstr = ""
    outstr += str(sum([len(x) for x in system.values()])) + " angstroem\n"
    outstr += system.cell.get_boundary_condition("angstroem") + "\n"

    for frag in system.values():
        for at in frag:
            pos = at.get_position("angstroem", cell=system.cell)
            outstr += at.sym + " " + " ".join([str(x) for x in pos]) + "\n"

    ofile.write(outstr)


def _process_atom(line, include_chain=False):
    """
    This processes a line of a pdb file and extracts information about an atom.

    Returns:
      (BigDFT.Atoms.Atom, int, str): return the Atom on this line, its id in
      the pdb, and the fragment it belongs to. It may include the chain id.
    """
    from BigDFT.Atoms import Atom
    # Get the basic information about this atom
    sym = line[76:78].strip().capitalize()

    if "1-" in sym:
        sym = sym.replace("1-", "")
    if "1+" in sym:
        sym = sym.replace("1+", "")

    fullname = line[12:16]

    name = fullname.strip()

    xpos = float(line[30:38])
    ypos = float(line[38:46])
    zpos = float(line[46:54])

    at = Atom({"sym": sym, "r": [xpos, ypos, zpos], "name": name,
               "units": "angstroem"})

    # Get the atom id for building the lookup table
    atid = int(line[6:11])

    # Information about the residue
    resname = line[17:20]
    chain = line[20:22].strip()
    resid = str(int(line[22:26]))
    fragid = resname+":"+resid
    if include_chain:
        fragid = chain+'-'+fragid

    return at, atid, fragid


class XYZReader():
    """
    A class which can be used to read from xyz files.

    This class should behave like a standard file, which means you can use
    it in ``with`` statements, and use the ``next`` command to iterate over it.

    Attributes:
      closed (bool): check if a file is open or not.
      units (str): the units that the xyz file was in.
      natoms (int): the number of atoms in the xyz file.
      cell (list): a list of floats describing the cell.

    Args:
      filename (str): the file to read from. You can also specify a molecule
        that might be in the database.
    """

    def __init__(self, filename):
        self.filename = filename
        self.natoms = None
        self.cell = None
        self.atoms_positions = []

    def open(self):
        self.__enter__()

    def __enter__(self):
        from os.path import join, abspath, dirname
        from BigDFT.UnitCells import UnitCell

        try:
            self._handle = open(self.filename, "r")
        except IOError as e:  # see if the file is in the database
            dirXYZ = join(dirname(__file__), "Database", "XYZs")
            filename = abspath(join(dirXYZ, self.filename + ".xyz"))
            try:
                self._handle = open(filename, "r")
            except IOError:  # Raise the error from the original path
                raise e

        # Read atoms and units.
        line = next(self._handle)
        self.natoms = int(line.split()[0])
        try:
            self.units = line.split()[1]
        except IndexError:
            self.units = "angstroem"

        # Read boundary condition.
        line = next(self._handle).split()
        if len(line) == 0:
            self.cell = UnitCell()
        elif line[0] == "wire":
            self.cell = UnitCell([float("inf"), float("inf"),
                                  float(line[3])], units=self.units)
        elif line[0] == "surface":
            self.cell = UnitCell([float(line[1]), float("inf"),
                                  float(line[3])], units=self.units)
        elif line[0] == "periodic":
            self.cell = UnitCell([float(x) for x in line[1:4]],
                                 units=self.units)
        else:
            self.cell = UnitCell()

        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        from BigDFT.Atoms import Atom
        line = next(self._handle)
        if self.natoms == len(self.atoms_positions):
            raise StopIteration
        split = line.split()
        sym = split[0]
        position = [float(x) for x in split[1:4]]
        this_pos = Atom({'sym': sym, 'r': position, "units": self.units})
        self.atoms_positions.append(this_pos)
        return this_pos

    def __iter__(self):
        return self

    def close(self):
        if self.natoms != len(self.atoms_positions):
            raise IOError('The number of atoms is not consistent with the'
                          ' number of lines')
        self.__exit__()

    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        self._handle.close()

    @property
    def closed(self):
        return self._handle.closed


class XYZWriter():
    """
    A class for writing XYZ files.

    This class should behave like a standard file, which means you can use
    it in ``with`` statements and write.

    Args:
      filename (str): the file to write to.
      natoms (int): how many atoms we will write.
      units (str): the units of the file. Defaults to angstroem.
      cell (BigDFT.UnitCells.UnitCell): the unit cell.
    """

    def __init__(self, filename, natoms, units="angstroem", cell=None):
        from BigDFT.UnitCells import UnitCell
        self.filename = filename
        self.natoms = natoms
        self.units = units
        if cell is None:
            self.cell = UnitCell()
        else:
            self.cell = cell

    def open(self):
        self.__enter__()

    def __enter__(self):
        self._handle = open(self.filename, "w")
        self._handle.write(str(self.natoms) + " ")
        self._handle.write(self.units)
        self._handle.write("\n")

        # The unit cell
        self._handle.write(self.cell.get_boundary_condition(self.units))
        self._handle.write("\n")

        return self

    def write(self, atomdict):
        """
        Write an atom to the file.

        Args:
          atom (dict): a dictionary describing the atom.
        """
        from BigDFT.Atoms import Atom
        at = Atom(atomdict)
        self._handle.write(at.sym + " ")
        pos = at.get_position(self.units)
        self._handle.write(" ".join([str(x) for x in pos]))
        self._handle.write("\n")

    def close(self):
        self.__exit__()

    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        self._handle.close()

    @property
    def closed(self):
        return self._handle.closed


def _example():
    """Test the XYZ Module"""
    from BigDFT.Systems import System
    from BigDFT.Fragments import Fragment
    from BigDFT.UnitCells import UnitCell
    file = "Si4"

    safe_print("First let's try reading an XYZ file.")
    atom_list = []
    with XYZReader(file) as reader:
        safe_print(reader.closed)
        for at in reader:
            atom_list.append(at)
    safe_print(reader.closed)
    safe_print(atom_list)
    safe_print()

    safe_print("Now let's try writing an XYZ file.")
    safe_print()
    with XYZWriter("test.xyz", len(atom_list),
                   units=reader.units) as writer:
        safe_print(writer.closed)
        for at in atom_list:
            writer.write(at)
    safe_print(writer.closed)
    safe_print()
    with open("test.xyz") as ifile:
        for line in ifile:
            safe_print(line, end='')
    safe_print()

    safe_print("Print with various boundary conditions")
    with XYZWriter("test.xyz", len(atom_list), reader.units,
                   cell=UnitCell()) as writer:
        for at in atom_list:
            writer.write(at)
    with XYZReader("test.xyz") as ifile:
        print(ifile.cell.get_boundary_condition())

    with XYZWriter("test.xyz", len(atom_list), reader.units,
                   cell=UnitCell([5, 5, 5])) as writer:
        for at in atom_list:
            writer.write(at)
    with XYZReader("test.xyz") as ifile:
        print(ifile.cell.get_boundary_condition())

    with XYZWriter("test.xyz", len(atom_list), reader.units,
                   cell=UnitCell([5, float("inf"), 5])) as writer:
        for at in atom_list:
            writer.write(at)
    with XYZReader("test.xyz") as ifile:
        print(ifile.cell.get_boundary_condition())

    with XYZWriter("test.xyz", len(atom_list), reader.units,
                   cell=UnitCell([float("inf"), float("inf"), 5])) as writer:
        for at in atom_list:
            writer.write(at)
    with XYZReader("test.xyz") as ifile:
        print(ifile.cell.get_boundary_condition())
    safe_print()

    safe_print("Now let's demonstrate the pdb and mol2 writer")
    sys = System()
    sys["FRAG:0"] = Fragment(atom_list)
    with open("test.pdb", "w") as ofile:
        write_pdb(sys, ofile)
    with open("test.pdb") as ifile:
        for line in ifile:
            safe_print(line, end='')
    safe_print()

    with open("test.mol2", "w") as ofile:
        write_mol2(sys, ofile)
    with open("test.mol2") as ifile:
        for line in ifile:
            safe_print(line, end='')
    safe_print()


if __name__ == "__main__":
    _example()
