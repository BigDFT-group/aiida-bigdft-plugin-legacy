"""
This module defines some tools used for reading/writing XYZ files.

"""
from futile.Utils import write as safe_print


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
        try:
            self._handle = open(self.filename, "r")
        except IOError:  # see if the file is in the database
            dirXYZ = join(dirname(__file__), "Database", "XYZs")
            filename = abspath(join(dirXYZ, self.filename + ".xyz"))
            self._handle = open(filename, "r")
        line = next(self._handle)
        self.natoms = int(line.split()[0])
        try:
            self.units = line.split()[1]
        except IndexError:
            self.units = "angstroem"
        line = next(self._handle).split()
        if len(line) == 0:
            self.cell = None
        elif line[0] == "#":  # Comment line
            self.cell = None
        elif line[0] == "free":
            self.cell = None
        else:
            self.cell = [float(x) for x in line[1:]]

        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        line = next(self._handle)
        if self.natoms == len(self.atoms_positions):
            raise StopIteration
        split = line.split()
        sym = split[0]
        position = split[1:4]
        this_pos = {'sym': sym, 'r': position, "units": self.units}
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
      cell (list): the unit cell.
    """

    def __init__(self, filename, natoms, units="angstroem", cell=None):
        self.filename = filename
        self.natoms = natoms
        self.units = units
        self.cell = cell

    def open(self):
        self.__enter__()

    def __enter__(self):
        self._handle = open(self.filename, "w")
        self._handle.write(str(self.natoms) + " ")
        self._handle.write(self.units)
        self._handle.write("\n")

        # The unit cell
        self._handle.write(_xyz_bc_spec(self.cell))
        self._handle.write("\n")

        return self

    def write(self, atomdict):
        """
        Write an atom to the file.

        Args:
          atom (dict): a dictionary describing the atom.
        """
        from BigDFT.Atom import Atom
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


def _xyz_bc_spec(cell):
    """
    Defines the specification for expressing the Boundary Conditions starting
    from a cell vector.

    Args:
      cell (list): array of the (orthorhombic) cell. Should be 0.0 on
        directions with free BC. If None is given, the BC are assumed to
        be Free.

    Return:
     (str): Line of the xyz file specifying the bc
    """
    if cell is None:
        return "free"
    elif cell[1] == 0.0 and cell[2] != 0.0:
        return "surface " + str(cell[0]) + " 0.0 " + str(cell[2]) + " "
    elif cell[1] == 0.0 and cell[2] == 0.0:
        return "wire 0.0 0.0 " + cell[2] + " "
    else:
        return "periodic " + str(cell[0]) + " " + str(cell[1]) + \
               " " + str(cell[2]) + " "


if __name__ == "__main__":
    """Test the XYZ Module"""
    from os.path import join
    from os import system
    for file in ["Si4.xyz", "SiO.xyz"]:
        safe_print("First let's try reading an XYZ file.")
        atom_list = []
        with XYZReader(join("Database", "XYZs", file)) as reader:
            safe_print(reader.closed)
            for at in reader:
                atom_list.append(at)
        safe_print(reader.closed)
        safe_print(atom_list)
        safe_print()

        safe_print("Now let's try writing an XYZ file.")
        safe_print()
        with XYZWriter("test.xyz", len(atom_list), cell=reader.cell,
                       units=reader.units) as writer:
            safe_print(writer.closed)
            for at in atom_list:
                writer.write(at)
        safe_print(writer.closed)
        system("cat test.xyz")
        safe_print()

    safe_print("Using the explicit open and close.")
    reader = XYZReader("test.xyz")
    reader.open()
    print(next(reader))
    reader.close()
