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
      filename (str): the file to read from.
    """

    def __init__(self, filename):
        self.filename = filename
        self.natoms = None
        self.cell = None

    def open(self):
        self.__enter__()

    def __enter__(self):
        self._handle = open(self.filename, "r")
        line = next(self._handle)
        self.natoms = int(line.split()[0])
        self.units = line.split()[1]
        line = next(self._handle).split()
        if len(line) == 0:
            self.cell = None
        elif line[0] == "#": # Comment line
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
        split = line.split()
        sym = split[0]
        position = split[1:]

        return {'sym': sym, 'r': position, "units": self.units}

    def __iter__(self):
        return self

    def close(self):
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
        from Atom import Atom
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
      cell (list): array of the (orthorhombic) cell. Should be 0.0 on directions
        with free BC. If None is given, the BC are assumed to be Free.

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

    # for atom in reader:
    #     safe_print(atom)

# class XYZfile():
#     """
#     .. |filename_docs| replace::
#          The file which will be created. If None, the file will be
#          eventually dumped in :class:~`sys.stdout`.
#
#     A class associated to a xyz input file as processed by BigDFT
#
#     :param filename: |filename_docs|
#     :type filename: string
#     :param units: The units of measure of the positions. Allowed avlues are
#       'atomic' or 'angstroem'
#     :type units: string
#     """
#
#     def __init__(self, filename=None, units='atomic'):
#         self.filename = filename
#         self.lines = []
#         self.units = units
#         self.fac = 1.0
#         if units == 'angstroem':
#             self.fac = AU_to_A
#
#     def append(self, array, basename='', names=None, attributes=None):
#         """
#         Add lines to the file position list
#
#         :param array: list of the atomic positions
#         :type array: list of  float triples
#         :param basename: base for the name of the atoms
#         :type basename: string
#         :param names: list of atom names. Will be appended to `basename` if the latter is present
#         :type names: list of strings
#         :param attributes: list of further attributes to be associated to each of the atoms.
#             Will be serialized close to each of the atomic positions
#         :type attributes: list of dictionaries
#         """
#         nm = basename
#         for i, r in enumerate(array):
#             if names is not None:
#                 nm = basename + names[i]
#             line = str(nm)
#             for t in r:
#                 line += ' ' + str(self.fac * t)
#             if attributes is not None:
#                 line += ' ' + str(attributes[i])
#             self.lines.append(line + '\n')
#
#     def dump(self, position='w'):
#         """
#         Dump the file on the file system if filename has been provided,
#         otherwise dump on sys.stdout.
#
#         :param position: filename position statement. Only menaingful for a file dumping.
#         :type position: char
#         """
#         import sys
#         f = sys.stdout
#         if self.filename is not None:
#             f = open(self.filename, position)
#         f.write(str(len(self.lines)) + ' ' + str(self.units) + '\n')
#         f.write('# xyz dump \n')
#         # then the positions
#         for l in self.lines:
#             f.write(l)
#         if self.filename is not None:
#             f.close()
#
#
# def open_xyz(filename, nat, unit, comment, position='a'):
#     import sys
#     f = sys.stdout
#     if filename is not None:
#         f = open(filename, position)
#     if (position != 'a'):
#         f.write(str(nat) + ' ' + str(unit) + '\n')
#         f.write(comment + '\n')
#     return f
#
#
# def close_xyz(f, filename):
#     if filename is not None:
#         f.close()
#
#
# def dump_xyz_positions(f, array, basename='', names=None):
#     nm = basename
#     for i, r in enumerate(array):
#         if names is not None:
#             nm = basename + names[i]
#         f.write(str(nm) + ' ' + str(r[0]) + ' ' +
#                 str(r[1]) + ' ' + str(r[2]) + '\n')
#
#
#
#
# def dump_xyz(array, basename='', units='atomic', names=None, filename=None, position='a', comment=None, cell=None):
#     """
#     Create a BigDFT xyz filename. Duplicates the meaning of the :class:`XYZfile` class.
#
#     :param filename: |filename_docs|
#     :type filename: string
#
#     .. todo::
#        Remove the duplication of operations in favour or the class.
#        Move the related operation into a lower-level module.
#     """
#     cmt = xyz_bc_spec(cell)
#     cmt += comment if comment is not None else '# xyz dump with basename "' + basename + '"'
#     f = open_xyz(filename, len(array), units, cmt, position)
#     dump_xyz_positions(f, array, basename=basename, names=names)
#     close_xyz(f, filename)
