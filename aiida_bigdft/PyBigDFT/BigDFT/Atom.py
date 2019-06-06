'''
This module defines the atom class, which is a class which contains very
general descriptions of a singel atom.
'''
try:
    from collections.abc import MutableMapping
except:
    from collections import MutableMapping
from futile.Utils import write as safe_print

AU_to_A = 0.52917721092
MULTIPOLE_ANALYSIS_KEYS = ['q0', 'q1', 'q2', 'sigma', 'multipole character']
PROTECTED_KEYS = MULTIPOLE_ANALYSIS_KEYS + ["frag"] + ["r", "units"]


def nzion(sym):
    """
    Returns the charge of the nucleus associated with an atomic symbol.

    Args:
      sym (str): the atomic symbol of this atom.

    Returns:
      (float): the charge of that nucleus.
    """
    return _nzion[sym]


_nzion = {"H":  1.0, "He": 2.0, \
          "Li": 1.0, "Be": 2.0, "B":  3.0, "C":  4.0,
          "N":  5.0, "O":  6.0, "F":  7.0, "Ne": 8.0}


class Atom(MutableMapping):
    """
    Defines a wrapper for atoms.

    An atom may have many quantities associated with it. These quantities
    are get and set in a dictionary like fashion, allow this atom to
    dynamically hold whatever data you need. But we still wrap it in a class
    so that we can have some common operations for it, as well as so we
    can maintain suitable units.

    It is this class's responsibility to extract the main properties of an
    atom (position, symbol) from the dictionary.

    Args:
      data (dict):
        A dictionary of miscellaneous values to associate with this atom.
    """

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))

    def dict(self):
        """
        Convert to a dictionary.
        """
        return self.store

    def get_external_potential(self, units="bohr"):
        """
        Transform the atom into a dictionary ready to be put as external
        potential.
        """
        from numpy import ndarray
        return_dict = {}
        return_dict["sym"] = self.sym
        return_dict["r"] = self.get_position(units)
        for k in MULTIPOLE_ANALYSIS_KEYS:
            val = self[k]
            if isinstance(val, ndarray):
                return_dict[k] = list(val)
            else:
                return_dict[k] = val

        return return_dict

    @property
    def q0(self):
        """
        Provides the charge of the atom
        """
        charge = self.get('q0')
        if charge is not None:
            charge = charge[0]
        return charge

    @property
    def q1(self):
        """
        Provides the dipole of the atom
        """
        import numpy as np
        dipole = self.get('q1')  # they are (so far) always given in AU
        if dipole is not None:
            dipole = np.array([dipole[2], dipole[0], dipole[1]])
        return dipole

    def set_multipole(self, mp):
        """
        Given another atom or a dictionary, this sets the multipole related
        values of this with those values.

        Args:
          mp (dict/Atom): an atom which contains information about multipoles.
        """
        for key in MULTIPOLE_ANALYSIS_KEYS:
            self[key] = mp[key]

    @property
    def sym(self):
        sym = _GetSymbol(self.store)
        if sym == 'r':
            sym = self.store['sym']
        return sym

    def get_position(self, units="bohr"):
        """
        Returns the position of the atom in the desired units.

        Args:
          units (str): the units to return the position in. Default is bohr.

        Returns:
          An array of position values.
        """
        from numpy import array

        # Grab the position from the store
        if 'r' in self.store:
            pos = self.store['r']
        else:
            pos = self.store[self.sym]
        pos = array([float(x) for x in pos])

        # Make sure the units are correct
        if _IsAngstroem(self):
            pos /= AU_to_A
        if _IsAngstroem(units):
            pos *= AU_to_A

        return [float(x) for x in pos]

    def set_position(self, new_pos, units="bohr"):
        """
        Set the position of the atom.

        Args:
          new_pos (list): a list of floats defining the new position.
          units(str): the units of the new position being passed. Default is
            bohr.
        """
        from numpy import array
        # Convert the input to the right units.
        pos = array(new_pos)
        if _IsAngstroem(units):
            pos /= AU_to_A
        if _IsAngstroem(self):
            pos *= AU_to_A
        pos = [x for x in pos]

        # Insert
        if 'r' in self.store:
            self.store['r'] = pos
        else:
            self.store[self.sym] = pos
        pass

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        if key == self.sym:
            raise ValueError("You can't delete the symbol from an atom.")
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

    def __eq__(self, other):
        """
        Compare two atoms. They are equal if they have the same position and
        symbol.

        other (dict, Atom): the atom (or something that can be cast to one)
          to compare with.
        """
        from numpy.linalg import norm
        from numpy import array

        sym1 = self.sym
        pos1 = array(self.get_position())

        othercomp = Atom(other)
        pos2 = array(othercomp.get_position())
        sym2 = othercomp.sym

        return norm(pos1 - pos2) < 1e-3 and sym1 == sym2


def _GetSymbol(atom):
    """
    Provides the key which defines the element of the of atom.

    Arguments:
      atom (dict): a dictionary describing the atom.
    Returns:
      (str): the symbol the atom.
    """
    ks = atom.keys()
    if 'sym' in ks:
        return atom['sym']

    for k in ks:
        if k not in PROTECTED_KEYS and type(atom[k]) == type([]):
            if len(atom[k]) == 3:
                return k

    raise ValueError


def _IsAngstroem(units):
    """
    Checks if a string or atom has angstroem as its units.

    Args:
      units: either a string, or an ``Atom``.
    """
    if isinstance(units, Atom):
        check = units.store.get("units")
        if not check:
            return False
    else:
        check = units
    return check == "angstroem" or check == "angstroemd0"


if __name__ == "__main__":
    """Test the atom module"""
    safe_print("Access the full data")
    test_atom = Atom({'r': [1.0, 0.0, 0.0], 'sym': "He", 'units': 'bohr'})
    safe_print(dict(test_atom))
    # Access the derived data
    safe_print(test_atom.sym)
    safe_print(test_atom.get_position())
    safe_print(test_atom.get_position('angstroem'))
    safe_print()

    safe_print("Create a new atom with different units")
    new_atom = Atom({
        'r': [float(x) for x in test_atom.get_position('angstroem')],
        'sym': test_atom.sym, 'units': 'angstroem'})
    safe_print("Are these atoms equal?")
    safe_print(new_atom == test_atom)
    safe_print()

    safe_print("Now other times we get an array that looks more like this")
    test_atom = Atom(He=[1.0, 0.0, 0.0], units='bohr')
    safe_print(dict(test_atom))
    safe_print("But everything else works as expected")
    safe_print(test_atom.sym)
    safe_print(test_atom.get_position())
    safe_print(new_atom == test_atom)
    safe_print()

    safe_print("The atom can be used as a dict for adding new properties.")
    test_atom["frag"] = "ANA"
    for key, value in test_atom.items():
        safe_print(key, value)
    safe_print()
    safe_print("And if we update the dictionary position or symbol,")
    safe_print("everything else reacts with suitable caution.")
    test_atom["He"] = [-1.0, 0.0, 0.0]
    safe_print(dict(test_atom))
    safe_print(test_atom.get_position('angstroem'))

    safe_print("There is a protection to prevent you from deleting the")
    safe_print("symbol of an atom from the dictionary.")
    try:
        del test_atom["frag"]
        del test_atom["He"]
    except ValueError as v:
        safe_print(v)
    safe_print(dict(test_atom))
    safe_print()

    safe_print("But you can change the symbol if you are working with the")
    safe_print("other representation.")
    safe_print(dict(new_atom))
    new_atom["sym"] = "Na"
    safe_print(new_atom.sym)
    safe_print(dict(new_atom))
    safe_print()

    safe_print("One final check of the atom comparison")
    new_atom["units"] = "bohr"
    new_atom["r"] = [-1.0, 0.0, 0.0]
    safe_print(new_atom.sym, new_atom.get_position())
    safe_print(test_atom.sym, test_atom.get_position())
    safe_print(new_atom == test_atom)
    safe_print()

    safe_print("We can also update the position")
    safe_print(test_atom.get_position())
    safe_print(dict(test_atom))
    test_atom.set_position([1.0, 1.0, 1.0], units="angstroem")
    safe_print(test_atom.get_position(units="angstroem"))
    safe_print(dict(test_atom))
    safe_print()
