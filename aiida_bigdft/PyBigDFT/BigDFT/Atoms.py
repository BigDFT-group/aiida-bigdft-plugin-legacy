'''
This module defines the atom class, which is a class which contains very
general descriptions of a single atom.
'''
try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping
from futile.Utils import write as safe_print

#: Conversion between Atomic Units and Bohr
AU_to_A = 0.52917721092
#: A list of valid keys for describing a multipole.
MULTIPOLE_ANALYSIS_KEYS = ['q0', 'q1', 'q2', 'sigma', 'multipole character']


def atomic_number(sym):
    """
    Returns the atomic number associated with the given symbol.

    Args:
      sym (str): the atomic symbol of this atom.

    Returns:
      (int): the atomic number
    """
    return _atomic_number[sym]


def atomic_weight(sym):
    """
    Returns the weight number associated with the given symbol.

    Args:
      sym (str): the atomic weight of this atom.

    Returns:
      (int): the atomic weight
    """
    return _atomic_weight[sym]


def number_to_symbol(number):
    """
    Returns the symbol of atoms with a given atomic number.

    Args:
      number (int): the atomic number to lookup.

    Returns:
      (str): the atomic symbol with the given number.
    """
    return [x for x in _atomic_number if _atomic_number[x] == number][0]


_atomic_number = {"H": 1, "D": 1, "He": 2,
                  "Li": 3, "Be": 4, "B": 5, "C": 6,
                  "N": 7, "O": 8, "F": 9, "Ne": 10,
                  "Na": 11, "Mg": 12, "Al": 13, "Si": 14,
                  "P": 15, "S": 16, "Cl": 17, "Ar": 18,
                  "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
                  "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
                  "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
                  "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
                  "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38,
                  "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
                  "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46,
                  "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
                  "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
                  "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
                  "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62,
                  "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 56,
                  "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
                  "Lu": 71, "Hf": 72, "W": 74, "Re": 75,
                  "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
                  "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
                  "Po": 84, "At": 85, "Rn": 86}

_atomic_weight = {"H": 1.00794, "D": 2.0141, "He": 4.003,
                  "Li": 6.941, "Be": 9.012182, "B": 10.811,
                  "C": 12.0107, "N": 14.00674, "O": 15.9994,
                  "F": 18.9984032, "Ne": 20.1797, "Na": 22.989770,
                  "Mg": 24.3050, "Al": 26.981538, "Si": 28.0855,
                  "P": 30.973761, "S": 32.066, "Cl": 35.4527,
                  "Ar": 39.948, "K": 39.0983, "Ca": 40.078,
                  "Sc": 44.955910, "Ti": 47.867, "V": 50.9415,
                  "Cr": 51.9961, "Mn": 54.938049, "Fe": 55.845,
                  "Co": 58.933200, "Ni": 58.6934, "Cu": 63.546,
                  "Zn": 65.39, "Ga": 69.723, "Ge": 72.61,
                  "As": 74.92160, "Se": 78.96, "Br": 79.904, "Kr": 83.80,
                  "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585, "Zr": 91.224,
                  "Nb": 92.90638, "Mo": 95.94, "Tc": 98, "Ru": 101.07,
                  "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411,
                  "In": 114.818, "Sn": 118.710, "Sb": 121.760, "Te": 127.60,
                  "I": 126.90447, "Xe": 131.29, "Cs": 132.90545196,
                  "Ba": 37.327, "La": 138.90547, "Ce": 140.116,
                  "Pr": 140.90766, "Nd": 144.242, "Pm": None, "Sm": 150.36,
                  "Eu": 151.964, "Gd": 157.25, "Tb": 158.925354, "Dy": 162.500,
                  "Ho": 164.930328, "Er": 167.259, "Tm": 168.934218,
                  "Yb": 173.045, "Lu": 174.9668, "Hf": 178.486,
                  "Ta": 180.94788, "W": 183.84, "Re": 186.207, "Os": 190.23,
                  "Ir": 192.217, "Pt": 195.084, "Au": 196.966570,
                  "Hg": 200.592, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98040,
                  "Po": None, "At": None, "Rn": None}

_nzion_default_psp = {"H":  1.0, "D":  1.0, "He": 2.0,
                      "Li": 1.0, "Be": 2.0, "B":  3.0, "C":  4.0,
                      "N": 5.0, "O": 6.0, "F":  7.0, "Ne": 8.0,
                      "Na": 1.0, "Mg": 2.0, "Al": 3.0, "Si": 4.0,
                      "P":  5.0, "S":  6.0, "Cl": 7.0, "Ar": 8.0,
                      "Cu": 11.0, "Zn": 12.0, "Br": 7.0}


class Atom(MutableMapping):
    """
    Defines a wrapper for atoms.

    An atom may have many quantities associated with it. These quantities
    are get and set in a dictionary like fashion, allowing an atom to
    dynamically hold whatever data you need. However, we still wrap it in a
    class so that we can have some common operations for it, as well as so we
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

    @property
    def atomic_number(self):
        """
        The atomic number of this atom.
        """
        return atomic_number(self.sym)

    @property
    def atomic_weight(self):
        """
        The atomic number of this atom.
        """
        return atomic_weight(self.sym)

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
            if k in self:
                val = self[k]
                if isinstance(val, ndarray):
                    return_dict[k] = list(val)
                else:
                    return_dict[k] = val

        return return_dict

    def serialize(self, units='bohr'):
        """
        Transform the atom in a dictionary that can be employed for
        the construction of dataframes or pandas series.

        Args:
            units (str): the units for the positions
        Returns:
            dict: the serialized dictionary
        """
        xyz = ['x', 'y', 'z']
        atdict = {}
        for key, val in self.dict().items():
            if key == self.sym:
                atdict['sym'] = self.sym
                reval = self.get_position(units=units)
                atdict.update({t+'_coord': reval[i]
                               for i, t in enumerate(xyz)})
                q_ion = self.get('nzion')  # nzion(self.sym)
                if q_ion is not None:
                    atdict['zion'] = q_ion
                mchar = self.dict().get('multipole character')
                if mchar is not None:
                    q0 = self.dict()['q0'][0]
                    if q_ion is not None:
                        if mchar == 'gross':
                            atdict['qel_0'] = q0
                        else:
                            atdict['qel_0'] = q0 - q_ion
                    else:
                        if mchar == 'net':
                            atdict['qel_0'] = q0
            elif key == 'r':
                fac = AU_to_A if self.dict()['units'] == 'angstroem' else 1.0
                atdict.update({key+'_'+str(i): t/fac
                               for i, t in enumerate(val)})
            elif key == 'units':
                atdict['units'] = units
            elif isinstance(val, list):
                atdict.update({key+'_'+str(i): t
                               for i, t in enumerate(val)})
            else:
                atdict[key] = val
            try:
                atdict["nel"] = self.nel
            except Exception:
                atdict["nel"] = None
        return atdict

    @property
    def is_link(self):
        """
        Whether or not this atom is a link atom or not.
        """
        try:
            return self["link_atom"]
        except KeyError:
            return False

    @is_link.setter
    def is_link(self, v):
        self["link_atom"] = v

    @property
    def nel(self):
        """
        The number of electrons in this atom.
        """
        if "nzion" in self:
            return self["nzion"]
        elif self.sym in _nzion_default_psp:
            return _nzion_default_psp[self.sym]
        else:
            raise Exception("Number of electrons not set for this atom",
                            "either explicitly set the nzion key or ",
                            "try something like set_electrons_from_log")

    @nel.setter
    def nel(self, val):
        self["nzion"] = val

    @property
    def q0(self):
        """
        Provides the charge of the atom.
        """
        charge = self.get('q0')
        if charge is not None:
            charge = charge[0]
        return charge

    @property
    def q1(self):
        """
        Provides the dipole of the atom.
        """
        import numpy as np
        dipole = self.get('q1')  # they are (so far) always given in AU
        if dipole is not None:
            dipole = np.array([dipole[2], dipole[0], dipole[1]])
        return dipole

    def set_multipole(self, mp, correct_charge=True):
        """
        Given another atom or a dictionary, this sets the multipole related
        values of this with those values.

        Todo:
          Arrive at a standard that avoids having to do the charge
          correction here.

        Args:
          mp (dict): a dictionary which contains information about multipoles.
          correct_charge (bool): currently there is an inconsistency in
            terms of gross charge, and this corrects it.
        """
        from copy import deepcopy

        for key in MULTIPOLE_ANALYSIS_KEYS:
            if key in mp:
                self[key] = deepcopy(mp[key])

        # Correct the charge
        if correct_charge and "q0" in self and "nzion" in mp:
            self["q0"][0] += mp["nzion"]
            self['multipole character'] = 'net'

    def get_force(self):
        """
        Returns the force on the atom in the desired units.

        Returns:
          An array of position values.
        """
        return self.get("force")

    def set_force(self, force):
        """
        Given an atom or a dictionary, this sets the force.

        Args:
          force (list): a list of force values.
        """
        self.store["force"] = force

    @property
    def sym(self):
        sym = _GetSymbol(self.store)
        if sym == 'r':
            sym = self.store['sym']
        return sym

    @sym.setter
    def sym(self, v):
        if 'sym' in self.store:
            self.store['sym'] = v
        else:
            sym = self.sym
            val = self.store[sym]
            del self[sym]
            self[v] = val

    def get_position(self, units="bohr", cell=None):
        """
        Returns the position of the atom in the desired units.

        Args:
          units (str): the units to return the position in. Default is bohr.
          cell (list): the vectors of the unit cell (3x3). If passed, the
            minimum image convention is enforced.

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
        if IsAngstroem(self):
            pos /= AU_to_A
        if IsAngstroem(units):
            pos *= AU_to_A

        # Enforce minimum image convention
        if cell is not None:
            pos = cell.minimum_image(pos, units)

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
        if IsAngstroem(units):
            pos /= AU_to_A
        if IsAngstroem(self):
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

        # Upcast to an Atom
        othercomp = Atom(other)

        # Compare Symbols
        sym1 = self.sym.title()
        sym2 = othercomp.sym.title()
        if sym1 != sym2:
            return False

        # Compare position
        pos1 = array(self.get_position())
        pos2 = array(othercomp.get_position())

        return norm(pos1 - pos2) < 1e-3


def _GetSymbol(atom):
    """
    Provides the key which defines the element of the of atom.

    Arguments:
      atom (dict): a dictionary describing the atom.
    Returns:
      (str): the symbol the atom.
    Raises:
       ValueError: atom with no symbol
    """
    ks = atom.keys()
    if 'sym' in ks:
        sym = atom['sym']
        if len(sym) == 0:  # fallback in case of no symbol
            sym = atom['name'][:1]
        return sym

    for k in ks:
        if k.title() in _atomic_number and isinstance(atom[k], list):
            if len(atom[k]) == 3:
                return k

    raise ValueError


def IsAngstroem(units):
    """
    Checks if a string or atom has angstroem as its units.

    Args:
      units: either a string, or a (BigDFT.Atoms.Atom).
    """
    if hasattr(units, 'store'):
        check = units.store.get("units")
        if not check:
            return False
    else:
        check = units
    return check == "angstroem" or check == "angstroemd0"


def _example():
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


if __name__ == "__main__":
    _example()
