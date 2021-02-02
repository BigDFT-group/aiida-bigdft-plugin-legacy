"""
This module contains data structures for describing fragments. Fragments are
orders lists of atoms.
"""

from futile.Utils import write as safe_print
try:
    from collections.abc import MutableSequence
except ImportError:
    from collections import MutableSequence


def distance(i, j, cell=None):
    """
    Distance between fragments, defined as distance between center of mass

    Args:
      i (Fragment): first fragment.
      j (Fragment): second fragment.
      cell (array): an array describing the unit cell.

    Returns:
      (float): the distance between centers of mass.
    """
    import numpy
    vec = i.centroid - j.centroid
    if cell:
        per = numpy.where(cell > 0.0)
        for i, p in enumerate(per):
            vec -= cell[i] * int(round(vec[i] / cell[i]))
    return numpy.sqrt(numpy.dot(vec, vec.T))


def pairwise_distance(i, j, cell=None):
    """
    Distance between fragments, as defined by the distance between their two
    nearest atoms.

    Args:
      i (Fragment): first fragment.
      j (Fragment): second fragment.
      cell (array): an array describing the unit cell.

    Returns:
      (float): the pairwise distance between fragments.
    """
    if not cell:
        from scipy.spatial.distance import cdist
        pos_i = [at.get_position() for at in i]
        pos_j = [at.get_position() for at in j]

        dmat = cdist(pos_i, pos_j)

        dist = dmat.min()

    else:
        dist = distance(Fragment([i[0]]), Fragment([j[0]]), cell)
        for at_i in i:
            for at_j in j:
                new_d = distance(Fragment([at_i]), Fragment([at_j]), cell)
                if new_d < dist:
                    dist = new_d

    return dist


def lineup_fragment(frag):
    """
    Align the principal axis of inertia of the fragments along the
    coordinate axis. Also shift the fragment such as its centroid is zero.

    Args:
      (BigDFT.Fragments.Fragment): the fragment to transform.

    Returns:
      (BigDFT.Fragments.Fragment): the transformed fragment.
    """
    from numpy.linalg import eig

    # Shift the centroid to the origin.
    Shift = Translation(frag.centroid)
    Shift.invert()
    new_frag = Shift.dot(frag)

    # Align the principal axis of inertia are on the coordinate axis
    Imat = new_frag.ellipsoid()
    w, v = eig(Imat)
    Redress = Rotation(v.T)
    return Redress.dot(new_frag)


class RotoTranslation():
    """
    Define a transformation which can be applied to a fragment. This
    rotation is defined by giving this class two fragments, and the
    rototranslation between those fragments is automatically computed.

    Args:
      frag1 (BigDFT.Fragments.Fragment): the first position.
      frag2 (BigDFT.Fragments.Fragment): the second position.
    """
    def __init__(self, frag1, frag2):
        try:
            from BigDFT import wahba
            from numpy import matrix

            # Setup each fragment as a matrix of positions
            pos1 = matrix([at.get_position() for at in frag1])
            pos2 = matrix([at.get_position() for at in frag2])

            # Compute transformation
            self.R, self.t, self.J = wahba.rigid_transform_3D(pos1, pos2)
        except Exception as e:
            safe_print('Error', e)
            self.R, self.t, self.J = (None, None, 1.0e10)

    def dot(self, frag):
        """
        Apply the rototranslations on a fragment.

        Args:
          frag (BigDFT.Fragments.Fragment): the fragment to rototranslate.

        Return:
          (BigDFT.Fragments.Fragment): the rototranslated fragment.
        """
        from BigDFT import wahba as w
        from numpy import matrix, array
        from copy import deepcopy

        # Convert to a position matrix
        pos = matrix([at.get_position() for at in frag])

        # Apply roto-translation
        if self.t is None:
            res = w.apply_R(self.R, pos)
        elif self.R is None:
            res = w.apply_t(self.t, pos)
        else:
            res = w.apply_Rt(self.R, self.t, pos)

        # Copy back to the fragment
        newfrag = deepcopy(frag)
        for i in range(0, len(newfrag)):
            newfrag[i].set_position([float(x) for x in array(res[i, :])[0]])

        return newfrag

    def invert(self):
        """
        Computes the inverse rototranslation.
        """
        self.t = -self.t
        if self.R is not None:
            self.R = self.R.T


class Translation(RotoTranslation):
    """
    This class defines a simple translation.

    Args:
      t (list): the vector describing the translation.
    """
    def __init__(self, t):
        import numpy
        self.R = None
        self.t = numpy.mat(t).reshape(3, 1)
        self.J = 0.0


class Rotation(RotoTranslation):
    """
    This class defines a simple rotation.

    Args:
      t (list): the vector describing the rotation.
    """
    def __init__(self, R):
        self.t = None
        self.R = R
        self.J = 0.0


def interpolate_fragments(A, B, steps, extrapolation_steps=0):
    """
    Given two fragments A and B, this generates a list of Fragments
    that interpolate between A and B in a specified number of steps.

    Args:
      A (BigDFT.Fragments.Fragment): starting fragment.
      B (BigDFT.Fragments.Fragment): ending fragment.
      steps (int): the number of steps to take between A and B.
      extrapolation_steps (int): optionally, we can extrapolate a number of
        steps beyond B on the same trajectory.

    Returns:
      (list): a list of fragments interpolating between A and B including
      A and B.
    """
    from numpy import matrix, array
    from BigDFT.wahba import interpolate_points
    from copy import deepcopy

    pos1 = matrix([at.get_position() for at in A])
    pos2 = matrix([at.get_position() for at in B])

    point_list = interpolate_points(pos1, pos2, steps, extrapolation_steps)

    frag_list = []
    for i in range(0, steps+extrapolation_steps+2):
        new_frag = deepcopy(A)
        for j in range(0, len(new_frag)):
            new_frag[j].set_position([float(x)
                                      for x in array(point_list[i][j])[0]])
        frag_list.append(new_frag)

    return frag_list


class Fragment(MutableSequence):
    """
    A fragment is a list of atoms in a system. Fragment might have quantities
    associated to it, like its electrostatic multipoles (charge, dipole, etc.)
    and also geometrical information (center of mass, principla axis etc.). A
    Fragment might also be rototranslated and combined with other moieteies to
    form a :class:`BigDFT.Systems.Systems`.

    Args:
      atomlist (list): list of atomic dictionaries defining the fragment
      xyzfile (BigDFT.IO.XYZReader): an XYZ file to read from.
      posinp (dict): the posinp style dictionary from a logfile/input file.
      astruct (dict): a BigDFT atomic structure style dictionary.
      system (BigDFT.Systems.System): a BigDFT system, esssentially this
        reduces many fragments into a single fragment.

    .. todo::
       Define and describe if this API is also suitable for solid-state
       fragments

    """

    def __init__(self, atomlist=None, xyzfile=None, posinp=None, astruct=None,
                 system=None):
        from BigDFT.Atoms import Atom
        self.atoms = []

        if system is not None:
            self._system_to_fragment(system)
            return

        # insert atoms.
        if atomlist:
            for atom in atomlist:
                self.append(Atom(atom))
        elif xyzfile:
            with xyzfile:
                for line in xyzfile:
                    self.append(Atom(line))
        elif posinp:
            units = posinp.get('units', 'angstroem')
            for atom in posinp['positions']:
                self.append(Atom(atom, units=units))
        elif astruct:
            units = astruct.get('units', 'angstroem')
            rshift = astruct.get('Rigid Shift Applied (AU)', [0.0, 0.0, 0.0])
            for atom in astruct["positions"]:
                self.append(Atom(atom, units=units))
            self.translate([-1.0*x for x in rshift])

        # Values
        self.q1 = None
        self.q2 = None
        self.frozen = None
        self.conmat = None

    def __len__(self):
        return len(self.atoms)

    def __delitem__(self, index):
        self.atoms.__delitem__(index)

    def insert(self, index, value):
        from BigDFT.Atoms import Atom
        self.atoms.insert(index, Atom(value))

    def __setitem__(self, index, value):
        from BigDFT.Atoms import Atom
        self.atoms.__setitem__(index, Atom(value))

    def __getitem__(self, index):
        # If they ask for only one atom, then we return it as an atom.
        # but if it's a range we return a ``Fragment`` with those atoms in it.
        if isinstance(index, slice):
            return Fragment(atomlist=self.atoms.__getitem__(index))
        else:
            return self.atoms.__getitem__(index)

    def __add__(self, other):
        from copy import deepcopy
        rval = deepcopy(self)
        rval += other
        return rval

    def __radd__(self, other):
        return self if other == 0 else self.__add__(other)

    def __eq__(self, other):
        """
        Compare two fragments. They are equal if all the atoms of the fragments
        are identical.

        other (Fragment): the fragment to compare with.
        """
        ok = True
        for at1, at2 in zip(self, other):
            ok = ok and at1 == at2
            if not ok:
                break
        return ok

    def serialize(self, name, units='bohr'):
        """
        Transform the fragment in a list that can be employed for
        the construction of dataframes or pandas series.

        Args:
            name (str): the name of the fragment
            units (str): the units for the positions
        Returns:
            list: the serialized fragment
        """
        positions = []
        for at in self:
            atdict = {'frag': name}
            atdict.update(at.serialize(units=units))
            positions.append(atdict)
        return positions

    @property
    def centroid(self):
        """
        The center of a fragment.
        """
        from numpy import mean, ravel
        pos = [at.get_position() for at in self]
        return ravel(mean(pos, axis=0))

    def center_of_charge(self):  # , zion):
        """
        The charge center which depends both on the position and net charge
        of each atom.
        """
        from numpy import array
        cc = array([0.0, 0.0, 0.0])
        qtot = 0.0
        for at in self:
            # netcharge = at.q0
            # zcharge = zion[at.sym]
            elcharge = at.nel  # zcharge - netcharge
            cc += elcharge * array(at.get_position())
            qtot += elcharge
        return cc / qtot

    @property
    def q0(self):
        """
        Provides the global monopole of the fragments given as a sum of the
        monopoles of the atoms.
        """
        if len(self) == 0:
            return None
        return [sum(filter(None, [at.q0 for at in self]))]

    def d0(self, center=None):
        """
        Fragment dipole, calculated only from the atomic charges.

        Args:
          center (list): the center of charge of the fragment.
            If this is not present, the centroid is used.
        """
        from numpy import zeros, array
        # one might added a treatment for non-neutral fragments
        # but if the center of charge is used the d0 value is zero
        if center is not None:
            cxyz = center
        else:
            cxyz = self.center_of_charge()  # centroid

        d0 = zeros(3)
        found = False
        for at in self:
            if at.q0 is None:
                found = False
                break
            found = True
            d0 += at.q0 * (array(at.get_position()) - array(cxyz))
        if found:
            return d0
        else:
            return None

    def d1(self, center=None):
        """
        Fragment dipole including the atomic dipoles.

        Args:
          center (list): the center of charge of the fragment.
            If this is not present, the centroid is used.
        """
        from numpy import zeros
        d1 = zeros(3)
        dtot = self.d0(center)
        if dtot is None:
            return dtot

        found = False
        for at in self:
            if at.q1 is not None:
                found = True
                d1 += at.q1

        if found:
            return d1 + dtot
        else:
            return None
        pass

    def ellipsoid(self, center=0.0):
        """
        Todo: define the ellipsoid.
        """
        import numpy as np
        Imat = np.mat(np.zeros(9).reshape(3, 3))
        for at in self:
            rxyz = np.array(at.get_position()) - center
            Imat[0, 0] += rxyz[0]**2  # rxyz[1]**2+rxyz[2]**2
            Imat[1, 1] += rxyz[1]**2  # rxyz[0]**2+rxyz[2]**2
            Imat[2, 2] += rxyz[2]**2  # rxyz[1]**2+rxyz[0]**2
            Imat[0, 1] += rxyz[1] * rxyz[0]
            Imat[1, 0] += rxyz[1] * rxyz[0]
            Imat[0, 2] += rxyz[2] * rxyz[0]
            Imat[2, 0] += rxyz[2] * rxyz[0]
            Imat[1, 2] += rxyz[2] * rxyz[1]
            Imat[2, 1] += rxyz[2] * rxyz[1]
        return Imat

    def get_external_potential(self, units="bohr", charge_offset=False):
        """
        Transform the fragment information into a dictionary ready to be
        put as an external potential.

        Args:
          units (str): the units of the external potential.

        Returns:
          (dict): a dictionary describing the external potential for use in
          an input file.
        """
        pot = [at.get_external_potential(units) for at in self]

        return pot

    def get_net_force(self):
        """
        Returns the net force on a fragment in Ha/Bohr.

        Returns:
          (list) Three values which describe the net force.
        """
        from numpy import array
        ret_val = array([0.0, 0.0, 0.0])
        for at in self:
            ret_val += array(at.get_force())
        return [float(x) for x in ret_val]

    @property
    def qcharge(self):
        """
        The net charge on a fragment.
        """
        netcharge = self.q0[0] if self.q0 is not None else 0
        for at in self:
            netcharge += at.nel
        return netcharge

    @property
    def nel(self):
        """
        The number of valence electrons of the atoms of the fragment
        """
        return sum([at.nel for at in self])

    def rotate_on_axis(self, angle, axis, units="radians"):
        """
        Rotate a fragment along a specific axis.

        Args:
          angle (float): angle to rotate along the axis.
          axis (list): a list of floats defining the vector to rotate along.
          units (str): either radians or degrees.
        """
        from math import pi, cos, sin
        from numpy import mat, array
        from numpy.linalg import norm
        from copy import deepcopy

        # Deal with the units.
        if units == "degrees":
            angle_value = angle * pi/180
        elif units == "radians":
            angle_value = angle
        else:
            raise ValueError("Units must be degrees or radians")

        # Translate back to the origin.
        centroid = self.centroid
        self.translate(-1.0 * centroid)

        # Normalize the unit vector.
        vec = array(axis) / norm(axis)

        # Build the rotation
        c = cos(angle_value)
        s = sin(angle_value)
        t = 1 - c
        x = vec[0]
        y = vec[1]
        z = vec[2]

        rot = mat([
                   [t*x*x + c,    t*x*y - z*s, t*x*z + y*s],
                   [t*x*y + z*s,  t*y*y + c,   t*y*z - x*s],
                   [t*x*z - y*s,  t*y*z + x*s, t*z*z + c]
                  ])

        # Rotate
        rt = Rotation(rot)
        rself = rt.dot(self)
        for i in range(0, len(self)):
            self[i] = deepcopy(rself[i])

        # Translate back
        self.translate(centroid)

    def rotate(self, x=None, y=None, z=None, units="radians"):
        """
        Rotate the fragment.

        Args:
          x (float): angle to rotate on the x axis.
          y (float): angle to rotate on the y axis.
          z (float): angle to rotate on the z axis.
          units (str): either radians or degrees.
        """
        from math import cos, sin, pi
        from numpy import mat, identity
        from copy import deepcopy

        # Deal with the units.
        if units == "degrees":
            if x:
                xval = x * pi / 180
            if y:
                yval = y * pi / 180
            if z:
                zval = z * pi / 180
        elif units == "radians":
            xval = x
            yval = y
            zval = z
        else:
            raise ValueError("Units must be degrees or radians")

        # Translate back to the origin.
        centroid = self.centroid
        self.translate(-1.0 * centroid)

        # Build the rotation
        rot = identity(3)
        if x:
            rx = mat([
                [1.0,        0.0,       0.0],
                [0.0,        cos(xval), -sin(xval)],
                [0.0,        sin(xval),  cos(xval)]
            ])
            rot = rx.dot(rot)
        if y:
            ry = mat([
                [cos(yval), 0.0,        sin(yval)],
                [0.0,        1.0,        0.0],
                [-sin(yval), 0.0,        cos(yval)]
            ])
            rot = ry.dot(rot)
        if z:
            rz = mat([
                [cos(zval), -sin(zval), 0.0],
                [sin(zval),  cos(zval), 0.0],
                [0.0,        0.0,       1.0]
            ])
            rot = rz.dot(rot)

        # Rotate
        rt = Rotation(rot)
        rself = rt.dot(self)
        for i in range(0, len(self)):
            self[i] = deepcopy(rself[i])

        # Translate back
        self.translate(centroid)

    def translate(self, vec):
        """
        Translate the fragment along the vector provided.

        Args:
          vec (list): a list of x, y, z values describing the translation (AU).
        """
        from copy import deepcopy
        rt = Translation(vec)

        trans = rt.dot(self)
        for i in range(0, len(self)):
            self[i] = deepcopy(trans[i])

    def _system_to_fragment(self, sys):
        # Convert to one big system.
        self.__init__()
        for frag in sys.values():
            for at in frag:
                self += [at]

        # Compute the connectivity, first the lookup table.
        if sys.conmat is not None:
            counter = 0
            lookup = {}
            for fragid, frag in sys.items():
                for i, at in enumerate(frag):
                    lookup[(fragid, i)] = counter
                    counter = counter + 1

            # Now the actual connectivity matrix.
            self.conmat = {}
            for fragid, frag in sys.items():
                for i in range(0, len(frag)):
                    t1 = (fragid, i)
                    t2list = sys.conmat[fragid][i]
                    self.conmat[lookup[t1]] = {}
                    for t2, bo in t2list.items():
                        self.conmat[lookup[t1]][lookup[t2]] = bo


def _example():
    """Example of using fragments"""
    from BigDFT.IO import XYZReader, XYZWriter
    from copy import deepcopy

    safe_print("Read in an xyz file and build from a list.")
    atom_list = []
    with XYZReader("SiO") as reader:
        for at in reader:
            atom_list.append(at)
    frag1 = Fragment(atomlist=atom_list)
    for at in frag1:
        safe_print(at.sym, at.get_position())
    safe_print("Centroid", frag1.centroid)
    safe_print()

    safe_print("Build from an xyz file directly.")
    reader = XYZReader("Si4")
    frag2 = Fragment(xyzfile=reader)
    for at in frag2:
        safe_print(at.sym, at.get_position())
    safe_print()

    safe_print("We can combine two fragments with +=")
    frag3 = deepcopy(frag1)
    frag3 += frag2
    for at in frag3:
        safe_print(at.sym, at.get_position())
    safe_print("Length of frag3", len(frag3))
    safe_print()

    safe_print("Since we can iterate easily, we can also write easily.")
    with XYZWriter("test.xyz", len(frag3), "angstroem") as writer:
        for at in frag3:
            writer.write(at)
    with open("test.xyz") as ifile:
        for line in ifile:
            safe_print(line)
    safe_print()

    safe_print("We can also extract using the indices")
    safe_print(dict(frag3[0]))
    sub_frag = frag3[1:3]
    for at in sub_frag:
        safe_print(dict(at))
    safe_print()


if __name__ == "__main__":
    _example()
