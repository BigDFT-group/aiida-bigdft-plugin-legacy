"""
This module is related to the usage of BigDFT with Fragment-related Quantities.
Input as well as Logfiles might be processed with the classes and methods
provided by it.

The main two classes for this module are :class:`BigDFT.Fragments.System` and
:class:`BigDFT.Fragments.Fragment`. A System is a named collection of
fragments, and a Fragment is a list of atoms. Thus a System behaves much like
a dictionary, whereas a Fragment behaves more like a list. All of the basic
dictionary and list like operations can be applied.
"""

from futile.Utils import write as safe_print
try:
    from collections.abc import MutableMapping, MutableSequence
except ImportError:
    from collections import MutableMapping, MutableSequence


def GetFragTuple(fragid):
    """
    Fragment ids should have the form: "NAME:NUMBER". This splits the fragment
    into the name and number value.

    Args:
      fragid (str): the fragment id string.

    Return:
      (tuple): fragment name, fragment number
    """
    return (fragid.split(":")[0], fragid.split(":")[1])


def distance(i, j, cell=None):
    """
    Distance between fragments, defined as distance between center of mass

    Args:
      i (Fragment): first fragment.
      j (Fragment): second fragment.
      cell (array): an array describing the unit cell.
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
    """
    dist = distance(Fragment([i[0]]), Fragment([j[0]]), cell)
    for at_i in i:
        for at_j in j:
            new_d = distance(Fragment([at_i]), Fragment([at_j]), cell)
            if new_d < dist:
                dist = new_d
    return dist


class Lattice():
    """
    Defines the fundamental objects to deal with periodic systems
    """

    def __init__(self, vectors):
        self.vectors = vectors

    def grid(self, origin=[0.0, 0.0, 0.0], extremes=None, radius=None):
        """
        Produces a set of translation vectors from a given origin
        """
        import numpy as np
        transl = []
        g = [[], [], []]  # the grid of discrete translations
        if extremes is not None:
            # print extremes
            for i, d in enumerate(extremes):
                for k in range(d[0], d[1] + 1):
                    g[i].append(k)
            # print g
            for i in g[0]:
                arri = np.array(self.vectors[0]) * i
                for j in g[1]:
                    arrj = np.array(self.vectors[1]) * j + arri
                    for k in g[2]:
                        arrk = np.array(self.vectors[2]) * k + arrj
                        vect = np.array(origin) + arrk
                        app = True
                        if radius is not None:
                            app = np.linalg.norm(arrk) < radius
                        if app:
                            transl.append(vect)
        return transl


class RotoTranslation():
    """
    Define a transformation which can be applied to a group of atoms. This
    rotation is defined by giving this class two positions, and the
    rototranslation between those points is auomatically computed.

    Args:
      pos1 (list): the first position.
      pos2 (list): the second position.
    """
    def __init__(self, pos1, pos2):
        try:
            import wahba
            self.R, self.t, self.J = wahba.rigid_transform_3D(pos1, pos2)
        except Exception as e:
            safe_print('Error', e)
            self.R, self.t, self.J = (None, None, 1.0e10)

    def dot(self, pos):
        """
        Apply the rototranslations on the set of positions provided by pos

        Args:
          pos (list): the position to rototranslate.
        """
        import wahba as w
        if self.t is None:
            res = w.apply_R(self.R, pos)
        elif self.R is None:
            res = w.apply_t(self.t, pos)
        else:
            res = w.apply_Rt(self.R, self.t, pos)
        return res

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


class Fragment(MutableSequence):
    """
    A fragment is a list of atoms in a system. Fragment might have quantities
    associated to it, like its electrostatic multipoles (charge, dipole, etc.)
    and also geometrical information (center of mass, principla axis etc.). A
    Fragment might also be rototranslated and combined with other moieteies to
    form a :class:`BigDFT.Fragments.System`.

    Args:
      atomlist (list): list of atomic dictionaries defining the fragment
      xyzfile (XYZReader): an XYZ file to read from.
      posinp (dict): the posinp style dictionary from a logfile/input file.

    .. todo::
       Define and describe if this API is also suitable for solid-state
         fragments

    """

    def __init__(self, atomlist=None, xyzfile=None, posinp=None):
        from BigDFT.Atom import Atom
        self.atoms = []

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

        # Values
        self.purity_indicator = None
        self.q0 = None
        self.q1 = None
        self.q2 = None

    def __len__(self):
        return len(self.atoms)

    def __delitem__(self, index):
        self.atoms.__delitem__(index)

    def insert(self, index, value):
        from Atom import Atom
        self.atoms.insert(index, Atom(value))

    def __setitem__(self, index, value):
        from Atom import Atom
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

    @property
    def centroid(self):
        """
        The center of a fragment.
        """
        from numpy import mean, ravel
        pos = [at.get_position() for at in self]
        return ravel(mean(pos, axis=0))

    def center_of_charge(self, zion):
        """
        The charge center which depends both on the position and net charge
        of each atom.
        """
        from numpy import array
        cc = array([0.0, 0.0, 0.0])
        qtot = 0.0
        for at in self:
            netcharge = at.q0
            zcharge = zion[at.sym]
            elcharge = zcharge - netcharge
            cc += elcharge * array(at.get_position())
            qtot += elcharge
        return cc / qtot

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
            cxyz = self.centroid

        d0 = zeros(3)
        found = False
        for at in self:
            if at.q0 is not None:
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

    def fragment_transformation(self, frag2):
        """
        Returns the transformation among fragments if it exists.

        frag2 (BigDFT.Fragments.Fragment): the fragment to transform between.

        Returns:
          (numpy.matrix) : the 3x3 transformation matrix.
        """
        from numpy import mat
        pos1 = mat(self.centroid)
        pos2 = mat(frag2.centroid)
        return RotoTranslation(pos1, pos2)

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

    def line_up(self):
        """
        Align the principal axis of inertia of the fragments along the
        coordinate axis. Also shift the fragment such as its centroid is zero.
        """
        from numpy.linalg import eig
        Shift = Translation(self.centroid)
        Shift.invert()
        self.transform(Shift)
        # now the centroid is zero
        Imat = self.ellipsoid()
        w, v = eig(Imat)
        Redress = Rotation(v.T)
        self.transform(Redress)
        # now the principal axis of inertia are on the coordinate axis

    @property
    def qcharge(self):
        """
        The net charge on a fragment.
        """
        from Atom import nzion
        netcharge = self.q0[0]
        for at in self:
            if "nzion" in at:
                zcharge = at.qcharge["nzion"]
            else:
                zcharge = nzion(at.sym)
            netcharge += zcharge
        return netcharge

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
        self.transform(Rotation(rot))

        # Translate back
        self.translate(centroid)

    def translate(self, vec):
        """
        Translate the fragment along the vector provided.

        vec (list): a list of x, y, z values describing the translation (AU).
        """
        from numpy import mat
        trans = RotoTranslation(mat([0, 0, 0]), mat(vec))
        self.transform(trans)

    def transform(self, Rt):
        """
        Apply a rototranslation of the fragment positions.

        Rt (BigDFT.Fragments.RotoTranslation): the transformation to apply.
        """
        from numpy import mat
        for at in self:
            pos = mat(at.get_position("bohr"))
            pos2 = Rt.dot(pos)
            at.set_position(pos2.A1)


def system_from_log(log, fragmentation=None):
    """
    This function returns a :class:`~BigDFT.Fragment.System` class out of a
    logfile. If the logfile contains information about fragmentation and atomic
    multipoles, then the system is created accordingly.
    Otherwise, the fragmentation scheme is determined by the fragmentation
    variable.

    Args:
       log (Logfile): the logfile of the QM run. In general must have been done
           with Linear Scaling formalism.
       fragmentation (str): the scheme to be used for the fragmentation in the
           case if not provided internally by the logfile.
           The possible values are ``atomic`` and ``full``, in which case the
           system as as many fragments as the number of atoms, or only one
           fragment, respectively.
    Returns:
        (Fragments.System) The instance of the class containing fragments.
    """
    name = log.log.get('run_name', 'FULL') + ':0'
    posinp = log.log['posinp']
    full_system = System()
    full_system[name] = Fragment(posinp=posinp)

    # provide the atomic information on the system
    if hasattr(log, 'electrostatic_multipoles'):
        full_system.set_atom_multipoles(log)
        if hasattr(log, 'forces'):
            full_system.set_atom_forces(log)

    # now we may defragment the system according to the provided scheme
    if fragmentation == 'full':
        return full_system
    elif fragmentation == 'atomic':
        atomic_system = System()
        for iat, at in enumerate(full_system[name]):
            atomic_system['ATOM:' + str(iat)] = Fragment([at])
        return atomic_system
    else:
        frag_dict = {}
        for iat, tupl in enumerate(zip(posinp['positions'],
                                       full_system[name])):
            at, obj = tupl
            fragid = at.get('frag', 'ATOM:' + str(iat))
            if isinstance(fragid, list):
                fragid = str(fragid[0]) + ':' + str(fragid[1])
            if fragid not in frag_dict:
                frag_dict[fragid] = [obj]
            else:
                frag_dict[fragid].append(obj)
        frag_system = System()
        for fragid in frag_dict:
            frag_system[fragid] = Fragment(frag_dict[fragid])
        return frag_system


class System(MutableMapping):
    """
    A system is defined by a collection of Fragments.

    It might be given by one single fragment
    """

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))

    def dict(self):
        """
        Convert to a dictionary.
        """
        return self.store

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

    @property
    def centroid(self):
        """
        Center of mass of the system
        """
        from numpy import mean
        return mean([frag.centroid for frag in self.values()], axis=0)

    @property
    def central_fragment(self):
        """
        Returns the fragment whose center of mass is closest to the centroid

        Returns:
          (str): the name of the fragment.
          (Fragment): the fragment object
        """
        import numpy as np
        CMs = [frag.centroid for frag in self.values()]
        idx = np.argmin([np.dot(dd, dd.T) for dd in (CMs - self.centroid)])
        return self.keys()[idx], self.values()[idx]

    def get_external_potential(self, units="bohr", charge_offset=False):
        """
        Transform the system information into a dictionary ready to be
        put as an external potential.

        Args:
          units (str): the units of the external potential.
          charge_offset (bool): by default the external potential ignores the
            counter charge from the protons. Setting this to true adds the
            positive charge to the potential.
        """
        ret_dict = {"units": units}
        ret_dict["values"] = []
        for frag in self.values():
            ret_dict["values"].extend(
                frag.get_external_potential(units, charge_offset))
        ret_dict["global monopole"] = sum(
            x["q0"][0] for x in ret_dict["values"])
        return ret_dict

    def get_nearest_fragment(self, target):
        """
        Given a fragment id in a system, this computes the nearest fragment.

        Args:
          target (str): the fragment to find the nearest neighbor of.

        Returns:
          (str): the id of the nearest fragment.
        """
        from scipy.spatial import KDTree

        # Setup the KD Tree for distance lookup
        poslist = []
        frag_lookup = []
        for fragid, frag in self.items():
            if fragid == target:
                continue
            for at in frag:
                poslist.append(at.get_position())
                frag_lookup.append(fragid)
        tree = KDTree(poslist)

        # Find the nearest fragment with a query of the tree.
        ndist, nearest = tree.query([x.get_position()
                                     for x in self[target]])
        return frag_lookup[nearest[0]]

    def get_net_force(self):
        """
        Returns the net force on a system in Ha/Bohr.

        Returns:
          (list): Three values which describe the net force.
        """
        from numpy import array
        ret_val = array([0.0, 0.0, 0.0])
        for frag in self.values():
            ret_val += array(at.get_net_force())
        return [float(x) for x in ret_val]

    def get_posinp(self, units='angstroem'):
        """
        Provide the dictionary which has to be passed to the ``posinp`` value
        of the :meth:`run` method of  the
        :class:`~BigDFT.Calculators.SystemCalculator` class instance.

        Args:
           units (str): The units of the file. May be "angstroem" or "bohr".
        """
        pos = [{at.sym: at.get_position(units),
                'frag': list(GetFragTuple(frag))}
               for frag in self for at in self[frag]]
        return {'units': units, 'positions': pos}

    def plot_purity(self, ax, minval=0):
        """
        Plot the purity values of the fragments in this system.
        This assumes they have already been set.

        Args:
          axs: the axs we we should plot on.
        """
        ax.set_xlabel("Fragment", fontsize=12)
        ax.set_ylabel("Purity Values", fontsize=12)
        ax.set_xticks(range(len(self.keys())))
        ax.set_yscale("log")

        pvals = [abs(frag.purity_indicator) for frag in self.values()]
        ax.plot([v for _, v in sorted(zip(self.keys(), pvals),
                                      key=lambda x:
                                      int(GetFragTuple(x[0])[1]))
                 if v > minval], 'x', markersize=12,
                color="black")

        ax.set_xticklabels([x for x, v in sorted(zip(self.keys(), pvals),
                                                 key=lambda x:
                                                 int(GetFragTuple(x[0])[1]))
                            if v > minval], rotation=90)

    @property
    def q0(self):
        """
        Provides the global monopole of the system given as a sum of the
        monopoles of the atoms.
        """
        if len(self) == 0:
            return None
        return [sum(filter(None, [frag.q0[0] for frag in self.values()]))]

    @property
    def qcharge(self):
        """
        The total qcharge of a system.
        """
        return sum([frag.qcharge for frag in self.values()])

    def rename_fragments(self):
        """
        This procedure automatically names the fragments in a system.

        Returns:
          (System): the same system, with the automatic naming scheme.
        """
        rnsys = System()
        for i, fragid in enumerate(self):
            rnsys["FRAG:"+str(i)] = self[fragid]
        return rnsys

    def set_atom_multipoles(self, logfile, correct_charge=True):
        """
        After a run is completed, we have a set of multipoles defined on
        each atom. This routine will set those values on to each atom
        in the system.

        Args:
          logfile (Logfiles.Logfile): logfile with the multipole values.
          correct_charge (bool): currently there is an inconsistency in
            terms of gross charge, and this corrects it.
        """
        mp = logfile.electrostatic_multipoles
        for pole in mp["values"]:
            pole["units"] = mp["units"]
        lookup = self.compute_matching(mp["values"])

        # Assign
        for fragid, frag in self.items():
            for i, at in enumerate(frag):
                idx = lookup[fragid][i]
                if idx >= 0:
                    at.set_multipole(mp["values"][idx], correct_charge)

    def set_atom_forces(self, logfile):
        """
        After a run is completed, we have the forces on each atom in the
        logfile. This routine will set those values to each atom in this sytem.

        Args:
          logfile (Logfiles.Logfile): logfile with the forces.
        """
        # We will use the multipoles to help figure out which force value
        # is associated with which atom
        mp = logfile.electrostatic_multipoles
        lookup = self.compute_matching(mp["values"])

        # Assign forces
        forces = logfile.forces
        for fragid, frag in self.items():
            for i, at in enumerate(frag):
                idx = lookup[fragid][i]
                if idx >= 0:
                    at.set_force(forces[idx].values()[0])

    def write_fragfile(self, filename, logfile):
        """
        Write out the file needed as input to the fragment analysis routine.

        Args:
          filename (str): name of the file to write to.
          logfile (Logfiles.Logfile): the log file this calculation is based
            on (to ensure a matching order of atoms).
        """
        from yaml import dump

        # Extract the indices
        mp = logfile.electrostatic_multipoles
        lookup = self.compute_matching(mp["values"])
        outlist = []
        for fragid, frag in self.items():
            outlist.append([])
            for i, at in enumerate(frag):
                idx = lookup[fragid][i]
                outlist[-1].append(idx + 1)

        # Write
        with open(filename, "w") as ofile:
            dump(outlist, ofile)

    def compute_matching(self, atlist, shift=None):
        """
        Frequently we are passed a list of atom like objects from which we
        need to extract data and assign it to a system. However, a system
        can potentially store those atoms in any order, and may not have
        the same set of atoms. This helper routine creates a mapping between
        this list view, to the dictionary view of the system class.

        Args:
          atlist (list): a list of atom like objects.
          shift (list): if the positions in atlist are shifted by some constant
            vector you can specify that here.

        Returns:
          (dict): a mapping from a system to indices in the atom list. If
            an atom is not in the list, an index value of -1 is assigned.
        """
        from Atom import Atom
        from numpy import array
        from scipy.spatial import KDTree

        # Convert everything to pure positions to avoid overhead.
        poslist = [array(Atom(x).get_position("bohr")) for x in atlist]
        if shift is not None:
            poslist = [x - array(shift) for x in poslist]
        tree = KDTree(poslist)

        # Seach for the mapping values
        mapping = {}
        for fragid, frag in self.items():
            mapping[fragid] = []
            for at in frag:
                atpos = array(at.get_position("bohr"))
                ndist, nearest = tree.query(atpos)
                mapping[fragid].append(nearest)

        return mapping


if __name__ == "__main__":
    from XYZ import XYZReader, XYZWriter
    from os.path import join
    from os import system
    from copy import deepcopy

    safe_print("Read in an xyz file and build from a list.")
    atom_list = []
    with XYZReader(join("Database", "XYZs", "SiO.xyz")) as reader:
        for at in reader:
            atom_list.append(at)
    frag1 = Fragment(atomlist=atom_list)
    for at in frag1:
        safe_print(at.sym, at.get_position())
    safe_print("Centroid", frag1.centroid)
    safe_print()

    safe_print("Build from an xyz file directory.")
    reader = XYZReader(join("Database", "XYZs", "Si4.xyz"))
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
    system("cat test.xyz")
    safe_print()

    safe_print("We can also extract using the indices")
    safe_print(dict(frag3[0]))
    sub_frag = frag3[1:3]
    for at in sub_frag:
        safe_print(dict(at))
    safe_print()

    safe_print("Now we move on to testing the system class.")
    safe_print("We might first begin in the easiest way.")
    sys = System(frag1=frag1, frag2=frag2)
    for at in sys["frag1"]:
        safe_print(dict(at))
    for at in sys["frag2"]:
        safe_print(dict(at))
    safe_print()

    safe_print("What if we want to combine two fragments together?")
    sys["frag1"] += sys.pop("frag2")
    for at in sys["frag1"]:
        safe_print(dict(at))
    safe_print("frag2" in sys)
    safe_print()

    safe_print("What if I want to split a fragment by atom indices?")
    temp_frag = sys.pop("frag1")
    sys["frag1"], sys["frag2"] = temp_frag[0:3], temp_frag[3:]
    for at in sys["frag1"]:
        safe_print(dict(at))
    for at in sys["frag2"]:
        safe_print(dict(at))
    safe_print()

    safe_print("Construct a system from an XYZ file.")
    fname = join("Database", "XYZs", "BH2.xyz")
    sys2 = System(frag1=Fragment(xyzfile=XYZReader(fname)))

    safe_print("Split it to fragments")
    sys2["frag1"], sys2["frag2"] = sys2["frag1"][0:1], sys2["frag1"][1:]

    safe_print("And write to file")
    with XYZWriter("test.xyz", len(frag3), "angstroem") as writer:
        for fragid, frag in sys2.items():
            for at in frag:
                writer.write(at)
    system("cat test.xyz")
