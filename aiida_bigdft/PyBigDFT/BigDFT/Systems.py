"""
This module contains the System level class of PyBigDFT. Systems are named
collections of fragments, and represent a complete system for simulation.
"""
from futile.Utils import write as safe_print

try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping


class System(MutableMapping):
    """
    A system is defined as a named collection of fragments. You can manipulate
    a system as if it were a standard python dictionary, however it also has
    helper routines for performing operations on the full system.
    """

    def __init__(self, *args, **kwargs):
        from BigDFT.UnitCells import UnitCell
        self.store = dict()
        self.update(dict(*args, **kwargs))
        self.conmat = None
        self.cell = UnitCell()

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

    def __eq__(self, other):
        """
        Compare two systems. They are equal if all the fragments of
            the systems are identical.

        other (System): the fragment to compare with
        """
        ok = True
        for fr1 in self:
            ok = ok and self[fr1] == other[fr1]
            if not ok:
                break
        return ok

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
        return list(self.keys())[idx], list(self.values())[idx]

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
        # ret_dict["global monopole"] = sum(
        #     x["q0"][0] for x in ret_dict["values"])
        return ret_dict

    def get_k_nearest_fragments(self, target, k, cutoff=None,
                                return_type='list'):
        """
        Given a fragment id in a system, this computes the nearest fragment.

        Args:
          target (str): the fragment to find the nearest neighbor of.
          k (int): the number of fragments to look for.
          cutoff (float): will only return fragments with a certain cutoff.
          return_type (str): 'list' or 'dict'
        Returns:
          (list, dict): the ids of the nearest fragments, or their distances
              as values in the case in which 'dict' is provided in the
              `return_type` argument
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

        # Find the nearest fragments with a query of the tree.
        targetpost = [x.get_position() for x in self[target]]
        if cutoff is not None:
            ndist, nearest = tree.query(targetpost, k=k,
                                        distance_upper_bound=cutoff)
        else:
            ndist, nearest = tree.query(targetpost, k=k)

        # We now have the nearest atom to each atom in this fragment.
        # Next we combine this information and extract the closest
        # fragments.

        if k == 1:
            ndist = [ndist]
            nearest = [nearest]

        distdict = {}
        for i in range(0, len(nearest)):
            for idx, dist in zip(nearest[i], ndist[i]):
                try:
                    fragidx = frag_lookup[idx]
                except IndexError:
                    # kdtree returns invalid indices if it can't find enough
                    # points.
                    continue
                if fragidx not in distdict:
                    distdict[fragidx] = dist
                elif distdict[fragidx] > dist:
                    distdict[fragidx] = dist

        # Extract the k smallest values.
        minlist = []
        mindict = {}
        for i in range(0, k):
            if len(distdict) == 0:
                break
            key = min(distdict, key=distdict.get)
            minlist.append(key)
            mindict[key] = distdict[key]
            del distdict[key]

        if return_type == 'list':
            return minlist
        elif return_type == 'dict':
            return mindict

    def get_nearest_fragment(self, target):
        """
        Given a fragment id in a system, this computes the nearest fragment.

        Args:
          target (str): the fragment to find the nearest neighbor of.

        Returns:
          (str): the id of the nearest fragment.
        """
        return self.get_k_nearest_fragments(target, k=1)[0]

    def get_net_force(self):
        """
        Returns the net force on a system in Ha/Bohr.

        Returns:
          (list): Three values which describe the net force.
        """
        from numpy import array
        ret_val = array([0.0, 0.0, 0.0])
        for frag in self.values():
            ret_val += array(frag.get_net_force())
        return [float(x) for x in ret_val]

    def get_posinp(self, units='angstroem'):
        """
        Provide the dictionary which has to be passed to the ``posinp`` value
        of the :meth:`run` method of  the
        :class:`~BigDFT.Calculators.SystemCalculator` class instance.

        Args:
           units (str): The units of the file. May be "angstroem" or "bohr".
        """
        pos = []
        for fragid, frag in self.items():
            for at in frag:
                atdict = {at.sym: at.get_position(units)}
                atdict["frag"] = list(GetFragTuple(fragid))
                if frag.frozen:
                    atdict["Frozen"] = frag.frozen
                pos.append(atdict)

        result = {}
        result["units"] = units
        result["positions"] = pos
        result["cell"] = self.cell.get_posinp(units)

        return result

    def serialize(self, order=None, units='bohr'):
        """
        Transform the system in a list that can be employed for
        the construction of dataframes or pandas series.
        Args:
            order (list): list of fragments to serialize in order
            units (str): the units for the positions
        Returns:
            list: the serialized system as well as a lookup dictionary
                that contains the order of the fragment atoms in the series
        """
        if order is None:
            order = self
        positions = []
        lookup = {}
        iat = 0
        for key in order:
            frag = self[key]
            lookup[key] = [iat+i for i in range(len(frag))]
            iat += len(frag)
            positions += frag.serialize(name=key, units=units)
        return positions

    @property
    def df(self):
        if not hasattr(self, '_df'):
            self._df = self.to_dataframe()
        return self._df

    def to_dataframe(self, **kwargs):
        """
        Convert the system into a dataframe, from the
        `py:meth:~System.serialize` method.

        Args:
            **kwargs: arguments to be passed to
                `py:meth:System.serialize` method
        """
        from pandas import DataFrame as DF
        df = DF(self.serialize(**kwargs))
        validate_dataframe_representation(self, df)
        return df

    def dataframe_slicing(self, df=None):
        """
        Define a dictionaries of two-element lists indicating the slicing
        (start and end points) of the corresponding fragment in the order
        provided by the dataframe

        Args:
            df (Dataframe): associated to the system

        Returns:
            dict:  dictionary of the slicings of the fragments in the dataframe
        """
        if df is None:
            df = self.df
        sl = fragment_slicing_from_system_dataframe(df)
        check_slicing(self, sl)
        return sl

    @property
    def PointParticles(self):
        """
        Transform the system into a `py:class:~PointParticles.PointParticles`
        object
        """
        from BigDFT.PointParticles import PointParticles as PP
        if not hasattr(self, '_PP'):
            self._PP = PP(**point_particle_objects(self.df))
        return self._PP

    @property
    def electrostatic_interactions(self):
        """Dictionary of the Electrostatic interactions between fragments.
        """
        return self.PointParticles.Eel_dict(self.dataframe_slicing())

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

    def rename_fragments(self, fragment_mapping=None):
        """
        This procedure automatically names the fragments in a system.

        Args:
           fragment_mapping (dict) : Dictionary containing list
                 of fragments that are additionally added to each
                 of the original system's fragments.

        Returns:
          BigDFT.Systems.System: the same system, with the automatic naming
          scheme.
        """
        rnsys = System()
        for i, fragid in enumerate(self):
            if fragment_mapping is None:
                tuplek = "FRAG:"+str(i)
            else:
                additional_list = fragment_mapping[fragid]
                if not isinstance(additional_list, list):
                    tuplek = additional_list
                elif len(additional_list) == 0:
                    tuplek = fragid
                else:
                    # tuplek = (fragid, ) + tuple((f for f in additional_list))
                    tuplek = fragid + '+'+'+'.join(
                             [f for f in additional_list])
            rnsys[tuplek] = self[fragid]
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
        from BigDFT.Fragments import Fragment
        posinp = logfile.log.get('posinp')
        if posinp is not None and not isinstance(posinp, str):
            atlist = Fragment(posinp=posinp)
        else:
            atlist = Fragment(astruct=logfile.astruct)
        lookup = self.compute_matching(atlist)

        # Assign forces
        try:
            forces = logfile.forces
        except AttributeError:
            forces = logfile.astruct["forces"]
        for fragid, frag in self.items():
            for i, at in enumerate(frag):
                idx = lookup[fragid][i]
                if idx >= 0:
                    at.set_force(list(forces[idx].values())[0])

    def compute_matching(self, atlist, shift=None, check_matching=True):
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
          check_matching (bool): if set to True, this will raise an error
            if we can't match all of the atoms in the system.

        Returns:
          (dict): a mapping from a system to indices in the atom list. If
            an atom is not in the list, an index value of -1 is assigned.
        """
        from BigDFT.Atoms import Atom
        from numpy import array
        from scipy.spatial import KDTree

        # Convert everything to pure positions to avoid overhead.
        poslist = [array(Atom(x).get_position("bohr", cell=self.cell))
                   for x in atlist]
        if shift is not None:
            poslist = [x - array(shift) for x in poslist]
        tree = KDTree(poslist)

        # Seach for the mapping values
        mapping = {}
        for fragid, frag in self.items():
            mapping[fragid] = []
            for at in frag:
                atpos = array(at.get_position("bohr", cell=self.cell))
                ndist, nearest = tree.query(atpos)
                if check_matching and ndist > 0.01:
                    raise ValueError("Unable to match atom" + str(dict(at)))
                mapping[fragid].append(nearest)

        return mapping

    def display(self, colordict=None, field_vals=None, cartoon=False):
        """
        Display the system using the inline visualizer of Py3DMol
        """
        from BigDFT.Visualization import InlineVisualizer
        viz = InlineVisualizer(400, 300)
        viz.display_system(self, colordict=colordict, field_vals=field_vals,
                           cartoon=cartoon)
        return viz

    def set_electrons_from_log(self, log):
        """
        BigDFT uses pseudopotentials, so this will extract from a logfile
        the number of actual electrons modelled for each atom in the system.

        Args:
           log (Logfiles.Logfile) : A BigDFT run logfile.
        """
        electrons = {}
        for key in log.log:
            if "psppar" not in key:
                continue
            electrons[key.split(".")[1]] = log.log[key]['No. of Electrons']

        for frag in self.values():
            for at in frag:
                at.nel = electrons[at.sym]

    def set_logfile_info(self, log):
        """
        Include the information of the logfile in the fragment quantities.

        Args:
           log (Logfiles.Logfile) : A BigDFT run logfile.
        """
        # provide the atomic information on the system
        if hasattr(log, 'electrostatic_multipoles'):
            self.set_atom_multipoles(log)
        if hasattr(log, 'forces'):
            self.set_atom_forces(log)
        self.set_electrons_from_log(log)

    def ase_potential_energy(self, ase_calculator):
        """
        Given a ASE calculator, calculates the potential energy
        of the system.

        Args:
            ase_calculator (ase.calculators.calculator.Calculator): ASE
              calculator.

        Returns:
            float: the potential energy, in Hartree
        """
        from Fragments import Fragment
        from ASEInterop import ase_potential_energy as asepot
        bigfrag = Fragment(system=self)
        return asepot(bigfrag, ase_calculator)

    def reform_superunits(self, mapping):
        """
        Creates a new system from the provided mapping

        Args:
            mapping(dict): dictionary of the form {newfrag: [frag1,frag2]}
                defining the mapping between the old fragments and the new

        Returns:
            BigDFT.Systems.System: a new system from the remapping
        """
        newsys = System()

        for newfrag, fraglist in mapping.items():
            newsys[newfrag] = sum(self[f] for f in fraglist)

        return newsys

    def update_positions_from_dict(self, posinp):
        """
        Update the atomic positions of a system from a posinp dictionary.

        This method only works if the order of atoms match.

        Args:
            posinp (dict): a posinp dictionary.
        """
        from BigDFT.Atoms import Atom
        units = posinp.get("units", "angstroem")
        i = 0
        for frag in self.values():
            for at in frag:
                at2 = Atom(posinp["positions"][i], units=units)
                at.set_position(at2.get_position())
                i += 1


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


def fragment_slicing_from_system_dataframe(df):
    """
    Define the slicing tuple needed to identify the fragment blocks
    into a system dataframe

    Args:
        df (Dataframe): the system dataframe

    Returns:
        dict: dictionary of the slicing obtained in form of  [start,end] list
    """
    current = df['frag'][0]
    slicing = {current: [0]}
    for i, f in enumerate(df['frag']):
        if f not in slicing:
            slicing[current].append(i)
            slicing[f] = [i]
            current = f
    slicing[current].append(i+1)
    return slicing


def check_slicing(system, slicing):
    """
    Assert the validity of a system's slicing
    """
    for frag, ss in slicing.items():
        assert len(system[frag]) == (ss[1]-ss[0])


def point_particle_objects(df):
    """
    Return the dictionary of the point particle quantities from a
    Systems' dataframe
    """
    if 'x_coord' in df:
        X = df[['x_coord', 'y_coord', 'z_coord']].to_numpy()
    else:
        X = df[['r_0', 'r_1', 'r_2']].to_numpy()
    P = df[['q1_2', 'q1_0', 'q1_1']].to_numpy()
    if 'qel_0' not in df:
        from BigDFT.Atoms import _nzion_default_psp as nzion
        from numpy import array
        df['zion'] = array([nzion[s] for s in df['sym']])
        # we assume net monopole definition
        df['qel_0'] = df['q0_0'] - df['nel']
    Q = df[['qel_0']].to_numpy()
    Z = df[['nel']].to_numpy()
    return {'X': X, 'P': P, 'Q': Q, 'Z': Z}


def _get_coordinate_keys(df):
    """Keys associated to atom coordinate in system's dataframe.
    """
    if 'x_coord' in df:
        return ['x_coord', 'y_coord', 'z_coord']
    else:
        return ['r_0', 'r_1', 'r_2']


def validate_dataframe_representation(sys, df):
    """
    Control if the units of the positions in the dataframe are in Bohr.
    Update the coordinate values if this is not so.

    Args:
        sys (BigDFT.System): System which provides the reference positions
        df (pandas.DataFrame): system's dataframe
    """
    import numpy as np
    from BigDFT.Atoms import AU_to_A
    from pandas import DataFrame
    assert all(df['units'] == 'bohr'), 'Dataframe units should be in Bohr'
    sl = sys.dataframe_slicing(df)
    keys = _get_coordinate_keys(df)
    all_coords = df[keys].to_numpy()
    for frag in sl:
        ist, ien = sl[frag]
        dfarr = all_coords[ist:ien]
        for at, coords in zip(sys[frag].atoms, dfarr):
            rxyz = np.array(at.get_position('bohr'))
            delta = rxyz - coords
            if np.linalg.norm(delta) > 1.e-3:  # this means angstroem
                dfarr /= AU_to_A
                delta = rxyz - coords
            assert np.linalg.norm(delta) <= 1.e-3, 'Atoms are not equivalent'
    df.update(DataFrame(all_coords, columns=keys))


class FragmentView():
    """
    The representation of a system in terms of fragments and
    groups of superunits.

    Args:
        purities (dict): dictionary of the purities of the system
        bond_order (dict): double dictionary of the bond orders of
            the superunits
        charges (dict): dictionary of the number of the electrons of
            the superunits
    """
    def __init__(self, purities, bond_orders, charges):
        self.purities = purities
        self.bond_orders = bond_orders
        self.charges = charges

    def __deepcopy__(self, memo):
        """
        Here we manually override deepcopy for performance reasons.
        """
        new_bo = {}
        new_pv = {}
        new_charges = {}
        for f1 in self.bond_orders:
            new_bo[f1] = {}
            for f2 in self.bond_orders[f1]:
                new_bo[f1][f2] = self.bond_orders[f1][f2]
            new_pv[f1] = self.purities[f1]
            new_charges[f1] = self.charges[f1]

        copy = FragmentView(new_pv, new_bo, new_charges)
        memo[id(self)] = copy
        return copy

    def refragment(self, mapping):
        newp, newbo = update_purity_and_bo(mapping, self.purities,
                                           self.bond_orders, self.charges)
        newchg = {}
        for refrag, remap in mapping.items():
            newchg[refrag] = sum(self.charges[f] for f in remap)
        return FragmentView(newp, newbo, newchg)

    def remove_fragment(self, fragid):
        """
        Remove a particular fragment from this view.

        Args:
          fragid (str): the id of the fragment to remove.
        """
        self.purities.pop(fragid)
        self.charges.pop(fragid)
        self.bond_orders.pop(fragid)

        for f in self.bond_orders:
            self.bond_orders[f].pop(fragid)


def select_from_view(view, targets):
    """
    Identify the fragments of the view that contain
    at least one of the targets

    Args:
        view (dict): if present, identifies the fragments that contain the
            relevant units
        targets (list): list of the fragments to search in the view
    Returns:
        list: fragments to select
    """
    reselect = []
    for frag in targets:
        for key, val in view.items():
            if frag in val and key not in reselect:
                reselect.append(key)
                break
    return reselect


def system_from_dict_positions(posinp, units='angstroem'):
    """
    Build a system from a set of positions from a dictionary whose yaml
    serialisation is compliant with the BigDFT yaml position format

    Args:
       posinp (list): list of the atomic specifications
    Returns:
       BigDFT.Systems.System: an instance of the system class.
          The employed fragment specification is specified in the file.
    """
    from BigDFT.Atoms import Atom
    from BigDFT.Fragments import Fragment
    sys = System()
    for iat, at in enumerate(posinp):
        frag = GetFragId(at, iat)
        if frag not in sys:
            sys[frag] = Fragment()
        sys[frag].append(Atom(at, units=units))
    return sys


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
        (BigDFT.Systems.System): The instance of the class containing
        fragments.
    """
    from BigDFT.Fragments import Fragment
    name = log.log.get('run_name', 'FULL') + ':0'

    full_system = System()
    posinp = log.log.get('posinp')
    if posinp is not None and not isinstance(posinp, str):
        full_system[name] = Fragment(posinp=posinp)
    else:
        full_system[name] = Fragment(astruct=log.astruct)

    full_system.set_logfile_info(log)

    # now we may defragment the system according to the provided scheme
    if fragmentation == 'full':
        return full_system
    elif fragmentation == 'atomic' or 'posinp' not in log.log:
        atomic_system = System()
        for iat, at in enumerate(full_system[name]):
            atomic_system['ATOM:' + str(iat)] = Fragment([at])
        return atomic_system
    else:
        posinp = log.log['posinp']
        frag_dict = {}
        for iat, tupl in enumerate(zip(posinp['positions'],
                                       full_system[name])):
            at, obj = tupl
            fragid = at.get('frag', 'ATOM:' + str(iat))
            if isinstance(fragid, list):
                fragid = ':'.join(map(str, fragid))
            if fragid not in frag_dict:
                frag_dict[fragid] = [obj]
            else:
                frag_dict[fragid].append(obj)
        frag_system = System()
        for fragid in frag_dict:
            frag_system[fragid] = Fragment(frag_dict[fragid])
        return frag_system


def plot_fragment_information(axs, datadict, colordict=None):
    """
    Often times we want to plot measures related to the different fragments
    in a system. For this routine, you can pass a dictionary mapping
    fragment ids to some kind of value. This routine takes care of the
    formatting of the axis to easily read the different fragment names.

    Args:
      axs (matplotlib.Axes): an axes object to plot on.
      datadict (dict): a dictionary from fragment ids to some kind of data
        value.
      colordict (dict): optionally, a dictionary from fragment ids to a
        color value.
    """
    # Sort by fragment id
    slabels = sorted(datadict.keys(),
                     key=lambda x: int(GetFragTuple(x)[1]))
    svalues = [datadict[x] for x in slabels]

    # Label the axis by fragments
    axs.set_xlabel("Fragment", fontsize=12)
    axs.set_xticks(range(len(datadict.keys())))
    axs.set_xticklabels(slabels, rotation=90)

    # Plot the actual values
    axs.plot(svalues, 'x', markersize=12, color='k')

    # Plot the colored values.
    if colordict:
        for i, key in enumerate(slabels):
            if key not in colordict:
                continue
            axs.plot(i, svalues[i], 'x', markersize=12, color=colordict[key])


def GetFragId(atdict, iat):
    """
    Obtain the fragment identifications from the atom description

    Args:
      atdict(dict): dictionary of the atom
      iat (int): position of the atom in the list

    Returns:
      str: fragment_id
    """
    fragid = atdict.get('frag', 'ATOM:' + str(iat))

    if all(':' in elem for elem in fragid):
        fragid = '+'.join(map(str, fragid))
    elif len(fragid) == 2 and ':' in str(fragid[1]):
        fragid = '-'.join(map(str, fragid))
    elif isinstance(fragid, list):
        fragid = ':'.join(map(str, fragid))
    return fragid


def GetFragTuple(fragid):
    """
    Fragment ids should have the form: "NAME:NUMBER" or "NAME-NUMBER". This
    splits the fragment into the name and number value.

    Args:
      fragid (str): the fragment id string.

    Return:
      (tuple): fragment name, fragment number
    """
    markers = [":", "-", "+"]
    if all(x not in fragid for x in markers):
        raise ValueError("Invalid format for fragment ID")
    if '+' in fragid:
        return tuple(fragid.split("+"))
    # elif '-' in fragid:
    #     return tuple(fragid.split("-"))
    elif isinstance(fragid, str):
        return tuple(fragid.split(":"))
    else:
        return tuple(fragid)


def copy_bonding_information(sys1, sys2):
    """
    This routine will take the bonding information of sys1 and copy it
    over to sys2. This routine requires that both systems have the exact
    same atoms in them. This is useful if you refragment a system and
    want to fix the bonding information.

    Args:
      sys1 (BigDFT.Systems.System): the system to copy bonding information
        from.
      sys2 (BigDFT.Systems.System): the system to copy information to.
    """
    # Check parameters
    if sys1.conmat is None:
        raise ValueError("The connectivity matrix of the first argument " +
                         "must be set")

    # Create at atom level connectivity matrix
    lookup = {}
    j = 0
    for fragid, frag in sys1.items():
        for i, at in enumerate(frag):
            lookup[(fragid, i)] = j
            j += 1

    atom_connectivity = {}
    for fid1, conn in sys1.conmat.items():
        for ati, links in enumerate(conn):
            i = lookup[(fid1, ati)]
            for (fid2, atj), v in links.items():
                j = lookup[(fid2, atj)]
                atom_connectivity[(i, j)] = v

    # Generate the matching
    atlist = []
    for fragid, frag in sys1.items():
        for at in frag:
            atlist.append(at)
    matching = sys2.compute_matching(atlist)
    rlookup = {}
    for fragid, frag in matching.items():
        for i, v in enumerate(frag):
            rlookup[matching[fragid][i]] = (fragid, i)

    # Copy over
    sys2.conmat = {}
    for fragid, frag in sys2.items():
        sys2.conmat[fragid] = []
        for i, at in enumerate(frag):
            sys2.conmat[fragid].append({})

    for (i, j), v in atom_connectivity.items():
        fragid1, at1 = rlookup[i]
        fragid2, at2 = rlookup[j]
        sys2.conmat[fragid1][at1][(fragid2, at2)] = v


def update_purity_and_bo(mapping, purity, bo, charges):
    """
    When merging fragments together, you will need to update the bond
    orders and purity values. This can be done by using
    `run_compute_purity` and `compute_bond_orders`, but this process is
    potentially slow. In this function, we use the old bond orders and
    purity values for update, which is much faster.

    Arguments:
        mapping (dict): a dictionary where the keys are the new fragments
          and the values are a list of old fragments that make up this
          new fragment.
        purity (dict): the old purity values of each fragment.
        bo (dict): the old bond orders of each fragment.
        charges (dict): the charges of each fragment (the sum of the
          number of electrons).

    Returns:
        (dict, dict): the new purity and bond order values.
    """
    # Purity loop
    new_purity = {}
    for f in mapping:
        new_purity[f] = 0
        old = mapping[f]

        # Sum the bond orders
        for o1 in old:
            for o2 in old:
                if o1 != o2:
                    new_purity[f] += bo[o1][o2]

        # Sum the un-normalized atomic purity values
        for o1 in old:
            new_purity[f] += purity[o1] * charges[o1] / 2.0

        # Normalize
        new_charge = sum([charges[x] for x in old])
        new_purity[f] *= 2.0 / new_charge

    # Bond order loop
    new_bo = {}
    for f1 in mapping:
        old1 = mapping[f1]
        new_bo[f1] = {}
        for f2 in mapping:
            old2 = mapping[f2]
            new_bo[f1][f2] = 0

            for o1 in old1:
                for o2 in old2:
                    new_bo[f1][f2] += bo[o1][o2]

    return new_purity, new_bo


def _example():
    """Example of using a system"""
    from BigDFT.IO import XYZReader
    from BigDFT.Fragments import Fragment

    safe_print("Read in some files for the fragments..")
    reader = XYZReader("SiO")
    frag1 = Fragment(xyzfile=reader)
    reader = XYZReader("Si4")
    frag2 = Fragment(xyzfile=reader)

    safe_print("Now we move on to testing the system class.")
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


if __name__ == "__main__":
    _example()
