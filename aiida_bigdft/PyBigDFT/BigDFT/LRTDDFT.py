
import numpy as np
from futile.Utils import write
HaeV = 27.21138386


def _occ_and_virt(log):
    """
    Extract the number of occupied and empty orbitals from a logfile
    """
    norb = log.log['Total Number of Orbitals']
    if log.log['Spin treatment'] == 'Averaged':
        norbv = log.evals[0].info[0]-norb
        return (norb,), (norbv,)
    elif log.log['Spin treatment'] == 'Collinear':
        mpol = log.log['dft']['mpol']
        norbu = int((norb+mpol)/2)
        norbd = norb-norbu
        norbvu = log.evals[0].info[0]-norbu
        norbvd = log.evals[0].info[0]-norbd
        return (norbu, norbd), (norbvu, norbvd)
    else:
        raise ValueError('Information for the orbitals to be implemented')


def transition_indexes(np, nalpha, indexes):
    """
    Returns the list of the indices in the bigdft convention that correspond
    to the couple iorb-ialpha with given spin.

    Args:
        np (tuple): (norbu,norbd) occupied orbitals: when of length 1 assumed
          spin averaged
        nalpha (tuple): (norbu, norbd)virtual orbitals: when of length 1
          assumed spin averaged
        indexes (list): list of tuples of (iorb,ialpha,ispin) desired indices
          in python convention (start from 0)
    """
    nspin = len(np)
    inds = []
    for iorb, ialpha, ispin in indexes:
        jspin = ispin if nspin == 2 else 0
        ind = ialpha+iorb*nalpha[jspin]  # local index of the spin subspace
        if ispin == 1:
            ind += np[0]*nalpha[0]  # spin 2 comes after spin one
        inds.append(ind)
    return inds


def _collection_indexes(np, nvirt_small):
    harvest = []
    for ispin in [0, 1]:
        jspin = ispin if len(np) == 2 else 0
        for ip in range(np[jspin]):
            for ialpha in range(nvirt_small[jspin]):
                harvest.append([ip, ialpha, ispin])
    return harvest


def _collection_indexes_iocc(iocc, nvirt, spin=None):
    """
    For each iocc and a selected spin provide the list of couples that are
    concerned up to nvirt
    If spin is none provide the list for all values of the spin
    """
    harvest = []
    for ispin in [0, 1]:
        jspin = ispin if len(nvirt) == 2 else 0
        if spin is not None and ispin != spin:
            continue
        for ialpha in range(nvirt[jspin]):
            harvest.append([iocc, ialpha, ispin])
    return harvest


class TransitionMatrix(np.matrix):
    """
    Matrix of Transition Quantities, might be either :class:`CouplingMatrix`
    or :class:`TransitionMultipoles`

    Args:
        matrix (matrix-like): data of the coupling matrix. If present also
          the number of orbitals should be provided.
        norb_occ (tuple): number of occupied orbitals per spin channnel.
          Compulsory if ``matrix`` is specified.
        norb_virt (tuple): number of empty orbitals per spin channnel.
          Compulsory if ``matrix`` is specified.
        log (Logfile): Instance of the logfile from which the coupling matrix
          calculation is performed. Automatically retrieves the ``norb_occ``
          and `norb_virt`` parameters. When ``log`` parameter is present the
          parameter ``matrix`` is ignored.

    Raises:
        ValueError: if the file of the coupling matrix indicated by ``log``
          does not exists
    """
    def __new__(cls, matrix=None, norb_occ=None, norb_virt=None, log=None):
        """
        Create the object from the arguments and return the ``self`` instance
        """
        import os
        if log is not None:
            datadir = log.log.get('radical', '')
            datadir = 'data-'+datadir if len(datadir) > 0 else 'data'
            cmfile = os.path.join(log.srcdir, datadir, cls._filename)
            if not os.path.isfile(cmfile):
                raise ValueError('The file "'+cmfile+'" does not exist')
            norb, norbv = _occ_and_virt(log)
            write('Loading data with ', norb, ' occupied and ',
                  norbv, ' empty states, from file "', cmfile, '"')
            try:
                import pandas as pd
                write('Using pandas:')
                mat = pd.read_csv(cmfile, delim_whitespace=True, header=None)
            except ImportError:
                write('Using numpy:')
                mat = np.loadtxt(cmfile)
            write('done')
        else:
            mat = matrix
        return super(TransitionMatrix, cls).__new__(cls, mat)

    def __init__(self, *args, **kwargs):
        """
        Perform sanity checks on the loaded matrix
        """
        log = kwargs.get('log')
        if log is not None:
            self.norb_occ, self.norb_virt = _occ_and_virt(log)
        else:
            self.norb_occ = kwargs.get('norb_occ')
            self.norb_virt = kwargs.get('norb_virt')
        assert(self.shape[0] == self._total_transitions())
        write("Shape is conformal with the number of orbitals")
        self._sanity_check()

    def _total_transitions(self):
        ntot = 0
        for no, nv in zip(self.norb_occ, self.norb_virt):
            ntot += no*nv
        if len(self.norb_occ) == 1:
            ntot *= 2
        return ntot

    def _subindices(self, norb_occ, norb_virt):
        for i, (no, nv) in enumerate(zip(norb_occ, norb_virt)):
            assert(no <= self.norb_occ[i] and nv <= self.norb_virt[i])
        harvest = _collection_indexes(norb_occ, norb_virt)
        return np.array(transition_indexes(norb_occ, self.norb_virt, harvest))

    def _sanity_check(self):
        pass


class CouplingMatrix(TransitionMatrix):
    """
    Casida Coupling Matrix, extracted from the calculation performed by BigDFT
    """
    _filename = 'coupling_matrix.txt'

    def _sanity_check(self):
        write('Casida Matrix is symmetric',
              np.allclose(self, self.T, atol=1.e-12))

    def subportion(self, norb_occ, norb_virt):
        """Extract a subportion of the coupling matrix.

        Returns a Coupling Matrix which is made by only considering the first
        ``norb_occ`` and ``norb_virt`` orbitals

        Args:
           norb_occ (tuple): new number of occupied orbitals. Must be lower
             that the instance value
           norb_virt (tuple): new number of virtual orbitals. Must be lower
             that the instance value
        """
        inds = self._subindices(norb_occ, norb_virt)
        mat = np.array([row[0, inds] for row in self[inds]])
        return CouplingMatrix(matrix=mat, norb_occ=norb_occ,
                              norb_virt=norb_virt)

    def diagonalize(self):
        """
        Diagonalize the Coupling Matrix

        Returns:
            (np.matrix, np.matrix):
            tuple of the Eigenvvalues and Eigenvectors of the coupling matrix,
              as returned by :meth:`numpy.linalg.eigh`. We perform the
              transpose of the matrix with eigenvectors to have them sorted as
              row vectors
        """
        write('Diagonalizing Coupling matrix of shape', self.shape)
        E2, C_E2 = np.linalg.eigh(self)
        write('Eigensystem solved')
        C_E2 = C_E2.T
        return E2, C_E2


class TransitionMultipoles(TransitionMatrix):
    """
    Transition dipoles, extracted from the calculation performed by BigDFT
    """
    _filename = 'transition_quantities.txt'

    def subportion(self, norb_occ, norb_virt):
        """Extract a subportion of the Transition Multipoles.

        Returns a set of transition multipoles which is made by only
        considering the first ``norb_occ`` and ``norb_virt`` orbitals

        Args:
           norb_occ (tuple): new number of occupied orbitals. Must be lower
             that the instance value
           norb_virt (tuple): new number of virtual orbitals. Must be lower
             that the instance value

        Returns:
           TransitionMultipoles: reduced transition multipoles
        """
        inds = self._subindices(norb_occ, norb_virt)
        mat = np.array(self[inds])
        return TransitionMultipoles(matrix=mat, norb_occ=norb_occ,
                                    norb_virt=norb_virt)

    def get_transitions(self):
        """
        Get the transition quantities as the dimensional objects which should
        contribute to the oscillator strengths.

        Returns:
            numpy.array: Transition quantities multiplied by the square root of
            the unperturbed transition energy
        """
        newdipole = []
        for line in self:
            newdipole.append(np.ravel(line[0, 0]*line[0, 1:]))
        return np.array(newdipole)


class TransitionDipoles(TransitionMultipoles):
    """
    Transition dipoles as provided in the version of the code < 1.8.0.
    Deprecated, to be used in some particular cases
    """
    _filename = 'transition_dipoles.txt'

    def get_transitions(self):
        return self


class Excitations():
    """LR Excited states of a system

    Definition of the excited states in the Casida Formalism

    Args:
       cm (CouplingMatrix): the matrix of coupling among transitions
       tm (TransitionMultipoles): scalar product of multipoles among
         transitions

    """

    def __init__(self, cm, tm):
        self.cm = cm
        self.tm = tm
        self.eigenvalues, self.eigenvectors = cm.diagonalize()

        # : array: transition quantities coming from the multipoles
        self.transitions = tm.get_transitions()

        scpr = np.array(np.dot(self.eigenvectors, self.transitions))

        #: array: oscillator strenghts components of the transitions defined
        #  as the square of  $\int w_a(\mathbf r) r_i $
        self.oscillator_strenghts = np.array([t**2 for t in scpr[:, 0:3]])

        # : array: average of all the components of the OS
        self.avg_os = np.average(self.oscillator_strenghts, axis=1)

        self.alpha_prime = 2.0*self.oscillator_strenghts / \
            self.eigenvalues[:, np.newaxis]
        """ array: elements of the integrand of the statical polarizability in
            the space of the excitations """

        self._indices_for_spin_comparison = \
            self._get_indices_for_spin_comparison()

        self.identify_singlet_and_triplets(1.e-5)

    def _get_indices_for_spin_comparison(self):
        inds = [[], []]
        inds0 = []
        # get the indices for comparison, take the minimum between the spins
        if len(self.cm.norb_occ) == 1:
            nocc = self.cm.norb_occ[0]
            nvirt = self.cm.norb_virt[0]
            nos = [nocc, nocc]
            nvs = [nvirt, nvirt]
        else:
            nocc = min(self.cm.norb_occ)
            nvirt = min(self.cm.norb_virt)
            nos = self.cm.norb_occ
            nvs = self.cm.norb_virt
        for ispin in [0, 1]:
            for a in range(nvirt):
                for p in range(nocc):
                    inds[ispin].append([p, a, ispin])
            for a in range(nvirt, nvs[ispin]):
                for p in range(nocc, nos[ispin]):
                    inds0.append([p, a, ispin])
        transA = transition_indexes(
            self.cm.norb_occ, self.cm.norb_virt, inds[0])
        transB = transition_indexes(
            self.cm.norb_occ, self.cm.norb_virt, inds[1])
        trans0 = transition_indexes(self.cm.norb_occ, self.cm.norb_virt, inds0)
        return transA, transB, trans0

    def spectrum_curves(self, omega, slice=None, weights=None):
        """Calculate spectrum curves.

        Provide the set of the curves associated to the weights. The resulting
        curves might then be used to draw the excitation spectra.

        Args:
            omega (array): the linspace used for the plotting, of shape
              ``(n,)``. Must be provided in Atomic Units
            slice (array): the lookup array that has to be considered. if Not
              provided the entire range is assumed
            weights (array): the set of arrays used to weight the spectra. Must
              have shape ``(rank,m)``, where ``rank`` is equal to the number
              of eigenvalues. If m is absent it is assumed to be 1. When not
              specified, it defaults to the average oscillator strenghts.

        Returns:
            array: a set of spectrum curves, of shape equal to ``(n,m)``,
            where ``n`` is the shape of ``omega`` and ``m`` the size of the
            second dimension of ``weights``.
        """
        if slice is None:
            oo = self.eigenvalues[:, np.newaxis] - omega**2
            wgts = weights if weights is not None else self.avg_os
        else:
            oo = self.eigenvalues[slice, np.newaxis] - omega**2
            oo = oo[0]
            wgts = weights if weights is not None else self.avg_os[slice]
        return np.dot(2.0/oo.T, wgts)

    def identify_singlet_and_triplets(self, tol=1.e-5):
        """
        Find the lookup tables that select the singlets and the triplets among
          the excitations

        Args:
           tol (float): tolerance to be applied to recognise the spin character
        """
        sings = []
        trips = []
        for exc in range(len(self.eigenvalues)):
            sing, trip = self.project_on_spin(exc, tol)
            if sing:
                sings.append(exc)
            if trip:
                trips.append(exc)

        if len(sings) > 0:
            self.singlets = (np.array(sings),)
            """array:  lookup table of the singlet excitations"""

        if len(trips) > 0:
            self.triplets = (np.array(trips),)
            """array:  lookup table of the triplet excitations"""

    def _project_on_occ(self, exc):
        """
        Project a given eigenvector on the occupied orbitals.
        In the spin averaged case consider all the spin indices nonetheless
        """
        norb_occ = self.cm.norb_occ
        norb_virt = self.cm.norb_virt
        pProj_spin = []
        for ispin, norb in enumerate(norb_occ):
            pProj = np.zeros(norb)
            for iorb in range(norb):
                harvest = _collection_indexes_iocc(
                    iorb, self.cm.norb_virt, spin=None if len(norb_occ) == 1
                    else ispin)
                inds = np.array(transition_indexes(
                    norb_occ, norb_virt, harvest))
                pProj[iorb] = np.sum(np.ravel(self.eigenvectors[exc, inds])**2)
            pProj_spin.append(pProj)
        return pProj_spin

    def project_on_spin(self, exc, tol=1.e-8):
        """
        Control if an excitation has a Singlet or Triplet character

        Args:
            exc (int): index of the excitation to be controlled

        Returns:
            tuple (bool,bool): couple of booleans indicating if the excitation
            is a singlet or a triplet respectively
        """
        A, B, zero = [np.ravel(self.eigenvectors[exc, ind])
                      for ind in self._indices_for_spin_comparison]
        issinglet = np.linalg.norm(A-B) < tol
        istriplet = np.linalg.norm(A+B) < tol
        return issinglet, istriplet
        # print (self.eigenvalues[exc], np.linalg.norm(A), np.linalg.norm(B),
        #        A-B,A+B, np.linalg.norm(zero))

    def _get_threshold(self, pProj_spin, th_energies, tol):
        """
        Identify the energy which is associated to the threshold of a given
        excitation. The tolerance is used to discriminate the component
        """
        ths = -1.e100
        for proj, en in zip(pProj_spin, th_energies):
            norb = len(en)
            pProj = proj.tolist()
            pProj.reverse()
            imax = norb-1
            for val in pProj:
                if val > tol:
                    break
                imax -= 1
            ths = max(ths, en[imax])
        return ths

    def split_excitations(self, evals, tol, nexc='all'):
        """Separate the excitations in channels.

        This methods classify the excitations according to the channel they
        belong, and determine if a given excitation might be considered as a
        belonging to a discrete part of the spectrum or not.

        Args:
            evals (BandArray): the eigenvalues as they are provided
              (for instance) from a `Logfile` class instance.
            tol (float): tolerance for determining the threshold
            nexc (int,str): number of excitations to be analyzed.
                If ``'all'`` then the entire set of excitations are analyzed.

        """
        self.determine_occ_energies(evals)
        self.identify_thresholds(self.occ_energies, tol, len(
            self.eigenvalues) if nexc == 'all' else nexc)

    def identify_thresholds(self, occ_energies, tol, nexc):
        """Identify the thresholds per excitation.

        For each of the first ``nexc`` excitations, identify the energy value
        of its corresponding threshold. This value is determined by projecting
        the excitation components on the occupied states and verifying that
        their norm for the highest energy level is below a given tolerance.

        Args:
           occ_energies (tuple of array-like): contains the list of the
             energies of the occupied states per spin channel
           tol (float): tolerance for determining the threshold
           nexc (int): number of excitations to be analyzed
        """
        # : Norm of the $w_p^a$ states associated to each excitation
        self.wp_norms = []
        threshold_energies = []
        for exc in range(nexc):
            proj = self._project_on_occ(exc)
            self.wp_norms.append(proj)
            threshold_energies.append(
                self._get_threshold(proj, occ_energies, tol))
        # : list: identified threshold for inspected excitations
        self.threshold_energies = np.array(threshold_energies)

        self.excitations_below_threshold = np.where(
            np.abs(self.threshold_energies) >= np.sqrt(
                self.eigenvalues[0:len(self.threshold_energies)]))
        """ array: Indices of the excitations which lie below their
                    corresponding threshold """

        self.excitations_above_threshold = np.where(
            np.abs(self.threshold_energies) <
            np.sqrt(self.eigenvalues[0:len(self.threshold_energies)]))
        """ array: Indices of the excitations which lie above their
                   corresponding threshold """

    def determine_occ_energies(self, evals):
        """
        Extract the occupied energy levels from a Logfile BandArray structure,
        provided the tuple of the number of occupied states

        Args:
            evals (BandArray): the eigenvalues as they are provided
              (for instance) from a `Logfile` class instance.
        """
        norb_occ = self.cm.norb_occ
        occ_energies = []
        # istart=0
        for ispin, norb in enumerate(norb_occ):  # range(len(norb_occ)):
            # istart:istart+norb_occ[ispin]]))
            occ_energies.append(np.array(evals[0][ispin][0:norb]))
            # istart+=norb_tot[ispin]
        # : array: energies of the occupied states out of the logfile
        self.occ_energies = occ_energies
        # : float: lowest threshold of the excitations. All excitations are
        #   discrete below this level
        self.first_threshold = abs(
            max(np.max(self.occ_energies[0]), np.max(self.occ_energies[-1])))

    def plot_alpha(self, **kwargs):
        """Plot the imaginary part.

        Plot the real or imaginary part of the dynamical polarizability.

        Keyword Arguments:
           real (bool): True if real part has to be plotted. The imaginary
             part is plotted otherwise
           eta (float): Value of the complex imaginary part. Defaults to 1.e-2.
           group (str): see :meth:`lookup`

           **kwargs:
               other arguments that might be passed to the :meth:`plot` method
               of the :mod:`matplotlib.pyplot` module.

        Returns:
             :mod:`matplotlib.pyplot`: the reference to
             :mod:`matplotlib.pyplot` module.
        """
        import matplotlib.pyplot as plt
        from futile.Utils import kw_pop
        emax = np.max(np.sqrt(self.eigenvalues))*HaeV
        kwargs, real = kw_pop('real', False, **kwargs)
        plt.xlim(xmax=emax)
        if real:
            plt.ylabel(r'$\mathrm{Re} \alpha$ (AU)', size=14)
        else:
            plt.ylabel(r'$\mathrm{Im} \alpha$', size=14)
            plt.yticks([])
        plt.xlabel(r'$\omega$ (eV)', size=14)
        if hasattr(self, 'first_threshold'):
            eps_h = self.first_threshold*HaeV
            plt.axvline(x=eps_h, color='black', linestyle='--')
        kwargs, eta = kw_pop('eta', 1.e-2, **kwargs)
        omega = np.linspace(0.0, emax, 5000)+2.0*eta*1j
        kwargs, group = kw_pop('group', 'all', **kwargs)
        slice = self.lookup(group)
        spectrum = self.spectrum_curves(omega, slice=slice)
        toplt = spectrum.real if real else spectrum.imag
        pltkwargs = dict(c='black', linewidth=1.5)
        pltkwargs.update(kwargs)
        plt.plot(omega*HaeV, toplt, **pltkwargs)
        return plt

    def lookup(self, group):
        """
        Identify the group of the excitations according to the argument

        Args:
            group (str):
               A string chosen between

               * ``"all"`` : provides the entire set of excitations
                 (:py:class:`None` instead of the lookup array)

               * ``"bt"`` : provides only the excitations below threshold

               * ``"at"`` : provides only the excitations above threshold

               * ``"singlets"`` : provides the index of the excitations that
                 have a singlet character

               * ``"triplets"`` : provides the index of the excitations that
                 have a triplet character

        """
        slice = None
        if group == 'bt':
            slice = self.excitations_below_threshold
        if group == 'at':
            slice = self.excitations_above_threshold
        if group == 'singlets':
            slice = self.singlets
        if group == 'triplets':
            slice = self.triplets
        return slice

    def plot_excitation_landscape(self, **kwargs):
        """
        Represent the excitation landscape as splitted in the excitations class

        Args:
            **kwargs:
               keyword arguments to be passed to the `pyplot` instance.
               The ``xlabel``, ``ylabel`` as well as ``xlim`` are already set.

        Returns:
            :mod:`matplotlib.pyplot`: the reference to :mod:`matplotlib.pyplot`
              module.

        Example:
           >>> ex=Excitations(cm,tm)
           >>> ex.split_excitations(evals=...,tol=1.e-4,nexc=...)
           >>> ex.plot_excitation_landscape(title='Excitation landscape')
        """
        import matplotlib.pyplot as plt
        Emin = 0.0
        Emax = np.max(np.sqrt(self.eigenvalues))*HaeV
        for level in self.occ_energies[0]:
            eng_th = level*HaeV
            plt.plot((Emin, eng_th), (level, level),
                     '--', c='red', linewidth=1)
            plt.plot((eng_th, Emax), (level, level), '-', c='red', linewidth=1)
            plt.scatter(abs(eng_th), level, marker='x', c='red')
        ind_bt = self.excitations_below_threshold
        exc_bt = np.sqrt(self.eigenvalues)[ind_bt]
        lev_bt = self.threshold_energies[ind_bt]
        plt.scatter(HaeV*exc_bt, lev_bt, s=16, marker='o', c='black')
        ind_at = self.excitations_above_threshold
        exc_at = np.sqrt(self.eigenvalues)[ind_at]
        lev_at = self.threshold_energies[ind_at]
        plt.scatter(HaeV*exc_at, lev_at, s=14, marker='s', c='blue')
        plt.xlabel('energy (eV)')
        plt.ylabel('Threshold energy (Ha)')
        plt.xlim(xmin=Emin-1, xmax=Emax)
        for attr, val in kwargs.items():
            if type(val) == dict:
                getattr(plt, attr)(**val)
            else:
                getattr(plt, attr)(val)
        return plt

    def dos_dict(self, group='all'):
        """Dictionary for DoS creation.

        Creates the keyword arguments that have to be passed to the
        `meth:BigDFT.DoS.append` method of the `DoS` class

        Args:
           group (str): see :meth:`lookup`

        Returns:
            :py:class:`dict`: kewyord arguments that can be passed to the
            `meth:BigDFT.DoS.append` method of the :class:`DoS.DoS` class

        """
        ev = np.sqrt(self.eigenvalues)
        slice = self.lookup(group)
        if slice is not None:
            ev = ev[slice]
        return dict(energies=np.array([np.ravel(ev)]), units='AU')

    def dos(self, group='all', **kwargs):
        """Density of States of the Excitations.

        Provides an instance of the :class:`~BigDFT.DoS.DoS` class,
        corresponding to the Excitations instance.

        Args:
           group (str): see :meth:`lookup`

           **kwargs: other arguments that might be passed to the
             :class:`DoS.DoS` instantiation

        Returns:
            :class:`DoS.DoS`: instance of the Density of States class
        """
        from BigDFT.DoS import DoS
        kwa = self.dos_dict(group=group)
        kwa['energies'] = kwa['energies'][0]
        if hasattr(self, 'first_threshold'):
            kwa['fermi_level'] = self.first_threshold
        else:
            kwa['fermi_level'] = 0.0
        kwa.update(kwargs)
        return DoS(**kwa)

    def plot_Sminustwo(self, coord, alpha_ref=None, group='all'):
        """Inspect S-2 sum rule.

        Provides an handle to the plotting of the $S^{-2}$ sum rule, which
        should provide reference values for the static polarizability tensor.

        Args:
            coord (str): the coordinate used for inspection. May be ``'x'``,
              ``'y'`` or ``'z'``.
            alpha_ref (list): diagonal of the reference static polarizability
               tensor (for instance calculated via finite differences).
               If present the repartition of the contribution of the various
               groups of excitations is plotted.
            group (str): see :meth:`lookup`

        Returns:
            reference to :mod:`matplotlib.pyplot` module.

        """
        import matplotlib.pyplot as plt
        idir = ['x', 'y', 'z'].index(coord)
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('energy (eV)', size=14)
        plt.ylabel(r'$\alpha_{'+coord+coord+r'}$ (AU)', size=14)
        if alpha_ref is not None:
            plt.axhline(y=alpha_ref[idir], color='r', linestyle='--')
        if hasattr(self, 'first_threshold'):
            eps_h = abs(HaeV*self.first_threshold)
            plt.axvline(x=eps_h, color='black', linestyle='--')
        e = np.sqrt(self.eigenvalues)*HaeV
        w_ii = self.alpha_prime[:, idir]
        slice = self.lookup(group)
        if slice is not None:
            e = e[slice]
            w_ii = w_ii[slice]
        ax1.plot(e, np.cumsum(w_ii))
        ax2 = ax1.twinx()
        ax2.plot(e, w_ii, color='grey', linestyle='-')
        plt.ylabel(r'$w_{'+coord+coord+r'}$ (AU)', size=14)

        return plt


def get_alpha_energy(log, norb, nalpha):
    return log.evals[0][0][norb+nalpha-1]


def identify_contributions(numOrb, na, exc, C_E2):
    pProj = np.zeros(numOrb*2)
    for p in range(numOrb):
        for spin in [0, 1]:
            # sum over all the virtual orbital and spin
            for alpha in range(na):
                # extract the value of the index of C_E2
                elements = transition_indexes(
                    [numOrb], [na], [[p, alpha, spin]])
                for el in elements:
                    pProj[p+numOrb*spin] += C_E2[exc][el]**2
    pProj = pProj[0:numOrb]+pProj[numOrb:2*numOrb]  # halves the components
    return pProj


def get_p_energy(log, norb):
    return log.evals[0][0][0:norb]


def get_threshold(pProj, th_energies, th_levels, tol):
    norb = len(th_energies)
    pProj = pProj.tolist()
    pProj.reverse()
    imax = norb-1
    for val in pProj:
        if val > tol:
            break
        imax -= 1
    return [th_levels[imax], th_energies[imax]]
