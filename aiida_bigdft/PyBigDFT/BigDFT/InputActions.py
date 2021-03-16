"""Actions to define on the Input parameters.

This module defines some of the most common actions that a BigDFT user might
like to perform on the input file. Such module therefore sets some of the keys
of the input dictionary to the values needed to perform the operations.
Users might also inspire to the actions performed in order to customize the
runs in a different way. All the functions of this module have as first
argument ``inp``, the dictionary of the input parameters.

Many other actions are available in BigDFT code. This module only regroups the
most common. Any of these functionalities might be removed from the input file
by the :py:func:`remove` function.

Note:

   Any of the action of this module, including the :py:func:`remove` function,
   can be also applied to an instance of the
   :py:class:`BigDFT.Inputfiles.Inputfile` class, by removing the first
   argument (``inp``). This adds extra flexibility as the same method may be
   used to a dictionary instance or to a BigDFT input files.
   See the example :ref:`input_action_example`.

We now list the available methods, in order of category.

Basic Setup and common functionalities
--------------------------------------------

.. autosummary::

   set_xc
   set_hgrid
   set_rmult
   set_mesh_sizes
   optimize_geometry
   spin_polarize
   charge
   charge_and_polarize
   set_symmetry
   apply_electric_field
   add_empty_scf_orbitals
   extract_virtual_states
   set_electronic_temperature
   set_implicit_solvent
   set_dispersion_correction

Self-Consistent-Field setup and algorithms
-------------------------------------------

.. autosummary::

   set_scf_convergence
   set_random_inputguess
   set_linear_scaling
   set_scf_method
   set_kernel_guess
   optimize_kernel
   optimize_coefficients
   optimize_support_functions

Input-Output and restart
------------------------

.. autosummary::

   write_orbitals_on_disk
   read_orbitals_from_disk
   write_density_on_disk
   change_data_directory
   connect_run_data
   write_support_function_matrices
   write_cubefiles_around_fermi_level
   write_support_functions_on_disk


Setting atomic-based information, other functionalities
-------------------------------------------------------

.. autosummary::

   set_atomic_positions
   calculate_dipole
   use_gpu_acceleration
   set_orbital_occupancy
   calculate_tddft_coupling_matrix
   set_external_potential
   set_psp_directory
   set_psp_file
   load
   remove

Note:

   Each of the actions here **must** have default value for the arguments
   (except the input dictionary ``inp``). This is needed for a good behaviour
   of the function `remove`.

We list here the extended documentation in alphabetic order.
"""

from futile.Utils import dict_set

__set__ = dict_set
"""func: Action function.

This is the pointer to the set function, useful to modify the action with the
undo method

"""


def __undo__(inp, *subfields):
    """
    Eliminate the last item of the subfields as provided to dict_set
    """
    from futile.Utils import push_path
    # remove the last key until the parent is empty
    lastkey = -1
    tmp = {}
    while len(subfields) > -lastkey and tmp == {}:
        keys = subfields[:lastkey]
        tmp, k = push_path(inp, *keys)
        tmp.pop(k)
        lastkey -= 1


def remove(inp, action):
    """Remove action from the input dictionary.

    Remove an action from the input file, thereby restoring the **default**
    value, as if the action were not specified.

    Args:
       inp (dict): dictionary to remove the action from.
       action (func): one of the actions of this module. It does not need to be
       specified before, in which case it produces no effect.

    Example:
       >>> from Calculators import SystemCalculator as C
       >>> code=C()
       >>> inp={}
       >>> set_xc(inp,'PBE')
       >>> write_orbitals_on_disk(inp)
       >>> log=code.run(input=inp) # perform calculations
       >>> remove(inp, write_orbitals_on_disk) #remove the action
       >>> read_orbitals_from_disk(inp)
       >>> # this will restart the scf from the previous orbitals
       >>> log2=code.run(input=inp)
    """
    global __set__
    __set__ = __undo__
    action(inp)
    __set__ = dict_set


def set_hgrid(inp, hgrids=0.4):
    """
    Set the wavelet grid spacing.

    Args:
       hgrid (float,list): list of the grid spacings in the three directions.
       It might also be a scalar, which implies the same spacing
    """
    __set__(inp, 'dft', 'hgrids', hgrids)


def set_scf_convergence(inp, gnrm='default', rpnrm='default'):
    """
    Set the tolerance acceptance level for stopping the self-consistent
    iterations. Useful both for LS and CS

    Args:
       gnrm (float): the tolerance level for the CS inner loop
       rpnrm (float): residue for the density/or potential norm. Useful both
           for CS and LS
    """
    __set__(inp, 'dft', 'gnrm_cv', gnrm)
    __set__(inp, 'lin_general', 'rpnrm_cv', rpnrm)
    __set__(inp, 'mix', 'rpnrm_cv', rpnrm)


def set_rmult(inp, rmult=None, coarse=5.0, fine=8.0):
    """
    Set the wavelet grid extension by modifying the multiplicative radii.

    Args:
       rmult (float,list): list of two values that have to be used for the
         coarse and the fine resolution grid. It may also be a scalar.
       coarse (float): if the argument ``rmult`` is not provided it sets the
         coarse radius multiplier
       fine (float): if the argument ``rmult`` is not provided it sets the fine
         radius multiplier
    """
    rmlt = [coarse, fine] if rmult is None else rmult
    __set__(inp, 'dft', 'rmult', rmlt)


def set_symmetry(inp, yes=True):
    """
    Set the symmetry detection for the charge density and the ionic forces and
    stressdef set_symmetry(inp,yes=True):

    Args:
       yes (bool): If ``False`` the symmetry detection is disabled
    """
    __set__(inp, 'dft', 'disablesym', not yes)


def set_linear_scaling(inp):
    """
    Activates the linear scaling mode
    """
    newid = 'linear'
    previous_ipid = inp.get('dft', 'False')
    if previous_ipid:
        previous_ipid = inp.get('inputpsiid', 'False')
    if previous_ipid == 2:
        newid = 102

    __set__(inp, 'dft', 'inputpsiid', newid)


def set_mesh_sizes(inp, ngrids=64):
    """
    Constrain the number of grid points in each direction.
    This is useful when performing periodic system calculations with variable
    cells which need to be compared each other. In this way the number of
    degrees of freedom is kept constant throughout the various simuilations.

    Args:
       ngrids (int,list): list of the number of mesh points in each direction.
         Might be a scalar.
    """
    __set__(inp, 'dft', 'ngrids', ngrids)


def spin_polarize(inp, mpol=1):
    """
    Add a collinear spin polarization to the system.

    Arguments:
       mpol (int): spin polarization in Bohr magneton units.
    """
    __set__(inp, 'dft', 'nspin', 2)
    __set__(inp, 'dft', 'mpol', mpol)


def charge(inp, charge=-1):
    """
    Charge the system

    Arguments:
        charge (int,float): value of the charge in units of *e* (the electron
          has charge -1). Also accept floating point numbers.
    """
    __set__(inp, 'dft', 'qcharge', charge)


def apply_electric_field(inp, elecfield=[0, 0, 1.e-3]):
    """
    Apply an external electric field on the system

    Args:
       electric (list, float): Values of the Electric Field in the three
         directions. Might also be a scalar.
    """
    __set__(inp, 'dft', 'elecfield', [e for e in elecfield])


def charge_and_polarize(inp):
    """
    Charge the system by removing one electron. Assume that the original
      system is closed shell, thus polarize.
    """
    charge(inp, charge=1)
    spin_polarize(inp, mpol=1)


def set_scf_method(inp, method='dirmin', mixing_on='density',
                   mixing_scheme='Pulay', ):
    """
    Set the algorithm for scf.

    Args:
       method (str): The algoritm chosen. Might be different for the cubic (CS)
         or linear scaling (LS) algorithm.
         * dirmin: Direct minimization approach (only CS)
         * mixing: Mixing scheme (only CS)
         * hybrid: outer loop of the LS mode
         * two_levels: two level of accuracy in the scf

       mixing_on (str): May be ``"density"`` or ``"potential"`` in the
         ``"mixing"`` case, decide to which quantity the mixing to be performed

       mixing_scheme (str): May be:

          * Pulay : DIIS mixing on the last 7 iterations

          * Simple: Simple mixing

          * Anderson: Anderson scheme

          * Anderson2: Anderson scheme based on the two pervious iterations

          * CG: Conjugate Gradient based on the minimum of the energy with
                respect of the potential

    Warning:
       Only the FOE method exhibit asymptotic linear scaling regime.

    Todo:
       Check if the linear scaling case needs another input variable for the
       mixing of the potential (density)

    """
    method.upper()
    if method == 'DIRMIN':
        __set__(inp, 'mix', 'iscf', 0)
        return
    if method == 'MIXING':
        iscf = 0
        if mixing_on == 'density':
            iscf += 10
        if mixing_scheme == 'Pulay':
            iscf += 7
        if mixing_scheme == 'Anderson':
            iscf += 3
        if mixing_scheme == 'Anderson2':
            iscf += 4
        if mixing_scheme == 'Simple':
            iscf += 2
        if mixing_scheme == 'CG':
            iscf += 5
        __set__(inp, 'mix', 'iscf', iscf)
        return
    if method == 'HYBRID':
        __set__(inp, 'lin_general', 'hybrid', True)
    if method == 'TWO_LEVELS':
        __set__(inp, 'lin_general', 'hybrid', False)


def add_empty_scf_orbitals(inp, norbs=10):
    """
    Insert ``norbs`` empty orbitals in the scf procedure

    Args:
       norbs (int): Number of empty orbitals

    Warning:
       In linear scaling case, this is only meaningful for the direct
       minimization approach.
    """
    __set__(inp, 'mix', 'norbsempty', norbs)
    __set__(inp, 'lin_general', 'extra_states', norbs)


def write_cubefiles_around_fermi_level(inp, nplot=1):
    """
    Writes the ``nplot`` orbitals around the fermi level in cube format

    Args:
      nplot (int): the number of orbitals to print around the fermi level.

    Warning:
       This is presently meaningful only for a empty states calculation.

    Warning:
       This would work only for the cubic scaling code at present.
    """
    __set__(inp, 'dft', 'nplot', nplot)


def write_orbitals_on_disk(inp, format='binary'):
    """
    Set the code to write the orbitals on disk in the provided format

    Args:
      format (str): The format to write the orbitals with. Accepts the strings:
         * 'binary'
         * 'text'
         * 'etsf' (requires etsf-io enabled)

    Todo:
      Verify if this option works for a linear scaling calulation.
    """
    fmt = format
    __set__(inp, 'output', 'orbitals', fmt)


def write_support_functions_on_disk(inp, format='text'):
    """Dump the support functions.

    Write the support function which are expressed in wavelet basis in the
    code as files at the end of the calculation.

    Args:
        format (str): the format of the data, can be 'text', 'ETSF' or 'binary'
             or one of the allowed codes of the `output_wf` key.
    """
    __set__(inp, 'lin_general', 'output_wf', format)


def write_support_function_matrices(inp, format='text'):
    """
    Write the matrices of the linear scaling formats.

    Args:
       format (str): The format to write the orbitals with. Accepts the
         strings:
          * 'binary'
          * 'text'

    Todo:
       Verify if the binary format is available and set the appropriate values
    """
    fmt = 0
    if format == 'text':
        fmt = 1
    elif format == 'binary':
        fmt = 4
    __set__(inp, 'lin_general', 'output_mat', fmt)


def set_atomic_positions(inp, posinp=None):
    """
    Insert the atomic positions as a part of the input dictionary
    """
    __set__(inp, 'posinp', posinp)


def read_orbitals_from_disk(inp, mode='cubic'):
    """
    Read the orbitals from data directory, if available.

    Args:
        mode (str): can be 'cubic' or 'linear'. In the first case,
            orbitals are read as KS objects, whereas in the
            second kernel (or coeffs) and support functions are
            expressed.
    """
    newid = 2
    previous_ipid = inp.get('dft', 'False')
    if previous_ipid:
        previous_ipid = inp.get('inputpsiid', 'False')
    if previous_ipid == 'linear' or previous_ipid == 100:
        newid = 102
    if mode == 'linear':
        newid = 102
    __set__(inp, 'dft', 'inputpsiid', newid)


def set_kernel_guess(inp, mode='kernel'):
    """ Method for guessing the kernel at restart.

    Args:
        mode (str): Guessing method. Can be:
           * kernel
           * coefficients
           * random
           * diag
           * weight
           * ao
           * full_kernel
           * full_coefficients
    """
    __set__(inp, 'lin_general', 'kernel_restart_mode', mode)


def set_random_inputguess(inp):
    """
    Input orbitals are initialized as random coefficients
    """
    __set__(inp, 'dft', 'inputpsiid', -2)


def set_electronic_temperature(inp, kT=1.e-3, T=0):
    """
    Define the electronic temperature, in AU (``kT``) or K (``T``)
    """
    TtokT = 8.617343e-5/27.21138505
    tel = TtokT*T if T != 0 else kT
    __set__(inp, 'mix', 'tel', tel)


def optimize_geometry(inp, method='FIRE', nsteps=50, betax=4.0, frac_fluct=1.0,
                      forcemax=0):
    """
    Optimize the geometry of the system

    Args:
       nsteps (int): maximum number of atomic steps.
       method (str): Geometry optimizer. Available keys:

          * SDCG:   A combination of Steepest Descent and Conjugate Gradient

          * VSSD:   Variable Stepsize Steepest Descent method

          * LBFGS:  Limited-memory BFGS

          * BFGS:   Broyden-Fletcher-Goldfarb-Shanno

          * PBFGS:  Same as BFGS with an initial Hessian from a force field

          * DIIS:   Direct inversion of iterative subspace

          * FIRE:   Fast Inertial Relaxation Engine, described by Bitzek et al.

          * SQNM:   Stabilized quasi-Newton minimzer
       betax (float): the step size for the optimization method.
          This stepsize is system dependent and it has therefore to be
          determined for each system.
       frac_fluct (float): Fraction of force fluctuations. Stop if
          fmax < forces_fluct*frac_fluct.
       forcemax (float): Max forces criterion when stop.
    """
    __set__(inp, 'geopt', 'method', method)
    __set__(inp, 'geopt', 'ncount_cluster_x', nsteps)
    __set__(inp, 'geopt', 'betax', betax)
    __set__(inp, 'geopt', 'frac_fluct', frac_fluct)
    __set__(inp, 'geopt', 'forcemax', forcemax)


def set_xc(inp, xc='PBE'):
    """
    Set the exchange and correlation approximation

    Args:
       xc (str): the Acronym of the XC approximation

    Todo:
       Insert the XC codes corresponding to ``libXC`` conventions
    """
    __set__(inp, 'dft', 'ixc', xc)


def write_density_on_disk(inp):
    """
    Write the charge density on the disk after the last scf convergence
    """
    __set__(inp, 'dft', 'output_denspot', 21)


def use_gpu_acceleration(inp):
    """
    Employ gpu acceleration when available, for convolutions and Fock operator

    Todo:
       Verify what happens when only one of the functionality is enabled at
       compile-time
    """
    __set__(inp, 'perf', 'accel', 'OCLGPU')
    __set__(inp, 'psolver', 'setup', 'accel', 'CUDA')


def set_scf_iterations(inp, nit=[50, 1]):
    """
    Set the number of the iteration per loop

    Args:
       nit (int,list): integer of the number of iterations. Might be a scalar
            or a list, up to length two. The first element of the list contains
            the number of iterations of the direct minimization loop. if
            ``nit`` is a scalar, only this contribution is taken into account.
            The second element is the number of subspace iterations where the
            hamiltonian is diagonalized in the subspace. For a LS calculation,
            the number of iteration correspond to the levels of the outer loop.
    """
    try:
        nlen = len(nit)
    except TypeError:
        nlen = 0
    if nlen >= 1:
        __set__(inp, 'dft', 'itermax', nit[0])
    if nlen == 2:
        __set__(inp, 'dft', 'nrepmax', nit[1])
    if nlen == 0:
        __set__(inp, 'dft', 'itermax', nit)
    __set__(inp, 'lin_general', 'nit', nit)


def change_data_directory(inp, name=''):
    """
    Modify the name of the ``data-`` directory.
    Useful to grab the orbitals from another directory than the run name

    Args:
       name (str): the name of the run
    """
    __set__(inp, 'radical', name)


def calculate_tddft_coupling_matrix(inp, tda=False, rpa=True, fxc=True):
    """
    Perform a Casida TDDFT coupling matrix extraction.

    Args:
       tda (bool): when ``True``, Tamm-Dancoff approximation is used for the
         extraction of the coupling matrix
       rpa (bool): when ``False``, the calculation of the RPA term (the linear
         response of the hartree potential) is switched off
       fxc (bool): when ``False``, the calculation of the fxc term (the linear
         response of the XC operator) is switched off.

    Note:
       The arguments ``fxc`` and ``rpa`` should not be simultaneously
       ``False``.

    Warning:
       Presently the LR-TDDFT casida fxc is only available for LDA
       functionals in ABINIT flavour.
    """
    approach = 'TDA' if tda else 'full'
    __set__(inp, 'tddft', 'tddft_approach', approach)
    if rpa and fxc:
        output = 'complete'
    elif rpa:
        output = 'rpa'
    else:
        output = 'fxc'
    __set__(inp, 'output', 'coupling_matrix', output)


def extract_virtual_states(inp, nvirt=8, davidson=False, norbv=None,
                           itermax_virt=150):
    """
    Extract a given number of empty states **after** the scf cycle.

    Args:
       davidson (bool): If set to ``True`` activates davidson calculation,
           otherwise Trace Minimization of the Hamiltonian is employed.
       norbv (int): Defines the total size of the virtual subspace,
           which may be larger than nvirt.
    """
    if norbv is None:
        norbv = nvirt
    nv = norbv if davidson else -norbv
    __set__(inp, 'dft', 'norbv', nv)
    __set__(inp, 'dft', 'nvirt', nvirt)
    __set__(inp, 'dft', 'itermax_virt', itermax_virt)


def connect_run_data(inp, log=None):
    """
    Associate the data of the run of a given logfile to the input
    by retrieving the data directory name of the logfile.

    Args:
       log (Logfile): instance of a Logfile class

    """
    if log is None:
        change_data_directory(inp)  # no effect
    else:
        ll = log if len(log) == 0 else log[0]
        change_data_directory(inp, ll.log['radical'])


def calculate_dipole(inp):
    """
    Extract the dipole moment from the total charge density.

    Note:
      This function is useful for the linear scaling setup as the cubic
      scaling approach always calculates the charge density multipoles.
    """
    __set__(inp, 'lin_general', 'calc_dipole', True)


def set_external_potential(inp, mm_pot):
    """
    Set the external potential to which the system is submitted outside the QM
    region

    Args:
       mm_pot (dict): dictionary of the external potential which contains the
          information on the counter-ions
    """
    __set__(inp, 'dft', 'external_potential', mm_pot)


def set_implicit_solvent(inp, itermax=20, minres=0.0001):
    """
    Add an acqueous implicit solvent around the system with the
        soft-sphere cavity method

    Args:
       itermax(int): maximum number of iteration of the
           Generalized Poisson Solver
       minres(float): minimum residue of the CG method to achieve
           GPS convergence
    """
    __set__(inp, 'psolver', 'environment', 'cavity', 'soft-sphere')
    __set__(inp, 'psolver', 'environment', 'itermax', itermax)
    __set__(inp, 'psolver', 'environment', 'minres', minres)


def set_dispersion_correction(inp):
    """
    Add Grimme's D3 correction to the energy and the forces

    Warning:
        Only works with Free Boundary conditions
    """
    __set__(inp, 'dft', 'dispersion', 5)


def set_psp_file(inp, filename=None, element=None):
    """
    Employ the given PSP file for the provided element

    Args:
       filename (str): the path of the psp file
       element(str): the element symbol for the PSP.
           Employs the psp name if absent
    """
    if element is None:
        element = filename.split('.')[-1]
    __set__(inp, 'psppar.'+element, filename)


def set_psp_directory(inp, directory='.'):
    """
    Employs all the "psppar.*" files of the given directory as pseudopotentials

    Args:
       directory (str): path of the psppar directory
    """
    from futile.Utils import file_list
    for psp in file_list(directory, prefix='psppar',
                         include_directory_path=True):
        set_psp_file(inp, filename=psp)


def load(inp, profile='', profiles=[]):
    """Load a profile or list of profiles.

    Args:
       profile (str): profile to be loaded, if a single profile is employed.
       profiles (list): profiles list to be loaded in subsequent order.
    """
    if len(profile) > 0:
        prof = profile
    else:
        prof = profiles
    __set__(inp, 'import', prof)


def optimize_kernel(inp, method='DIAG', dynamic_convergence=True, nit=5,
                    rpnrm=1.e-10, alphamix=0.5):
    """Methods for scf of the density kernel.

    Strategies for the optimization of the kernel.

    Args:
        dynamic_convergence (bool): if False, the threshold for the
             convergence of the kernel are not dynamically adjusted
        nit (int): number of scf iterations
        rpnrm (float): convergence criterion, change in density/potential
        method (str): optimization method ('DIRMIN', 'DIAG', 'NTPOLY', 'FOE')
        alphamix (float): mixing parameter
    """
    if not dynamic_convergence:
        __set__(inp, 'perf', 'adjust_kernel_iterations', False)
        __set__(inp, 'perf', 'adjust_kernel_threshold', False)
        __set__(inp, 'lin_kernel', 'delta_prnm', -1)
    __set__(inp, 'lin_kernel', 'nit', nit)
    __set__(inp, 'lin_kernel', 'rpnrm_cv', rpnrm)
    __set__(inp, 'lin_kernel', 'linear_method', method)
    __set__(inp, 'lin_kernel', 'alphamix', alphamix)


def optimize_coefficients(inp, nit=1, gnrm=1.e-5):
    """Methods for scf of the kernel coefficients.

    Args:
        nit (int): number of scf iterations
        gnrm (float): convergence criterion on the residue
    """
    __set__(inp, 'lin_kernel', 'gnrm_cv_coeff', gnrm)
    __set__(inp, 'lin_kernel', 'nstep', nit)


def optimize_support_functions(inp, nit=1, gnrm=1.e-6, diis_history=0):
    """Methods for scf of the Support functions.

    Args:
        nit (int): number of scf iterations
        gnrm (float): convergence criterion on the residue
    """
    __set__(inp, 'lin_basis', 'gnrm_cv', gnrm)
    __set__(inp, 'lin_basis', 'nit', nit)
    __set__(inp, 'lin_basis', 'idsx', diis_history)


def set_orbital_occupancy(inp, avg={}, up={}, down={}, kpoint=1):
    """Control the occupation number of the KS orbitals.

    With this funtionality one can fix the value of the occupation numbers of
    the KS orbitals.

    Args:
       avg (dict): dictionary of the form {i: f} where i is an integer stating
           the non-default value of the occupation number, for spin averaged
           occupancy
       up (dict): same as `avg` but for the up spin channel
       down (dict): same as `avg` but for the down spin channel
       kpoint (int): label of the kpoint associated to the occupancy
    """
    def orb_dict(dict):
        return {'Orbital '+str(i): f for i, f in dict.items()}
    spin_dict = {}
    if len(up) > 0:
        spin_dict['up'] = orb_dict(up)
    if len(down) > 0:
        spin_dict['down'] = orb_dict(down)
    if len(avg) > 0:
        spin_dict = orb_dict(avg)
    occupation_dict = {'K point '+str(kpoint): spin_dict}
    __set__(inp, 'occupation', str(occupation_dict))
