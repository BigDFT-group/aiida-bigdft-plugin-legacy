"""
This module defines some classes to perform a calculation using BigDFT
using binding (GIBinding) or using system call (SystemCalculator).

"""

# In our case for the class SystemCalculator which uses system calls:
# * We define posinp (equivalent of Atoms)
# * We have a python dictionary for the parameter
# * We define a calculator (equivalent of BFGS which is an Optimizer
#   (a method to optimize))
# Then we perform the method run.
#
# For the class GIBinding using the Gobject Introspection bindings, two methods
# set and update are added.
#
# The goal is to have a light Calculator almost compatible with ASE
# (Atomic Simulation environment, see https://gitlab.com/ase/ase)
# .. todo::
# In a future we add our method to ASE which is at a higher level
# (workflow of simulations).
# :Example:
# >>> from ase import Atoms
# >>> from ase.optimize import BFGS
# >>>  from ase.calculators.nwchem import NWChem
# >>>  from ase.io import write
# >>>  h2 = Atoms('H2',
# >>>             positions=[[0, 0, 0],
# >>>                        [0, 0, 0.7]])
# >>>  h2.calc = NWChem(xc='PBE')
# >>>  opt = BFGS(h2, trajectory='h2.traj')
# >>>  opt.run(fmax=0.02)
# >>>  BFGS:   0  19:10:49    -31.435229     2.2691
# >>>  BFGS:   1  19:10:50    -31.490773     0.3740
# >>>  BFGS:   2  19:10:50    -31.492791     0.0630
# >>>  BFGS:   3  19:10:51    -31.492848     0.0023
# >>>  write('H2.xyz', h2)
# >>>  h2.get_potential_energy()  # ASE's units are eV and Ang
# >>>  -31.492847800329216
##

import os
from futile.Utils import write as safe_print
import BigDFT.Logfiles as Lf


class GIBinding():
    """
    Calculator for BigDFT from Gobject Introspection bindings.
    """

    def __init__(self):
        # Import bindings about BigDFT (if the bindings are not generated, do
        # not work at all)
        from gi.repository import BigDFT
        self.runObj = -1
        # MPI initialisation
        (ierr, self.iproc, self.nproc, igroup, ngroup) = BigDFT.lib_init(0)
        self.runObj = None

    def update(self, inputfile):
        # If the inputpsiid is not present in the inputfile
        # assumes that the user wants to do a restart
        from futile.Utils import dict_merge
        if "dft" in inputfile and "inputpsiid" in inputfile["dft"]:
            var = inputfile
        else:
            var = inputfile.copy()
            dict_merge(var, {'dft': {'inputpsiid': 1}})
        from gi.repository import BigDFT
        self.runObj.update(BigDFT.Dict(var))

    def run(self):
        self.out = self.runObj.calculate(self.iproc, self.nproc)
        return self.out

    def set(self, inputfile=None):
        from gi.repository import BigDFT
        if inputfile is None:
            var = {}
        else:
            var = inputfile
        # Free memory first
        self.out = None
        self.runObj = None
        self.runObj = BigDFT.Run.new_from_dict(BigDFT.Dict(var))

    def __del__(self):
        if self.runObj == -1:
            return
        # MPI finalisation.
        self.out = None
        self.runObj = None
        from gi.repository import BigDFT
        BigDFT.lib_finalize()


class Runner():
    """Run of something.

    This object is associated with the concept of execution of a action.
    It may be customized to be used inside workflows and datasets.
    The central functionality is the `run` method that can be customized on
    subclasses of `Runner`. In this object there are global and local options
    of a run method. All arguments passed at the instantiation are stored as
    global options. For each call to `run`, these global options may updated by
    the arguments of the run call.

    Args:
        **kwargs: global options of the runner. Deepcopied in the dictionary
          returned by :meth:`global_options`.

    Example:

        >>> torun=Runner(args1='one',args2='two')
        >>> print(torun.global_options())
        {'args1':'one','args2':'two'}
        >>> print(torun.get_global_option('args1'))
        'one'

    """

    def __init__(self, **kwargs):
        import copy
        self._global_options = copy.deepcopy(kwargs)

    def global_options(self):
        """
        Get all global options dict.

        Returns:
            :py:class:`dict`: The dictionary of the global options in its
            current status
        """
        return self._global_options

    def get_global_option(self, key):
        """

        Get one key in global options

        Args:
           key (string): the global option key

        Returns:
            The value of the global options labelled by ``key``

        """
        return self._global_options[key]

    def update_global_options(self, **kwargs):
        """
        Update the global options by providing keyword arguments.

        Args:
           **kwargs: arguments to be updated in the global options
        """
        self._global_options.update(kwargs)

    def pop_global_option(self, key):
        """
        Remove a given global option from the global option dictionary

        Args:
           key (string): the global option key

        Returns:
           The value of the global option
        """
        self._global_options.pop(key)

    def _run_options(self, **kwargs):
        """
        Create a local dictionary for a specific run.
        It combines the present status of global option with the local
        dictionary of the run
        """
        import copy
        # First deepcopy from global_options and update from kwargs (warning: a
        # dictionary is not update)
        self.run_options = copy.deepcopy(self._global_options)
        """ :py:class`dict`: Local options of process_run.
         This dictionary can be accessed during the definition of the
         process_run method.
        It contains all the relevant keys for the definition of the runner.
        """

        self.run_options.update(kwargs)

    def run(self, **kwargs):
        """
        Run method of the class. It performs the following actions:

         * Constructs the local dictionary to be passed as ``**kwargs`` to the
           `process_run` function;
         * Calls the :meth:`pre_processing` method (intended to prepare some
           actions associated to the :meth:`process_run` method);
         * Calls :meth:`process_run` with the dictionary returned by
           :meth:`pre_processing` as  `**kwargs`;
         * Update such dictionary with the results returned by
           :meth:`process_run` and call :meth:`post_processing`;
         * Returns the object passed by the call to :meth:`post_processing`
           class method

        Developers are therefore expected to override :meth:`pre_processing`
        :meth:`process_run` and :meth:`post_processing`,
        when subclassing :class:`Runner`.

        """
        from futile.Utils import dict_merge
        self._run_options(**kwargs)
        run_args = self.pre_processing()
        run_results = self.process_run(**run_args)
        # safe_print('run_args',run_args,'run_results',run_results)
        dict_merge(dest=run_args, src=run_results)
        # safe_print('run_updated, again',run_args)
        return self.post_processing(**run_args)

    def pre_processing(self):
        """
        Pre-treat the keyword arguments and the options, if needed.

        Returns:
           :py:class:`dict`: dictionary of the pre-treated keyword arguments
           that have to be actually considered by process_run.
        """
        return {}

    def process_run(self, **kwargs):
        """
        Main item of the runner, defines the information that have to be
        post_processed by post_processing.

        Args:
          **kwargs (:py:class:`dict`): keyword arguments as returned from the
            :meth:`pre_processing` method.

        Returns:
          :py:class:`dict`:
               dictionary objects to be passed to post_processing, once the
               dictionary returned by :meth:`pre_processing` has been updated
        """
        return kwargs

    def post_processing(self, **kwargs):
        """
        Post-processing, take the arguments as they are provided by the update
        of :meth:`process_run` and :meth:`pre_processing` methods.

        Returns:
           The final object that each call to the :meth:`run` method is
           supposed to provide.
        """
        return None


class SystemCalculator(Runner):
    """Define a BigDFT calculator.

    Main calculator of BigDFT code. It performs :py:meth:`os.system` calls to
    the main ``bigdft`` executable in the ``$BIGDFT_ROOT`` directory. It is
    designed for two purposes:

      * Run the code in a workstation-based environment, for example within
        notebooks or scripts.
      * Run the code from a python script that is submitted to a batch
        scheduler in a potnentially large-scale supercomputer.

    For triggering the execution, this code gets two variables from the
        environment:

        * The value of ``OMP_NUM_THREADS`` to set the number of
          OMP_NUM_THREADS. If this variable is not present in the environment,
          :class:`SystemCalculator` sets it to the value provided by the
          ``omp`` keyword at initialization.

        * The value of ``BIGDFT_MPIRUN`` to define the MPI execution command.
          If absent, the run is executed simply by ``$BIGDFT_ROOT/bigdft``,
          followed by the command given by post-processing.

    Arguments:
         omp (int): number of OpenMP threads.
            It defaults to the $OMP_NUM_THREADS variable in the environment, if
            present, otherwise it fixes the run to 1 thread.
         mpi_run (str): define the MPI command to be used.
            It defaults to the value $BIGDFT_MPIRUN of the environment, if
            present. When using this calculator into a job submission script,
            the value of $BIGDFT_MPIRUN variable may be set appropriately to
            launch the bigdft executable.
         skip (bool): if ``True``, do not run the calculation if the
           corresponding logfile exists.
         verbose (bool): if ``True`` the class prints out informations about
           the operations that are being performed by the calculator
         dry_run (bool): check the input, estimate the memory but do not
           perform the calculation.
         dry_mpi (int): Number of MPI processes for the estimation of the
           memory when ``dry_run`` is ``True`` (not yet implemented)
         taskgroup_size (int): number of MPI processes of each of the
           taskgroup in the case of a runs_file.

    Warning:
       At the initialization, `SystemCalculator` checks if the environment
         variable $BIGDFT_ROOT is defined.
       This would mean (although not guarantee) that the environment has been
         properly set prior to the evaluation of the python command.
       Also, it checks that the executable file ``bigdft`` might be found in
         the ``$BIGDFT_ROOT/bigdft`` path.

    Example:
        >>> inpdict = { 'dft': { 'ixc': 'LDA' }} #a simple input file
        >>> study = SystemCalculator(omp=1)
        >>> logf = study.run(name="test",input=inpdict)
        Executing command:  $BIGDFT_MPIRUN <path_to_$BIGDFT_ROOT>/bigdft test



    Methods:

         run(name='',run_dir='.',outdir='',run_names='',input=None,posinp='posinp.xyz'):

                Run a calculation building the input file from a dictionary.

                Args:
                   name (str): naming scheme of the run i.e. <name>.yaml is
                        the input file and log-<name>.yaml the output one.
                        Data will then be written in the directory
                        `data-<name>.yaml`, unless the "radical" keyword is
                        specified in the input dictionary.
                   run_dir (str): specify the directory where bigdft will be
                   executed (the input and log file will be created in it)
                             It can be a recursive directory path.
                   outdir (str): specify the output directory for all data
                            coming from bigdft
                   run_names (str): File containing the list of the run ids
                            which have to be launched independently
                            (list in yaml format). This option is not
                            compatible with the ``name`` option.
                   input (:py:class:`dict`): give the input parameters (a
                   dictionary or a list of dictionary). If this parameter is
                           absent it is assumed that an inputfile named
                           `name`.yaml exists in the directory indicated by
                           `run_dir`
                   posinp (file or :py:class:`dict`): indicate the posinp file
                           (atomic position file). It can be either a path or a
                           dictionary in the yaml format.

                Returns:
                     Logfile: Instance of the logfile associated to the run.
                     associated to the run which has been just performed.
                     If the run failed for some reasons, the logfile seem not
                     existing or it cannot be parsed it returns `None`.

                Raises:
                     ValueError: if the logfile does not exists or is not
                     accessible, or if the posinp file does not exists

                Todo:
                     Set the return value of run in the case of a run_file. It
                     should be a list of Logfile classes
    """
    import os
    import shutil

    def __init__(self,
                 omp=os.environ.get('OMP_NUM_THREADS', '1'),
                 mpi_run=os.environ.get('BIGDFT_MPIRUN', ''),
                 dry_run=False, skip=False, verbose=True):
        # Use the initialization from the Runner class (so all options inside
        # __global_options)
        Runner.__init__(self, omp=str(omp), mpi_run=mpi_run,
                        dry_run=dry_run, skip=skip, verbose=verbose)
        # Verify if $BIGDFT_ROOT is in the environment
        assert 'BIGDFT_ROOT' in os.environ
        executable = os.path.join(os.environ['BIGDFT_ROOT'], 'bigdft')
        # the bigdft file should be present in the BIGDFT_ROOT directory
        assert os.path.isfile(executable)
        # Build the command setting the number of omp threads
        self.command = (
            self._global_options['mpi_run'] + ' ' + executable).strip()
        safe_print(
            'Initialize a Calculator with OMP_NUM_THREADS=%s and command %s' %
            (self._global_options['omp'], self.command))

    def pre_processing(self):
        # def run(self, name='', outdir='', run_name='', input={},
        #         posinp=None,**kwargs):
        """
        Process local run dictionary to create the input directory and identify
        the command to be passed

        Returns:
            :py:class:`dict`: dictionary containing the command to be passed to
            :meth:`process_run`

        """
        from futile.Utils import make_dict
        from futile import YamlIO as Y
        self._ensure_run_directory()
        # Create the input file (deepcopy because we modify it)
        inp = self.run_options.get('input', {})
        # from here onwards the local input is a dict and not anymore anothe
        # class
        local_input = make_dict(inp)
        # Add into the dictionary a posinp key
        posinp = self.run_options.get('posinp', None)
        if posinp is not None:
            local_input['posinp'] = self._posinp_dictionary_value(posinp)
        # Creating the yaml input file

        if 'input'in self.run_options:
            input_file = self._get_inputfilename()
        elif 'run_names' not in self.run_options:
            input_file = os.path.join(self.run_dir, 'default.yaml')

        if self.run_options['verbose']:
            safe_print('Creating the yaml input file "%s"' % input_file)
        Y.dump(local_input, filename=input_file)
        return {'command': self._get_command()}

    def process_run(self, command):
        """Finally launch the code.

        Routine associated to the running of the ``bigdft`` executable.

        Arguments:
           command (str): the command as it is set by the ``pre_processing``
             method.

        Returns:
           :py:class:`dict`: The dictionary containing `timedbg` and `logname`
             values to be passed to `post_processing` function
        """
        # check if the debug file will be updated (case of erroneous run)
        timedbg = self._get_debugfile_date()
        verbose = self.run_options['verbose']
        # Set the number of omp threads only if the variable is not present
        # in the environment
        if 'OMP_NUM_THREADS' not in os.environ:
            os.environ['OMP_NUM_THREADS'] = self.run_options['omp']
        if verbose:
            if self.run_dir != '.':
                safe_print('Run directory', self.run_dir)
            safe_print('Executing command: ', command)
        # Run the command
        os.system("cd " + self.run_dir + "; " + command)
        return {'timedbg': timedbg, 'logname': self._get_logname()}

    def post_processing(self, timedbg, logname, command):
        """
        Check the existence and the log file and return an instance logfile.

        Returns:
               (Logfile) Instance of the logfile associated to the run.
               associated to the run which has been just performed.
               If the run failed for some reasons, the logfile seem not
               existing or it cannot be parsed it returns `None`.

        Raises:
               ValueError: if the logfile does not exists or is not accessible

        Todo:
           Set the return value of run in the case of a run_file. It should be
           a list of Logfile classes

        """
        # verify that no debug file has been created
        if self._get_debugfile_date() > timedbg:
            verbose = self.run_options['verbose']
            if verbose:
                safe_print(
                    "ERROR: some problem occured during the execution of the"
                    " command, check the 'debug/' directory and the logfile")
                # the debug file is sane, we may print out the error message
                self._dump_debugfile_info()
            try:
                return Lf.Logfile(logname)
            except IOError:
                return None
        if os.path.exists(logname):
            from futile.Utils import file_time
            inputname = self._get_inputfilename()
            if file_time(logname) < file_time(inputname) and not \
               self.run_options['skip']:
                safe_print("ERROR: The logfile (", logname,
                           ") is older than the inputfile (", inputname, ").")
                return None
            else:
                return Lf.Logfile(logname)
        else:
            raise ValueError("The logfile (", logname, ") does not exist.")

    def _get_command(self):
        name = self.run_options.get('name', '')
        dry_run = self.run_options['dry_run']
        run_names = self.run_options.get('run_names', '')
        outdir = self.run_options.get('outdir', '')
        taskgroup_size = self.run_options.get('taskgroup_size', '')
        # Check if it is a dry run
        if dry_run:
            # Use bigdft-tool (do not use BIGDFT_MPIRUN because it is a python
            # script)
            command = os.path.join(
                os.environ['BIGDFT_ROOT'], 'bigdft-tool') + \
                ' -a memory-estimation -l'
            if name > 0:
                command += ' --name=' + name
        else:
            # Adjust the command line with options
            command = self.command
            if name:
                command += ' -n ' + name
            if run_names:
                command += ' -r ' + run_names
            if outdir:
                command += ' -d ' + outdir
            if taskgroup_size:
                command += ' -t ' + taskgroup_size
            if self.run_options['skip']:
                command += ' -s Yes'
        return command

    def _get_logname(self):
        import os
        outdir = self.run_options.get('outdir', '')
        name = self.run_options.get('name', '')
        logname = 'log-' + name + '.yaml' if name else 'log.yaml'
        if outdir:
            logname = os.path.join(outdir, logname)
        logname = os.path.join(self.run_dir, logname)
        return logname

    def _get_inputfilename(self):
        import os
        name = self.run_options.get('name', '')
        input_file = name + '.yaml' if name else 'input.yaml'
        return os.path.join(self.run_dir, input_file)

    def _ensure_run_directory(self):
        from futile.Utils import ensure_dir
        run_dir = self.run_options.get('run_dir', '.')
        # Restrict run_dir to a sub-directory
        if ("/" in run_dir or run_dir == ".."):
            raise ValueError(
                "run_dir '%s' where bigdft is executed must be a sub-directory"
                % run_dir)
        # Create the run_dir if not exist
        if ensure_dir(run_dir) and self.run_options['verbose']:
            safe_print("Create the sub-directory '%s'" % run_dir)

        self.run_dir = run_dir
        """Run directory.
        str: the directory where the inputfile has been copied to.
              Might be useful to associate to each of the calculation of a
              given run a different directory. Note that this is different than
              setting the ``outdir`` or the ``name`` arguments at it refers
              to the directory of the inputfile.
        Note:
            This is not a global property of the calculator, as the same
            calculator instance might be used for various workflows.
        """

    def _posinp_dictionary_value(self, posinp):
        """
        Create the dictionary value associated to posinp field

        Args:
          posinp (str, dict): path of the posinp file. Might be relative or
          absolute. Copied into `run_dir` if not existing. If it is a
          dictionary, it is a representation of the atomic position.

        Returns:
          str,dict: the value of the key ``posinp`` of the input file, if
          posinp is a string, otherwise the posinp dictionary
        """
        import os
        from futile.Utils import ensure_copy, make_dict
        if isinstance(posinp, dict):
            return make_dict(posinp)
        # Check if the file does exist
        if not os.path.isfile(posinp):
            raise ValueError(
                "posinp: The atomic position file '%s' does not exist"
                % posinp)
        posinpdict = posinp
        posinpfile = os.path.basename(posinp)
        # Copy the posinp if not identical
        cp_posinp = os.path.join(self.run_dir, posinpfile)
        copied = ensure_copy(src=posinp, dest=cp_posinp)
        if copied:
            posinpdict = posinpfile
            if self.run_options['verbose']:
                safe_print("Copy the posinp file '%s' into '%s'" %
                           (posinp, self.run_dir))
        return posinpdict

    def _get_debugfile_date(self):
        """
        Get the information about the debug time of the last file in the
        current directory
        """
        from futile.Utils import file_time
        return file_time(os.path.join(self.run_dir, 'debug',
                                      'bigdft-err-0.yaml'))

    def _dump_debugfile_info(self):
        from futile import YamlIO as Y
        debugfile = os.path.join(self.run_dir, 'debug', 'bigdft-err-0.yaml')
        if os.path.isfile(debugfile):
            debugdict = Y.load(debugfile, doc_lists=False)
            safe_print('The error occured is', self._get_error_key(debugdict))
            safe_print('Additional Info: ', debugdict['Additional Info'])

    def _get_error_key(self, debugdict):
        for key in debugdict:
            if 'Calling sequence' in key:
                continue
            if 'Global dictionary' in key:
                continue
            if 'Additional Info' in key:
                continue
            return key


# Test the calculators
if __name__ == '__main__':

    import yaml
    import matplotlib.pyplot as plt

    basicinput = """
#mode: {method: lj}
logfile: No
dft: { ixc: HF, nspin: 2}
posinp:
   positions:
   - {Be : [0.0, 0.0, 0.0]}#, IGSpin: -1}
   - {Be : [0.0, 0.0, 1.0]}#, IGSpin: 1}
#   properties: {format: yaml}
ig_occupation:
   Atom 1: {2s: {up: 1.0, down: 0.9}, 2p: {up: 0.0, down: 0.2} }
   Atom 2: {2s: {up: 0.9, down: 1.0}, 2p: {up: 0.2, down: 0.0} }

psppar.Be: {Pseudopotential XC: 11}
"""

    # Initialize the calculator
    study = GIBinding()
    inp = yaml.load(basicinput)
    study.set(inp)

    # Perform the first calculation
    out = study.run()
    safe_print('starting energy', out.eKS)
    energy = [out.eKS]
    pos = [1.0]

    # Perform a dissociation curve
    for i in range(10):
        sh = float(i + 1) * 0.02
        inp['posinp']['positions'][-1]['Be'][2] += sh
        study.update(inp)
        out = study.run()
        energy.append(out.eKS)
        pos.append(pos[-1] + sh)
        if study.iproc == 0:
            safe_print('iter', i, 'shift', sh, 'energy', out.eKS)
    out = None
    safe_print('End of the calculations')

    # Plot the dissociation curve
    if study.iproc == 0:
        plt.plot(pos, energy)
        plt.show()
