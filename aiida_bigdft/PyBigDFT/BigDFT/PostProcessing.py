"""A module for post processing BigDFT calculations.

"""


def _system_command(command, options):
    """
    Run the command as ``os.system (command + options)``

    Args:
       command (str): the actual command to run.
       options (str): the options to pass to the command.
    """
    from subprocess import call

    command_str = command + " " + options
    call(command_str, shell=True)


def _get_datadir(log):
    """
    Given a log file, this returns the path to the data directory.

    Args:
      log (Logfile): logfile from a BigDFT run.
    Returns:
      str: path to the data directory.
    """
    from os.path import join
    return join(log.srcdir, log.data_directory)


class BigDFTool(object):
    """

    This module defines the actions that can be performed using the
    ``bigdft-tool`` python script. Such a script sometimes invokes the
    ``memguess`` executable, which is executed serially. Sometimes the
    ``utilities`` main program is necessary to take care of more time
    consuming operations which require parallelism. We might also use
    the ``obabel`` wrapper provided by BigDFT.

    Args:
      omp (int): number of OpenMP threads.
        It defaults to the $OMP_NUM_THREADS variable in the environment,
        if present, otherwise it fixes the run to 1 thread.
      mpi_run (str): define the MPI command to be used.
        It defaults to the value $BIGDFT_MPIRUN of the environment, if present.
        When using this calculator into a job submission script, the value of
        $BIGDFT_MPIRUN variable may be set appropriately to launch the bigdft
        executable.
    """

    def __init__(self, omp=1, mpi_run=None):
        from futile.YamlArgparse import get_python_function
        from futile import YamlIO
        from os import environ
        from os.path import join
        from functools import partial
        from copy import deepcopy

        # Executables
        self.bigdftroot = environ['BIGDFT_ROOT']
        # environ['OMP_NUM_THREADS'] = str(omp)
        if not mpi_run:
            mpi_run = environ["BIGDFT_MPIRUN"]

        # Load the dictionary that defines all the operations.
        # This path should be easier to specify if you use the python
        # egg setup.
        with open(join(environ['PYTHONPATH'], "BigDFT", "Database",
                       "postprocess.yaml")) as ifile:
            db = YamlIO.load(stream=ifile)[0]

        # Create the subroutines of the class from the dictionary
        for action, spec in db.items():
            naction = action.replace("_", "-")
            nspec = deepcopy(spec)
            cmd = join(self.bigdftroot, nspec["tool"])
            if nspec["mpi"]:
                cmd = mpi_run + " " + cmd
            else:
                nspec["args"]["mpirun"] = {"default": "" + mpi_run + ""}
            nspec["args"]["action"] = {"default": "" + naction + ""}
            my_action = get_python_function(partial(self._invoke_command, cmd),
                                            action, nspec)
            setattr(self, action, my_action)

    def _invoke_command(self, command, **kwargs):
        from futile.Utils import option_line_generator
        _system_command(command, option_line_generator(**kwargs))

    def _run_fragment_multipoles(self, log, system=None, orbitals=None):
        """
        Performs the actual run of the fragment multipoles. This means we
        process the default parameters, override the parameters if they are
        specified, write all the necessary input files, and then run.

        Returns:
          (str): the name of the file containing the multipoles
        """
        from os.path import join, isfile
        from os import remove
        from inspect import getargspec
        from yaml import dump
        from futile.YamlIO import load

        # Convert the arguments of the function to a dictionary
        args, vargs, keywords, default = getargspec(self.fragment_multipoles)
        options = {a: d for a, d in zip(args, default)}

        # Use the logfile to determine the right matrix format
        format = log.log["lin_general"]["output_mat"]
        if format == 4:
            options["matrix_format"] = "parallel_mpi-native"
            for key in options:
                if ".txt" in options[key]:
                    options[key] = options[key].replace(".txt", ".mpi")

        # Replace the default directory with the appropriate one if it is
        # available
        if log.log["radical"]:
            data_dir = _get_datadir(log)
            for a, d in options.items():
                if a == "mpirun" or a == "action" or a == "matrix_format":
                    continue
                options[a] = join(data_dir, d)

        # Create the frag.yaml file from the provided system.
        if system:
            system.write_fragfile(options["fragment_file"], log)

        # Create the orbitals.yaml file from the provided orbitals.
        if orbitals is None:
            with open(options["orbital_file"], "w") as ofile:
                ofile.write(dump([-1]))

        if isfile(options["log_file"]):
            remove(options["log_file"])
        self.fragment_multipoles(**options)

        return load(options["log_file"], doc_lists=False)

    def set_fragment_multipoles(self, system, log):
        """
        Set the fragment multipoles of a system based on a run.

        The multipoles and purity values of the system are
        updated by this call.

        Args:
          system (System): instance of a System class, which defines the
          fragments we will use.
          log (Logfile): instance of a Logfile class
        """
        mp_data = self._run_fragment_multipoles(log, system)

        # Update multipoles and purity values.
        mp_dict = mp_data["Orbital occupation"][0]["Fragment multipoles"]

        for frag, fdata in zip(system.values(), mp_dict):
            frag.purity_indicator = float(fdata["Purity indicator"])
            frag.q0 = [float(x) for x in fdata["q0"]]
            frag.q1 = [float(x) for x in fdata["q1"]]
            frag.q2 = [float(x) for x in fdata["q2"]]

    def get_frag_indices(self, sys, log):
        """
        Compute a lookup table of matrix indices for each fragment in a
          system.

        Args:
          system (System): instance of a System class, which defines the
            fragmentation.
          log (Logfile): logfile from the run computed on this system.

        Returns:
          (dict): a mapping of fragment ids to lists of indices
        """
        from Spillage import MatrixMetadata
        from os.path import join

        data_dir = _get_datadir(log)
        metadatafile = join(data_dir, "sparsematrix_metadata.dat")
        metadata = MatrixMetadata(metadatafile)
        frag_indices = metadata.get_frag_indices(sys)

        return frag_indices

    def run_compute_purity(self, system, log, kxs=None, frag_list=None,
                           frag_indices=None):
        """
        Compute the purity values of the different system fragments.

        Note that this can also be computed using the fragment multipoles,
        but this provides an implementation for when you don't need those
        values.

        Args:
          system (System): instance of a System class, which defines the
            fragmentation.
          log (Logfile): logfile from the run computed on this system.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.
          frag_list (list): we can also only compute the purity values of some
            select fragments.
          frag_indices (list): the indices of the matrix associated with each
            fragment. This can be precomputed and passed.

        Returns:
          (dict): for each fragment id, what is the purity value.
        """
        # Handle the optional parameters
        if frag_list is None:
            working_sys = system
        else:
            working_sys = {k: system[k] for k in frag_list}
        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        # Get the metadata
        if frag_indices is None:
            fidx = self.get_frag_indices(working_sys, log)
        else:
            fidx = {x: frag_indices[x] for x in working_sys}

        # Compute the purity array
        charges = {}
        for fragid, frag in working_sys.items():
            charges[fragid] = sum(x.nel for x in frag)
        purity_array = _compute_purity_values(kxs, charges, fidx)

        # Set the purity values of the system
        for fragid in working_sys:
            system[fragid].purity_indicator = purity_array[fragid]

        return purity_array

    def get_matrix_kxs(self, log):
        """
        Computes the matrix K*S, the mulliken version of the density matrix,
          and loads it into memory.

        log (Logfile): instance of a Logfile class

        Returns:
          (scipy.sparse.csc_matrix): the matrix K*S
        """
        from os.path import join, isfile
        from os import environ
        from Spillage import serial_compute_puritybase
        from scipy.sparse import csc_matrix
        from scipy.io import mmread
        # Define the input files.
        data_dir = _get_datadir(log)
        sfile = join(data_dir, "overlap_sparse.txt")
        dfile = join(data_dir, "density_kernel_sparse.txt")

        # Convert to text format if necessary
        if log.log["lin_general"]["output_mat"] == 4:
            for fname in ["overlap_sparse.txt", "density_kernel_sparse.txt"]:
                infile = join(data_dir, fname.replace(".txt", ".mpi"))
                outfile = join(data_dir, fname)
                if not isfile(outfile):
                    self.convert_matrix_format(conversion="binary_to_bigdft",
                                               infile=infile, outfile=outfile)

        # Check whether to use python or bigpoly version
        # And then perform the computation of the inverse etc
        if isfile(join(environ['BIGDFT_ROOT'], "BigPoly")):
            options = {}
            options["action"] = "multiply_matrices"
            options["infile"] = sfile
            options["infile2"] = dfile
            options["outfile"] = join(data_dir, "kxs.mtx")

            if not isfile(options["outfile"]):
                self.compute_spillage(**options)

            # Read from file
            kxs = csc_matrix(mmread(options["outfile"]))
        else:  # Compute serially with python
            # First convert to ccs matrix format
            soutfile = join(data_dir, "overlap_sparse.ccs")
            if not isfile(soutfile):
                self.convert_matrix_format(conversion="bigdft_to_ccs",
                                           infile=sfile, outfile=soutfile)
            doutfile = join(data_dir, "density_kernel_sparse.ccs")
            if not isfile(doutfile):
                self.convert_matrix_format(conversion="bigdft_to_ccs",
                                           infile=dfile, outfile=doutfile)
            # Then compute with python
            kxs = serial_compute_puritybase(soutfile, doutfile)

        return kxs

    def run_compute_spillage(self, system, log, targetid):
        """
        Compute a measure of the spillage interaction between fragments.

        Args:
          system (System): instance of a System class, which defines the
            fragments we will use.
          log (Logfile): instance of a Logfile class
          targetid (str): which fragment to compute the spillage of all other
            fragments with. You can either set this or target.
        Returns:
          (tuple): target spillage value and for each fragment id, what the
            spillage value is.
        """
        from os.path import join, isfile
        from os import environ
        from Spillage import compute_spillage_values
        from Spillage import serial_compute_spillbase
        from scipy.sparse import csc_matrix
        from scipy.io import mmread

        # Define the input files.
        data_dir = _get_datadir(log)
        sfile = join(data_dir, "overlap_sparse.txt")
        hfile = join(data_dir, "hamiltonian_sparse.txt")

        # Convert to text format if necessary
        if log.log["lin_general"]["output_mat"] == 4:
            for fname in ["overlap_sparse.txt", "hamiltonian_sparse.txt"]:
                infile = join(data_dir, fname.replace(".txt", ".mpi"))
                outfile = join(data_dir, fname)
                if not isfile(outfile):
                    self.convert_matrix_format(conversion="binary_to_bigdft",
                                               infile=infile, outfile=outfile)
        # Get the metadata
        frag_indices = self.get_frag_indices(system, log)

        # Check whether to use python or bigpoly version
        # And then perform the computation of the inverse etc
        if isfile(join(environ['BIGDFT_ROOT'], "BigPoly")):
            options = {}
            options["action"] = "compute_spillage"
            options["infile"] = hfile
            options["infile2"] = sfile
            options["outfile"] = join(data_dir, "sinvxhfile.mtx")
            options["outfile2"] = join(data_dir, "sinvxh2file.mtx")

            self.compute_spillage(**options)

            # Read from file
            sinvxh = csc_matrix(mmread(options["outfile"]))
            sinvxh2 = csc_matrix(mmread(options["outfile2"]))
        else:
            # First convert to ccs matrix format
            soutfile = join(data_dir, "overlap_sparse.ccs")
            if not isfile(soutfile):
                self.convert_matrix_format(conversion="bigdft_to_ccs",
                                           infile=sfile, outfile=soutfile)
            houtfile = join(data_dir, "hamiltonian_sparse.ccs")
            if not isfile(houtfile):
                self.convert_matrix_format(conversion="bigdft_to_ccs",
                                           infile=hfile, outfile=houtfile)
            # Then compute with python
            sinvxh, sinvxh2 = serial_compute_spillbase(soutfile, houtfile)

        # Compute the spillage array
        spillage_tuple = compute_spillage_values(
            sinvxh, sinvxh2, frag_indices, targetid)

        return spillage_tuple

    def plot_spillage(self, axs, spillvals, colors=None, minval=0.0):
        """
        Plot the spillage values.

        Args:
          axs: the axs we we should plot on.
          spillvals (dict): a dictionary mapping fragments to spillage values.
          colors (dict): you can optionally pass a dictionary which sets
            the colors of the plot. The keys are floating point numbers and
            the values are colors, where all spillage values above that key
            are colored with the value.
        """
        from Fragments import GetFragTuple
        from numpy.ma import masked_less

        axs.set_xlabel("Fragment", fontsize=12)
        axs.set_ylabel("Spillage Values", fontsize=12)
        axs.set_yscale("log")

        svals = []
        labels = []
        for id, val in spillvals.items():
            if abs(val) > minval:
                svals.append(abs(val))
                labels.append(id)

        axs.set_xticks(range(len(labels)))
        sorted_labels = sorted(labels, key=lambda x: int(GetFragTuple(x)[1]))
        axs.set_xticklabels(sorted_labels, rotation=90)
        sorted_values = [v for _, v in sorted(zip(labels, svals),
                                              key=lambda x:
                                              int(GetFragTuple(x[0])[1]))]
        axs.plot(sorted_values, marker='o', color='black', linestyle='--')

        if colors:
            for thresh, col in sorted(colors.items()):
                axs.plot(masked_less(sorted_values, thresh),
                         color=col, marker='o', markersize=12, linestyle='--',
                         markeredgewidth='1.5', markeredgecolor='black')

    def auto_fragment(self, system, log, cutoff, verbose=False, rand=False,
                      refine=False):
        """
        Refragments a system based on the results of a previous calculation.

        The auto fragment protocol searches for the least pure fragment, and
        then merges it with its nearest neighbor. This is repeated until the
        purity values of the entire system are less than the provided cutoff.
        You can also set the rand parameter to have the protocol merge random
        non pure fragments instead.

        The refine option performs a refinement step after the fragmentation
        is complete. For each pure fragment, we break it up again, and then
        use the random merge approach within that fragment.

        Args:
          system (System): instance of a System class, which defines the
            fragments we will use.
          log (Logfile): instance of a Logfile class
          cutoff (float): the termination criteria. When the worst fragment is
            more pure than this cutoff, the calculation stops.
          verbose (bool): whether to print out as it proceeds.
          rand (bool): whether to randomly select fragments to merge or not.
          refine (bool): whether to do the additional refinement step.

        Returns:
          (System): a refragmented system where each fragment has a purity
            value better than the cutoff.
        """
        from copy import deepcopy
        from BigDFT.Fragments import System, Fragment
        from random import choice

        msys = deepcopy(system)
        kxs = self.get_matrix_kxs(log)

        frag_indices = self.get_frag_indices(msys, log)
        self.run_compute_purity(msys, log, kxs=kxs, frag_indices=frag_indices)
        pv = {k: v.purity_indicator for k, v in msys.items()}

        # Initial Fragmentation
        if verbose:
            print("Entering main loop")
        for i in range(0, len(system) - 1):
            if verbose and len(system) > 20:
                step = int(len(system) / 20.0)
                if i % step == 0:
                    print(str(5 * (i / step)) + "% done", len(msys))

            # Find the least pure fragment to modify
            target = min(pv, key=pv.get)
            target_purity = pv[target]
            if abs(target_purity) < cutoff:
                break

            if rand:
                target_list = [pid for pid in pv if abs(pv[pid]) > cutoff]
                target = choice(target_list)

            # Find the neighbor to add in
            addition = msys.get_nearest_fragment(target)

            # Combine
            msys[target] = msys[target] + deepcopy(msys[addition])
            frag_indices[target] = frag_indices[target] + \
                deepcopy(frag_indices[addition])
            del msys[addition]
            del pv[addition]
            del frag_indices[addition]

            self.run_compute_purity(msys, log, kxs=kxs, frag_list=[
                                    target], frag_indices=frag_indices)
            pv[target] = msys[target].purity_indicator

        # Give the fragments a name
        msys = msys.rename_fragments()

        # Refinement Loop
        if refine:
            msys2 = System()
            for fragid in msys:
                # Split up the subsystem
                subsys = System()
                for i, at in enumerate(msys[fragid]):
                    subsys["ATOM:"+str(i)] = Fragment([at])

                # Autofragment
                msubsys = self.auto_fragment(subsys, log, cutoff, rand=True)

                # Add those fragments to the main system
                for sfrag in msubsys.values():
                    sname = "FRAG:"+str(len(msys2)+1)
                    msys2[sname] = sfrag

            msys = msys2.rename_fragments()

        return msys

    def fragment_bond_order(self, sys, fragid1, fraglist, log, kxs=None):
        """
        Computes "bond order" between two fragments using the method of Mayer.
        For two atomic fragments, this would describe the bond multiplicity of
        the covalent bond.

        Args:
          sys (BigDFT.Fragments.System): the system containing the fragments
            of interest
          fragid1 (str): the first fragment.
          fraglist (str): a list of fragments to compute the bond order
            between.
          log (BigDFT.Logfiles.Logfile): the log describing a calculation.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.

        Return:
          (dict): a measure of the bond multiplicity for each fragment.
        """
        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        # Get the metadata
        frag_indices = self.get_frag_indices(sys, log)
        idx1 = frag_indices[fragid1]

        bond_orders = {}
        smat1 = kxs[:, idx1]
        for fragid2 in fraglist:
            idx2 = frag_indices[fragid2]
            smat2 = smat1[idx2, :]
            bond_orders[fragid2] = sum(x**2 for x in smat2.data)

        return bond_orders

    def fragment_population(self, sys, fragid, log, kxs=None):
        """
        Performs Mulliken population analysis on a fragment, in case charges
        haven't been computed by doing a multipole analysis.

        Args:
          sys (BigDFT.Fragments.System): the system containing the fragment
            of interest
          frag (str): the fragment to compute the charge of.
          log (BigDFT.Logfiles.Logfile): the log describing a calculation.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.

        Return:
          (float): the amount of charge on a fragment.
        """
        from numpy import trace

        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        # Get the metadata
        frag_indices = self.get_frag_indices(sys, log)
        indices = frag_indices[fragid]

        smat = kxs[:, indices]
        smat = smat[indices, :]

        return sum(x.nel for x in sys[fragid]) - trace(smat.todense())

    def compute_fragment_smiles(self, frag):
        """
        Computes the SMILES representation of a given fragment.

        This uses openbabel, so you need to either have written it on your
        path or to have built BigDFT with the openbabel option.

        Args:
          frag (BigDFT.Fragments.Fragment): the fragment to compute the
            representation of.

        Return:
          (str): the smiles representation of this molecule.
        """
        from BigDFT.XYZ import write_pdb
        from tempfile import NamedTemporaryFile
        from subprocess import check_output
        from os.path import join

        with NamedTemporaryFile() as ftemp:
            fname = ftemp.name

        write_pdb(fname, frag)
        return check_output([join(self.bigdftroot, "obabel"),
                             "-ipdb", fname, "-osmi"]).split()[0]

    def compare_smiles(self, smiles1, smiles2):
        """
        Computes a similarity score between two smiles representations.

        This uses openbabel, so you need to either have written it on your
        path or to have built BigDFT with the openbabel option.

        Args:
          smiles1 (str): the first smiles string.
          smiles2 (str): the second smiles string.

        Return:
          (float): the similarity score.
        """
        from subprocess import check_output
        from os.path import join

        with open("temp1.smi", "w") as ftemp:
            ftemp.write(smiles1)
            fname1 = ftemp.name
        with open("temp2.smi", "w") as ftemp:
            ftemp.write(smiles2)
            fname2 = ftemp.name

        output = check_output([join(self.bigdftroot, "obabel"),
                               fname1, fname2, "-ofpt"])
        output = [x for x in output.split("\n") if "Tanimoto" in x][0]
        return float(output.split()[-1])


def _compute_purity_values(kxs, charges, frag_indices):
    """
    Computes the actual purity values.

    Args:
      kxs (scipy.sparse.csc): K*S
      charges (dict): the charge to divide by for each framgent.
      frag_indices (dict): indices associated with each fragment.

    Returns:
      (dict): purity of each fragment.
    """
    from numpy import trace
    purity_values = {}

    for id in frag_indices:
        indices = frag_indices[id]
        smat = kxs[:, indices]
        smat = 0.5*smat[indices, :]
        purity_values[id] = 2.0 * trace(
            (smat.dot(smat) - smat).todense()) / charges[id]

    return purity_values
