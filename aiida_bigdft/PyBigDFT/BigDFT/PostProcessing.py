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

    def _auto_fragment_distance(self, system, log, cutoff, verbose, rand,
                                kxs, pv, frag_indices):
        from copy import deepcopy
        from random import choice

        msys = deepcopy(system)

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

        return msys

    def _auto_fragment_bondorder(self, system, log, cutoff, verbose, rand,
                                 kxs, pv, frag_indices):
        from copy import deepcopy
        from random import choice
        from numpy import infty

        msys = deepcopy(system)

        # Cutoff for how far to look for a fragment to add.
        distcut = 10.0

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

            # Find some nearby neighbors to potentially add
            k = 3*len(msys[target])
            addition_list = msys.get_k_nearest_fragments(target, k=k,
                                                         cutoff=distcut)

            # Find the neighbor to add in
            bo = self.fragment_bond_order(msys, [target], addition_list, log,
                                          kxs=kxs, frag_indices=frag_indices)
            bo = bo[target]

            # In case nothing was found, we widen the search distance.
            if len(bo) == 0:
                distcut = infty
                continue

            addition = max(bo, key=bo.get)

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

        return msys

    def auto_fragment(self, system, log, cutoff, verbose=False, rand=False,
                      refine=False, kxs=None, criteria="distance"):
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
          kxs (scipy.sparse.csr_matrix): the density times the overlap,
            if it has been precomputed.
          criteria (string): either distance or bondorder.

        Returns:
          (System): a refragmented system where each fragment has a purity
          value better than the cutoff.
        """
        from copy import deepcopy
        from BigDFT.Fragments import System, Fragment

        msys = deepcopy(system)
        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        frag_indices = self.get_frag_indices(msys, log)
        self.run_compute_purity(msys, log, kxs=kxs, frag_indices=frag_indices)
        pv = {k: v.purity_indicator for k, v in msys.items()}

        if criteria == "distance":
            msys = self._auto_fragment_distance(msys, log, cutoff, verbose,
                                                rand, kxs, pv, frag_indices)
        elif criteria == "bondorder":
            msys = self._auto_fragment_bondorder(msys, log, cutoff, verbose,
                                                 rand, kxs, pv, frag_indices)
        else:
            raise ValueError("Criteria must be distance or bondorder")

        # Give the fragments a name
        msys = msys.rename_fragments()

        # Refinement Loop
        if refine:
            if verbose:
                print("Entering refine loop")
            msys2 = System()
            for count, fragid in enumerate(msys):
                if verbose and len(msys) > 20:
                    step = int(len(msys) / 20.0)
                    if count % step == 0:
                        print(str(5 * (count / step)) + "% done", len(msys))
                # Split up the subsystem
                subsys = System()
                for i, at in enumerate(msys[fragid]):
                    subsys["ATOM:" + str(i)] = Fragment([at])

                # Autofragment
                msubsys = self.auto_fragment(subsys, log, cutoff, rand=True,
                                             kxs=kxs)

                # Add those fragments to the main system
                for sfrag in msubsys.values():
                    sname = "FRAG:" + str(len(msys2) + 1)
                    msys2[sname] = sfrag

            msys = msys2.rename_fragments()

        return msys

    def set_fragment_connectivity(self, frag, log, kxs=None):
        """
        Compute the connectivty of a given fragment using the fragment bond
        order, and sets that information to frag.conmat.

        Args:
          frag (BigDFT.Fragments.Fragment): the fragment to compute the
            connectivity of.
          log (Logfile): the Logfile of a calculation of the full system.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.
        """
        from BigDFT.Fragments import System
        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        # First convert to an atomic system
        atomsys = System()
        for i, at in enumerate(frag):
            atomsys["ATOM:"+str(i)] = [at]

        # Compute the bond order.
        bo = self.fragment_bond_order(atomsys, atomsys.keys(), atomsys.keys(),
                                      log, kxs=kxs)

        # Convert to the matrix format.
        conmat = {}
        for i in range(0, len(frag)):
            at1 = "ATOM:"+str(i)
            for j in range(0, len(frag)):
                at2 = "ATOM:"+str(j)
                if i == j:
                    continue
                v = round(0.5*(bo[at1][at2] + bo[at2][at1]))
                if v > 0:
                    if i not in conmat:
                        conmat[i] = {}
                    conmat[i][j] = v

        frag.conmat = conmat

    def set_system_connectivity(self, sys, log, kxs=None):
        """
        Compute the connectivty of a given fragment using the fragment bond
        order, and sets that information to frag.conmat.

        Args:
          sys (BigDFT.Fragments.System): the system to compute the
            connectivity of.
          log (Logfile): the Logfile of a calculation of the full system.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.
        """
        from scipy.spatial import KDTree
        from BigDFT.Fragments import System, Fragment, GetFragTuple

        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        sys.conmat = {}

        # Create a KD tree of atoms, and a lookup table back to their
        # fragments.
        poslist = []
        frag_lookup = []
        atom_lookup = []
        within_index_lookup = []
        for fragid, frag in sys.items():
            for i, at in enumerate(frag):
                poslist.append(at.get_position())
                frag_lookup.append(fragid)
                atom_lookup.append(at)
                within_index_lookup.append(i)
        tree = KDTree(poslist)

        # Loop over every fragment and atom.
        for fragid, frag in sys.items():
            sys.conmat[fragid] = []
            for at in frag:
                # First, we create a system which is an atom and its
                # nearest neighbors.
                # For these parameters, it's probably fair to say no
                # atom bonds to more than 5 others or to an atom more
                # than 6 bohr away.
                ndist, nearest = tree.query(at.get_position(), k=5+1,
                                            distance_upper_bound=6.0)
                nsys = System()
                nsys["TARGET:0"] = Fragment([at])
                for idx, dist in zip(nearest[1:], ndist[1:]):
                    if dist == float("inf"):
                        break
                    nsys["ATOM:"+str(idx)] = Fragment([atom_lookup[idx]])

                # Compute the bond order which is now an atomic bond order
                # between this atom and its neighbors.
                bo = self.fragment_bond_order(nsys, ["TARGET:0"],
                                              nsys.keys(), log, kxs=kxs)
                bo = bo["TARGET:0"]
                del bo["TARGET:0"]

                toadd = {}
                for k, v in bo.items():
                    rounded = round(v)
                    if rounded > 0:
                        kt = GetFragTuple(k)
                        atom_number = int(kt[1])
                        lookup_key = (frag_lookup[atom_number],
                                      within_index_lookup[atom_number])
                        toadd[lookup_key] = rounded
                sys.conmat[fragid].append(toadd)

    def create_layered_qmmm_sys(self, system, log, target, cutoff, layers,
                                criteria="bondorder", link_atoms=False,
                                kxs=None):
        """
        Creates a multilayered system suitable for QM/MM calculations.
        For each layer, a suitable QM region is built around it.

        Args:
          system (System): a System class, already broken up into fragments.
          log (Logfile): the Logfile of a calculation of the full system.
          target (str): the name of the fragment to treat as the target of
            the qm/mm run.
          cutoff (float): a cutoff value for fragment interactions.
          layers (int): the number of layeres
          criteria (str): how to determine which atoms are included in the
            QM region. Valid choices are "bondorder" and "distance".
          link_atoms (bool): whether to generate link atoms.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.

        Returns:
          (list): a list of Systems, one for each QM layer.
          (System): the MM region.

        """
        from copy import deepcopy
        from BigDFT.Fragments import System

        # First, we will build the layers using only keys.
        qm_list = []
        qm, mm = self.create_qmmm_system(system, log, target, cutoff,
                                         criteria, link_atoms=False, kxs=kxs)
        qm_list.append(list(qm.keys()))
        mm_keys = list(mm.keys())

        # To avoid duplicates
        used_list = deepcopy(list(qm.keys()))

        for l in range(1, layers):
            qm_list.append([])
            for subtarget in qm_list[l-1]:
                qm, mm = self.create_qmmm_system(system, log, subtarget,
                                                 cutoff, criteria,
                                                 link_atoms=False, kxs=kxs)
                qm_list[l] += qm.keys()

            # Remove duplicates
            qm_list[l] = [x for x in qm_list[l] if x not in used_list]
            qm_list[l] = list(set(qm_list[l]))
            used_list = list(set(used_list + qm_list[l]))

            mm_keys = [x for x in mm_keys if x not in qm_list[l]]

        # Next, we will convert those keys into systems.
        qmsys = []
        for l in range(0, layers):
            qmsys.append(System())
            for k in qm_list[l]:
                qmsys[l][k] = deepcopy(system[k])
        mmsys = System()
        for k in mm_keys:
            mmsys[k] = deepcopy(system[k])

        # Add the link atoms.
        if link_atoms:
            # One big system to add the link atoms.
            mergesys = System()
            for l in range(0, layers):
                for fragid, frag in qmsys[l].items():
                    mergesys[fragid] = frag
            qmsys = self.generate_link_atoms(system, mergesys, log, kxs)

            # Break the system up again by layers.
            for l in range(0, layers):
                for fragid in qmsys[l]:
                    qmsys[l][fragid] = mergesys[fragid]

        return qmsys, mmsys

    def create_qmmm_system(self, system, log, target, cutoff,
                           criteria="bondorder", link_atoms=False,
                           kxs=None):
        """
        Creates a system suitable for QM/MM calculations.

        Args:
          system (System): a System class, already broken up into fragments.
          log (Logfile): the Logfile of a calculation of the full system.
          target (str): the name of the fragment to treat as the target of
            the qm/mm run.
          cutoff (float): a cutoff value for fragment interactions.
          criteria (str): how to determine which atoms are included in the
            QM region. Valid choices are "bondorder" and "distance".
          link_atoms (bool): whether to generate link atoms.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.

        Returns:
          (System): the QM region.
          (System): the MM region.
        """
        from BigDFT.Fragments import System, pairwise_distance
        from copy import deepcopy

        qmsys = System()
        mmsys = deepcopy(system)

        if criteria == "bondorder":
            # Compute the spillage
            interaction = self.fragment_bond_order(
                system, [target], system.keys(), log, kxs=kxs)
            interaction = interaction[target]

            # Order the fragments by their spillage value
            sort_spill = sorted(interaction.keys(),
                                key=lambda x: interaction[x], reverse=True)

            # Iterate until we reach the cutoff
            cumsum = sum(interaction.values())
            for key in sort_spill:
                qmsys[key] = deepcopy(system[key])
                del mmsys[key]

                # check cutoff
                cumsum -= interaction[key]
                if cumsum < cutoff:
                    break

        elif criteria == "distance":
            interaction = {x: pairwise_distance(system[target], system[x])
                           for x in system.keys()}
            for fragid, frag in system.items():
                if interaction[fragid] < cutoff:
                    qmsys[fragid] = deepcopy(frag)
                    del mmsys[fragid]

        else:
            raise ValueError("Criteria must be either distance or bondorder")

        if link_atoms:
            if not kxs:
                kxs = self.get_matrix_kxs(log)
            qmsys = self.generate_link_atoms(system, qmsys, log, kxs)

        return qmsys, mmsys

    def fragment_bond_order(self, sys, fraglist1, fraglist2, log, kxs=None,
                            frag_indices=None):
        """
        Computes "bond order" between two sets of fragments using the method of
        Mayer. For two atomic fragments, this would describe the bond
        multiplicity of the covalent bond.

        Args:
          sys (BigDFT.Fragments.System): the system containing the fragments
            of interest
          fraglist1 (list): a list of fragments to compute the bond order
            of.
          fraglist2 (list): a list of fragments to compute the bond order
            between.
          log (BigDFT.Logfiles.Logfile): the log describing a calculation.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.
          frag_indices (dict): the matrix indices associated with each
            fragment.

        Return:
          (dict): a dictionary of dictionaries, mapping the bond order of each
          fragment in fraglist1 to each fragment in fraglist2.
        """
        # Get the metadata
        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        if frag_indices is None:
            frag_indices = self.get_frag_indices(sys, log)

        # Double loop over the pairs to compute.
        bond_orders = {}
        for fragid1 in fraglist1:
            bond_orders[fragid1] = {}
            idx1 = frag_indices[fragid1]
            smat1 = kxs[:, idx1]

            for fragid2 in fraglist2:
                idx2 = frag_indices[fragid2]
                smat2 = smat1[idx2, :]
                bond_orders[fragid1][fragid2] = sum(x**2 for x in smat2.data)

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

    def generate_link_atoms(self, fullsys, subsys, log, kxs, distcut=6.0):
        """
        This routine adds link atoms to a subsystem based on the bond
        order of a full system.

        Args:
          fullsys (BigDFT.Fragments.System): the full system that the subsystem
            is embedded into.
          subsys (BigDFT.Fragments.System): the embedded system which needs
            link atoms.
          log (BigDFT.Logfiles.Logfile): the log describing a calculation.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.
          distcut (float): this cutoff is the largest distance value we expect
            allow a bond to be.

        Returns:
          (BigDFT.Fragments.System): the subsystem with link atoms added.
        """
        from copy import deepcopy

        if fullsys.conmat is None:
            raise ValueError("Generating link atoms requires connectivity"
                             " information")

        linksys = deepcopy(subsys)

        # Loop over atoms and look for bonds running out of the QM system.
        for fragid in subsys:
            linklist = []
            for i in range(0, len(fullsys[fragid])):
                for ft, bv in fullsys.conmat[fragid][i].items():
                    if ft[0] in subsys.keys():
                        continue
                    if bv == 1:
                        newat = deepcopy(fullsys[ft[0]][ft[1]])
                        # Change the symbol to hydrogen
                        newat.sym = "H"
                        linklist.append(newat)
                    elif bv > 1:
                        print("Not yet implemented double/triple bonds.", bv)
                        raise NotImplementedError

            # Add those atoms to the fragment
            for link in linklist:
                linksys[fragid] += [link]

        return linksys

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
        from BigDFT.Spillage import compute_spillage_values
        from BigDFT.Spillage import serial_compute_spillbase
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
        from BigDFT.Spillage import MatrixMetadata
        from os.path import join

        data_dir = _get_datadir(log)
        metadatafile = join(data_dir, "sparsematrix_metadata.dat")
        metadata = MatrixMetadata(metadatafile)
        frag_indices = metadata.get_frag_indices(sys)

        return frag_indices

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
        from BigDFT.Spillage import serial_compute_puritybase
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
        smat = 0.5 * smat[indices, :]
        purity_values[id] = 2.0 * trace(
            (smat.dot(smat) - smat).todense()) / charges[id]

    return purity_values
