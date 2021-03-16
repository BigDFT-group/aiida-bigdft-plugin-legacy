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


def _get_target(pv, cutoff, rand):
    from numpy.random import choice
    # Find the least pure fragment to modify
    target = min(pv, key=pv.get)
    target_purity = pv[target]

    if abs(target_purity) >= cutoff and rand:
        target_list = [pid for pid in pv if abs(pv[pid]) > cutoff]
        plist = [abs(x) for x in pv.values() if abs(x) > cutoff]
        normalize = sum(plist)
        plist = [x/normalize for x in plist]
        target = choice(target_list, p=plist)
    return target, target_purity


def _get_distance_dict(tree, frag_lookup, target, target_frag, distcut, k):
    # Build a list of nearby atoms for consideration
    ndist, nearest = tree.query(target_frag, distance_upper_bound=distcut, k=k)
    ndist = [x for sub in ndist for x in sub]
    nearest = [x for sub in nearest for x in sub]

    # Update that into a list of nearby fragments, and keep track of
    # the closest fragment
    distdict = {}
    for idx, dist in zip(nearest, ndist):
        try:
            fragidx = frag_lookup[idx]
        except IndexError:
            # kdtree returns invalid indices if it can't find enough
            # points.
            continue
        if fragidx == target:
            continue
        if fragidx not in distdict:
            distdict[fragidx] = dist
        elif distdict[fragidx] > dist:
            distdict[fragidx] = dist
    return distdict


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
        from os import environ, pathsep
        from os.path import join
        from functools import partial
        from copy import deepcopy

        # Executables
        self.bigdftroot = environ.get('BIGDFT_ROOT')
        # environ['OMP_NUM_THREADS'] = str(omp)
        if not mpi_run:
            mpi_run = environ.get("BIGDFT_MPIRUN")

        # Load the dictionary that defines all the operations.
        # This path should be easier to specify if you use the python
        # egg setup.
        found = False
        for path in environ.get("PYTHONPATH", '').split(pathsep):
            try:
                with open(join(path, "BigDFT", "Database",
                               "postprocess.yaml")) as ifile:
                    db = YamlIO.load(stream=ifile)[0]
                    found = True
            except IOError:
                continue
        if not found and self.bigdftroot is not None:
            raise IOError("Couldn't find postprocess.yaml. " +
                          "Check your python path.")

        # Create the subroutines of the class from the dictionary
        if self.bigdftroot is not None:
            for action, spec in db.items():
                naction = action.replace("_", "-")
                nspec = deepcopy(spec)
                cmd = join(self.bigdftroot, nspec["tool"])
                if nspec["mpi"]:
                    cmd = mpi_run + " " + cmd
                else:
                    nspec["args"]["mpirun"] = {"default": "" + mpi_run + ""}
                    nspec["args"]["action"] = {"default": "" + naction + ""}
                    my_action = get_python_function(
                        partial(self._invoke_command, cmd), action, nspec)
                    setattr(self, action, my_action)

    def _auto_fragment(self, msys, cutoff, verbose, rand, view, criteria):
        from scipy.spatial import KDTree

        # Setup one big kdtree and reverse lookup.
        poslist = []
        posdict = {}
        frag_lookup = []
        # Note the iteration is not over msys because the view may be a subset.
        for fragid in view.purities:
            posdict[fragid] = []
            for at in msys[fragid]:
                pos = at.get_position()
                poslist.append(pos)
                posdict[fragid].append(pos)
                frag_lookup.append(fragid)
        tree = KDTree(poslist)

        fragment_mapping = {frag: [] for frag in view.purities}
        starting_length = len(fragment_mapping)

        # Cutoff for how far to look for a fragment to add.
        distcut = 10.0
        k = 10

        # Zero the internal bond order so we never merge a fragment with itself
        for fragid in view.bond_orders:
            view.bond_orders[fragid][fragid] = 0

        for i in range(0, starting_length - 1):
            if verbose and starting_length > 20:
                step = int(starting_length / 20.0)
                if i % step == 0:
                    print(str(5 * (i / step)) + "% done",
                          len(fragment_mapping))
            target, target_purity = _get_target(view.purities, cutoff, rand)
            if abs(target_purity) < cutoff:
                break

            distdict = _get_distance_dict(tree, frag_lookup, target,
                                          posdict[target], distcut, k)

            # In case nothing was found, we widen the search distance.
            if len(distdict.keys()) == 0:
                distcut *= 2
                k *= 2
                continue

            # Find the neighbor to add in
            if criteria == 'distance':
                addition = min(distdict, key=distdict.get)
            elif criteria == 'bondorder':
                bo_target = {x: view.bond_orders[target][x] for x in distdict}
                addition = max(bo_target, key=bo_target.get)
            else:
                raise ValueError("Criteria must be distance or bondorder")

            # New purity values.
            new_purity = view.bond_orders[target][addition] + \
                view.bond_orders[addition][target]
            new_purity += view.purities[target]*view.charges[target]/2.0
            new_purity += view.purities[addition]*view.charges[addition]/2.0
            new_purity *= 2.0 / (view.charges[target] + view.charges[addition])
            view.purities[target] = new_purity

            posdict[target] = posdict[target] + posdict[addition]

            # New bond order.
            for fragid in view.bond_orders:
                view.bond_orders[target][fragid] += \
                    view.bond_orders[addition][fragid]
                view.bond_orders[fragid][target] += \
                    view.bond_orders[fragid][addition]
            view.bond_orders[target][target] = 0

            # New charge.
            view.charges[target] += view.charges[addition]

            # update lookups
            for fragid in list(view.bond_orders):
                del view.bond_orders[fragid][addition]
            del view.bond_orders[addition]
            del posdict[addition]
            del view.purities[addition]
            del view.charges[addition]
            fragment_mapping[target] += [addition] + fragment_mapping[addition]
            del fragment_mapping[addition]

            # Update the fragment lookup
            frag_lookup = [target if x == addition else x for x in frag_lookup]

        return fragment_mapping

    def auto_fragment(self, system, view, cutoff, verbose=False,
                      rand=False, criteria="bondorder"):
        """
        Refragment a system based on the results of a previous calculation.

        The auto fragment protocol performs a greedy optimization using either
        the bond order or the distance between fragments as a guide. The
        purity values and the fragment bond orders are essential quantities
        for this calculation, so they can be passed if they are already
        cached.

        By using the `rand` keyword, you can trigger a stochastic
        refragmentation process. If this process is repeated many times, you
        may find a result that improves upon the greedy optimization approach.

        Args:
          system (BigDFT.Systems.System): the system to fragment.
          view (BigDFT.Systems.System): a view of the system.
          cutoff (float): the termination criteria. When the worst fragment is
            more pure than this cutoff, the calculation stops.
          verbose (bool): whether to print out as it proceeds.
          rand (bool): whether to activate the stochastic approach.
          criteria (string): either distance or bondorder.

        Returns:
          dict: a mapping from the old system to the new system where each
          fragment fullfills the purity criteria.
        """
        from copy import deepcopy

        # This deep copy is manually overrided for a view so this is not
        # too slow.
        new_view = deepcopy(view)

        # Auto Fragment.
        fragment_mapping = self._auto_fragment(system, cutoff, verbose,
                                               rand, new_view, criteria)

        # Rename
        remapping = {}
        for fragid, add_list in fragment_mapping.items():
            new_name = "+".join([fragid] + add_list)
            remapping[new_name] = [fragid] + add_list

        return remapping

    def fragment_mask_matrix(self, sys, mat, fragments, log):
        """
        Sometimes we don't want to store an entire matrix, just the parts
        related to some fragments of interest. This routine will mask out a
        matrix, keeping entries only related to the list of fragments
        provided.

        Args:
          sys (BigDFT.Systems.System): the system associated with the matrix.
          mat (scipy.sparse.csr_matrix): the matrix to mask.
          fragments (list): a list of fragment ids to keep.
          log (BigDFT.Logfiles.Logfile): the logfile associated with this
            matrix's calculation.
        Returns:
          (scipy.sparse.csr_matrix): the masked matrix.
        """
        from BigDFT.Systems import System
        from copy import deepcopy
        from scipy.sparse import dok_matrix

        subsys = System()
        for fragid in fragments:
            subsys[fragid] = deepcopy(sys[fragid])

        frag_indices = self.get_frag_indices(subsys, log)
        frag_indices = sum(frag_indices.values(), [])

        imat = dok_matrix(mat.shape)
        for i in frag_indices:
            imat[i, i] = 1

        return imat.dot(mat).dot(imat)

    def compute_fragment_dos(self, frag, log, ks_coeff, eigvals,
                             frag_indices=None, smat=None, assume_pure=False,
                             **kwargs):
        """
        Compute the partial density of states projected onto a given
        fragment.

        Args:
            sys (BigDFT.Fragments.Fragment): the fragment to project on to.
            log (BigDFT.Logfiles.Logfile): the log of the calculation.
            ks_coeff (scipy.sparse.csc_matrix): the matrix of eigenvectors.
            eigvals (list): a list of eigenvalues.
            frag_indices (list): list of indices associated with this
              fragment.
            smat (scipy.sparse.csc_matrix): the overlap matrix.
            assume_pure (bool): an optimization can be performed if we assume
              the target is pure.
            kwargs (dict): any extra options to pass to the DoS constructor.
        Returns:
          (BigDFT.DoS.DoS): a density of states object built using the partial
          density of states.
        """
        from BigDFT.Systems import System
        from BigDFT.DoS import DoS
        # Subsystem with just this fragment.
        subsys = System()
        subsys["FRAG:0"] = frag

        # Optional Argments
        if smat is None:
            smat = self.get_matrix_s(log)
        if frag_indices is None:
            frag_indices = self.get_frag_indices(subsys, log)["FRAG:0"]

        subsmat = smat[frag_indices, :]
        subkmat = ks_coeff[frag_indices, :].T

        vals = []
        for n in range(len(eigvals)):
            work = subkmat[n, :].dot(subsmat)
            if assume_pure:
                sumval = work[:, frag_indices].dot(ks_coeff[frag_indices, n])
            else:
                sumval = work[:, :].dot(ks_coeff[:, n])
            vals.append(sumval[0, 0])

        # Construct the dos object
        dos = DoS(energies=eigvals, norm=[vals], units='AU', **kwargs)
        return dos

    def create_layered_qmmm_system(self, system, target, pairwise_bo, cutoffs,
                                   criteria="bondorder", link_atoms=False):
        """
        Creates a multilayered system suitable for QM/MM calculations.
        For each layer, a suitable QM region is built around it.

        Args:
          system (BigDFT.Systems.System): a System class, already broken up
            into fragments.
          pairwise_bo (dict): pairwise bond order values between all fragments
            in the system.
          target (str): the name of the fragment to treat as the target of
            the qm/mm run.
          cutoffs (list): a list of cutoff value for fragment interactions.
            The number of layers is equal to the number of cutoffs.
          criteria (str): how to determine which atoms are included in the
            QM region. Valid choices are "bondorder" and "distance".
          link_atoms (bool): whether to generate link atoms.

        Returns:
          (list): a list of Systems, one for each QM layer.
          (System): the MM region.

        """
        from copy import deepcopy
        from BigDFT.Systems import System

        # First, we will build the layers using only keys.
        qm_list = []
        qm, mm = self.create_qmmm_system(system, target, pairwise_bo[target],
                                         cutoffs[0], criteria,
                                         link_atoms=False)
        qm_list.append(list(qm.keys()))
        mm_keys = list(mm.keys())

        # To avoid duplicates
        used_list = deepcopy(list(qm.keys()))

        for k in range(1, len(cutoffs)):
            qm_list.append([])
            for subtarget in qm_list[k-1]:
                qm, mm = self.create_qmmm_system(system,
                                                 subtarget,
                                                 pairwise_bo[subtarget],
                                                 cutoffs[k], criteria,
                                                 link_atoms=False)
                qm_list[k] += qm.keys()

            # Remove duplicates
            qm_list[k] = [x for x in qm_list[k] if x not in used_list]
            qm_list[k] = list(set(qm_list[k]))
            used_list = list(set(used_list + qm_list[k]))

            mm_keys = [x for x in mm_keys if x not in qm_list[k]]

        # Next, we will convert those keys into systems.
        qmsys = []
        for j in range(0, len(cutoffs)):
            qmsys.append(System())
            for k in qm_list[j]:
                qmsys[j][k] = deepcopy(system[k])
        mmsys = System()
        for k in mm_keys:
            mmsys[k] = deepcopy(system[k])

        # Add the link atoms.
        if link_atoms:
            # One big system to add the link atoms.
            mergesys = System()
            for k in range(0, len(cutoffs)):
                for fragid, frag in qmsys[k].items():
                    mergesys[fragid] = frag
            linksub, remove = self.generate_link_atoms(system, mergesys)

            # Insert the link atoms into the correct layers.
            for fragid, frag in linksub.items():
                for at in frag:
                    if at.is_link:
                        for k in range(0, len(cutoffs)):
                            if fragid in qmsys[k]:
                                qmsys[k][fragid] += [at]
                                break

        return qmsys, mmsys

    def create_qmmm_system(self, system, target, bond_order, cutoff,
                           criteria="bondorder", link_atoms=False):
        """
        Creates a system suitable for QM/MM calculations.

        Args:
          system (System): a System class, already broken up into fragments.
          target (str): the name of the fragment to treat as the target of
            the qm/mm run.
          bond_order (dict): bond order values between the target fragment and
            all other fragments in the system.
          cutoff (float): a cutoff value for fragment interactions.
          criteria (str): how to determine which atoms are included in the
            QM region. Valid choices are "bondorder" and "distance".
          link_atoms (bool): whether to generate link atoms.

        Returns:
          (System): the QM region.
          (System): the MM region.
        """
        from BigDFT.Systems import System
        from BigDFT.Fragments import pairwise_distance
        from copy import deepcopy

        qmsys = System()
        mmsys = deepcopy(system)

        if criteria == "bondorder":
            # Order the fragments by their spillage value
            sort_spill = sorted(bond_order.keys(),
                                key=lambda x: bond_order[x], reverse=True)

            # Iterate until we reach the cutoff
            cumsum = sum(bond_order.values())
            for key in sort_spill:
                qmsys[key] = deepcopy(system[key])
                del mmsys[key]

                # check cutoff
                cumsum -= bond_order[key]
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

        # Double check that the target is in qmsys.
        # (this can happen if the target is not a pure fragment)
        if target not in qmsys:
            qmsys[target] = system[target]

        if link_atoms:
            qmsys, remove = self.generate_link_atoms(system, qmsys)

        return qmsys, mmsys

    def fragment_bond_order(self, sys, fraglist1, fraglist2, log, kxs=None,
                            frag_indices=None):
        """
        Computes "bond order" between two sets of fragments using the method of
        Mayer. For two atomic fragments, this would describe the bond
        multiplicity of the covalent bond.

        Args:
          sys (BigDFT.Systems.System): the system containing the fragments
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
        from scipy.sparse import csr_matrix, csc_matrix
        # Get the metadata
        if kxs is None:
            kxs = self.get_matrix_kxs(log)

        if frag_indices is None:
            frag_indices = self.get_frag_indices(sys, log)

        slist1 = sorted(fraglist1, key=lambda x: min(frag_indices[x]))
        slist2 = sorted(fraglist2, key=lambda x: min(frag_indices[x]))
        bond_orders = units_quadratic_traces(slist1, slist2, frag_indices,
                                             csc_matrix(kxs),
                                             csr_matrix(kxs), 0.5, 0.5)

        return bond_orders

    def fragment_interaction_energy(self, sys, fraglist1, fraglist2, log,
                                    frag_indices=None, sinvh=None, kxs=None):
        """
        Compute the interaction energies between two sets of fragments.

        Args:
          fraglist1 (list): a list of fragments to compute the interaction
            energy of.
          fraglist2 (list): a list of fragments to compute the interaction
            energy between.
          log (BigDFT.Logfiles.Logfile): the log describing a calculation.
          frag_indices (dict): the matrix indices associated with each
            fragment.
          sinvh (scipy.sparse.csc_matrix): the matrix S^{-1}*H, which might be
            already computed to reduce I/O time.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.

        Returns:
          (float): the projected energy.
        """
        from scipy.sparse import csr_matrix, csc_matrix
        # Get the metadata
        if kxs is None:
            kxs = self.get_matrix_kxs(log)
        if sinvh is None:
            sinvh = self.get_matrix_sinvh(log)

        if frag_indices is None:
            frag_indices = self.get_frag_indices(sys, log)

        interaction_energy = units_quadratic_traces(fraglist1, fraglist2,
                                                    frag_indices,
                                                    csc_matrix(sinvh),
                                                    csr_matrix(kxs))
        return interaction_energy

    def fragment_population(self, sys, log, frag_indices=None, kxs=None):
        """
        Performs Mulliken population analysis on a fragment, in case charges
        haven't been computed by doing a multipole analysis.

        Args:
          sys (BigDFT.Systems.System): the system to compute the population of.
          log (BigDFT.Logfiles.Logfile): the log describing a calculation.
          frag_indices (dict): the matrix indices associated with each
            fragment.
          kxs (scipy.sparse.csc_matrix): the matrix K*S, which might be already
            computed to reduce I/O time.

        Return:
          (dict): a mapping from fragment ids to charges.
        """
        from numpy import trace

        if kxs is None:
            kxs = self.get_matrix_kxs(log)
        if frag_indices is None:
            frag_indices = self.get_frag_indices(sys, log)

        charges = {}
        for fragid in sys:
            smat = kxs[:, frag_indices[fragid]]
            smat = smat[frag_indices[fragid], :]
            charges[fragid] = sum(x.nel for x in sys[fragid]) - \
                trace(smat.todense())

        return charges

    def generate_link_atoms(self, fullsys, subsys, distcut=6.0):
        """
        This routine adds link atoms to a subsystem based on the bond
        order of a full system. Link atom positions are automatically adjusted
        based on the length of some standard bonds.

        Args:
          fullsys (BigDFT.Systems.System): the full system that the subsystem
            is embedded into.
          subsys (BigDFT.Systems.System): the embedded system which needs
            link atoms.
          distcut (float): this cutoff is the largest distance value we expect
            allow a bond to be.

        Returns:
          (BigDFT.Systems.System): the subsystem with link atoms added.
          (BigDFT.Systems.System): a system which has the atoms that were
            removed and replaced with link atoms.
        """
        from BigDFT.Systems import System
        from BigDFT.Fragments import Fragment
        from copy import deepcopy
        from numpy.linalg import norm
        from numpy import array
        from warnings import warn

        # Bond lengths for link atoms.
        bond_lengths = {"C": 2.0598, "N": 1.90862, "O": 1.81414, "F": 1.73855,
                        "P": 2.68341, "S": 2.53223, "Cl": 2.39995,
                        "Br": 2.66451, "I": 3.04246}

        if fullsys.conmat is None:
            raise ValueError("Generating link atoms requires connectivity"
                             " information")

        linksys = deepcopy(subsys)
        removesys = System()

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
                        newat.is_link = True
                        # If possible we adjust the position for a reasonable
                        # bond length
                        conat = fullsys[fragid][i]
                        if conat.sym in bond_lengths:
                            pos1 = array(conat.get_position("bohr"))
                            pos2 = array(newat.get_position("bohr"))
                            vec = pos2 - pos1
                            vec *= bond_lengths[conat.sym]/norm(vec)
                            newpos = [x + y for x, y in zip(pos1, vec)]
                            newat.set_position(newpos, units="bohr")
                        else:
                            warn(conat.sym+"bondlength unknown", UserWarning)
                        linklist.append(newat)
                    elif bv > 1:
                        print("Not yet implemented double/triple bonds.", bv)
                        raise NotImplementedError
                    if ft[0] not in removesys:
                        removesys[ft[0]] = Fragment()
                    removesys[ft[0]] += [fullsys[ft[0]][ft[1]]]

            # Add those atoms to the fragment
            for link in linklist:
                linksys[fragid] += [link]

        return linksys, removesys

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
        from scipy.sparse import csc_matrix

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
        purity_array = _compute_purity_values(csc_matrix(kxs), charges, fidx)

        return purity_array

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

    def get_ks_coeff(self, log):
        """
        Retrieve the Kohn-Sham coefficients and the eigenvalues in matrix
        form.

        log (Logfile): instance of a Logfile class

        Returns:
          (scipy.sparse.csc_matrix): the matrix of coefficients
          (list): list of eigenvalues.
        """
        from os.path import join
        from scipy.sparse import coo_matrix, csc_matrix

        # Define the input files.
        data_dir = _get_datadir(log)
        infile = join(data_dir, "KS_coeffs.bin")

        evals = []
        rows = []
        cols = []
        vals = []
        with open(infile) as ifile:
            nspin, ntmb, nfvctr = [int(x) for x in next(ifile).split()[:3]]

            # Values
            for i in range(0, nfvctr):
                evals.append(float(next(ifile).split()[0]))

            # Vectors
            for line in ifile:
                coeff, j, i = line.split()[:3]
                rows.append(int(j)-1)
                cols.append(int(i)-1)
                vals.append(float(coeff))

        mat = coo_matrix((vals, (rows, cols)), shape=(ntmb, nfvctr))

        return csc_matrix(mat), evals

    def get_matrix_kxs(self, log):
        """
        Computes the matrix K*S, the mulliken version of the density matrix,
        and loads it into memory.

        log (Logfile): instance of a Logfile class

        Returns:
          (scipy.sparse.csc_matrix): the matrix K*S
        """
        from BigDFT.Spillage import serial_compute_puritybase
        from scipy.sparse import csc_matrix
        from scipy.io import mmread, mmwrite
        from os.path import join

        def generate_kxs(kxsfile, k_and_s):
            kfile, sfile = k_and_s
            # Then compute with python
            kxs = serial_compute_puritybase(sfile.obj, kfile.obj)
            mmwrite(kxsfile, kxs)

        kxsfile = join(_get_datadir(log), 'kxs_sparse.mtx')

        try:
            kxs = self.derived_quantity(
                      kxsfile, log, ['density_kernel', 'overlap'],
                      generate_kxs)
        except IOError:
            kxs = self.derived_quantity(
                      kxsfile, log,
                      ['density_kernel_sparse', 'overlap_sparse'],
                      generate_kxs)

        # validate the node and read the matrix
        kxs.validate()

        return csc_matrix(mmread(kxsfile))

    def derived_quantity(self, name, log, files, generator):
        """
        Defines a quantity that can be derived from the ccs files
        in the files list. Requires the generator function.
        """
        from futile.Utils import Node, more_recent_than_parent
        from os.path import join, basename

        def more_recent_than_logfile(filename):
            logfile = join(log.srcdir, basename(log.label))
            return more_recent_than_parent(filename, logfile)

        def validate_all(h_and_s):
            ok = True
            for m in h_and_s:
                ok = ok and m.validate()
                if not ok:
                    break
            return ok

        def more_recent_than_parents(spillfile, node):
            ok = True
            for m in node:
                ok = ok and more_recent_than_parent(spillfile, m.obj)
                if not ok:
                    break
            return ok

        newfile = join(_get_datadir(log), name)
        # in case there is no need to rerun the generation
        newnode = Node(newfile, valid_if=more_recent_than_logfile)
        if newnode.valid:
            return newnode

        nodefiles = [self.ccs_node(log, f) for f in files]

        node = Node(nodefiles, valid_if=validate_all)

        newnode = Node(newfile)
        newnode.requires(node, generator=generator)
        newnode.valid_if(more_recent_than_parents)

        try:
            node.validate()
        except AttributeError:
            raise IOError()

        return newnode

    def ccs_node(self, log, matrixname):
        """
        Defines the  protocol to validate the ccs matrix file
        """
        from os.path import join, isfile
        from futile.Utils import Node, non_null_size

        data_dir = _get_datadir(log)
        mtxfile = join(data_dir, matrixname+'.mtx')
        if isfile(mtxfile):  # matrix market output
            return Node(mtxfile, valid_if=non_null_size)
        if log.log["lin_general"]["output_mat"] == 4:
            mpifile = join(data_dir, matrixname+".mpi")
        else:
            mpifile = None
        txtfile = join(data_dir, matrixname+".txt")
        ccsfile = join(data_dir, matrixname+".ccs")
        return file_dependency_queue(self, ccsfile, txtfile, mpifile=mpifile)

    def get_matrix_sinvh(self, log):
        """
        Computes the matrix S^{-1}*H, the mulliken version of the spillage
        matrix, and loads it into memory.

        log (Logfile): instance of a Logfile class

        Returns:
          (scipy.sparse.csc_matrix): the matrix S^{-1}*H
        """
        from BigDFT.Spillage import serial_compute_spillagebase
        from scipy.sparse import csc_matrix
        from scipy.io import mmread, mmwrite
        from os.path import join

        def generate_sinvh(spillfile, h_and_s):
            houtfile, soutfile = h_and_s
            # Then compute with python
            sinvh = serial_compute_spillagebase(soutfile.obj, houtfile.obj)
            mmwrite(spillfile, sinvh)

        spillfile = join(_get_datadir(log), 'sinvh_sparse.mtx')

        try:
            sinvh = self.derived_quantity(spillfile, log,
                                          ['hamiltonian', 'overlap'],
                                          generate_sinvh)
        except IOError:
            sinvh = self.derived_quantity(spillfile, log,
                                          ['hamiltonian_sparse',
                                           'overlap_sparse'],
                                          generate_sinvh)

        # validate the node and read the matrix
        sinvh.validate()

        return csc_matrix(mmread(spillfile))

    def _get_matrix(self, name, log):
        from os.path import join, isfile
        from BigDFT.Spillage import read_ccs
        from scipy.io import mmread
        from scipy.sparse import csc_matrix
        # Define the input files.
        data_dir = _get_datadir(log)

        # Check for mtx file
        mtxfile = join(data_dir, name+".mtx")
        if isfile(mtxfile):
            return csc_matrix(mmread(mtxfile))

        soutfile = join(data_dir, name+".ccs")

        # Convert if necessary
        if not isfile(soutfile):
            sfile = join(data_dir, name+".txt")
            # Convert to text format if necessary
            if log.log["lin_general"]["output_mat"] == 4:
                infile = join(data_dir, sfile.replace(".txt", ".mpi"))
                if not isfile(sfile):
                    self.convert_matrix_format(conversion="binary_to_bigdft",
                                               infile=infile, outfile=sfile)
            # Convert to CCS
            soutfile = join(data_dir, name+".ccs")
            self.convert_matrix_format(conversion="bigdft_to_ccs",
                                       infile=sfile, outfile=soutfile)

        return read_ccs(soutfile)

    def get_matrix_s(self, log):
        """
        Read the overlap matrix into memory.

        log (Logfile): instance of a Logfile class

        Returns:
          (scipy.sparse.csc_matrix): the matrix S
        """
        return self._get_matrix("overlap_sparse", log)

    def get_matrix_h(self, log):
        """
        Read the hamiltonian matrix into memory.

        log (Logfile): instance of a Logfile class

        Returns:
          (scipy.sparse.csc_matrix): the matrix H
        """
        return self._get_matrix("hamiltonian_sparse", log)

    def get_matrix_k(self, log):
        """
        Read the density matrix into memory.

        log (Logfile): instance of a Logfile class

        Returns:
          (scipy.sparse.csc_matrix): the matrix K
        """
        return self._get_matrix("density_kernel_sparse", log)

    def _invoke_command(self, command, **kwargs):
        from futile.Utils import option_line_generator
        _system_command(command, option_line_generator(**kwargs))


def file_dependency_queue(btool, ccsfile, txtfile, mpifile=None):
    from futile.Utils import Node, non_null_size, more_recent_than_parent

    def binary_to_bigdft(node, parent):
        btool.convert_matrix_format(conversion="binary_to_bigdft",
                                    infile=parent, outfile=node)

    def bigdft_to_ccs(node, parent):
        btool.convert_matrix_format(conversion="bigdft_to_ccs",
                                    infile=parent, outfile=node)

    if mpifile is not None:
        mpi = Node(mpifile, valid_if=non_null_size)
        txt = Node(txtfile)
        txt.requires(mpi, generator=binary_to_bigdft)
        txt.valid_if(more_recent_than_parent)
    else:
        txt = Node(txtfile, valid_if=non_null_size)

    ccs = Node(ccsfile)
    ccs.requires(txt, generator=bigdft_to_ccs)
    ccs.valid_if(more_recent_than_parent)

    return ccs


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


def superunits_purities(bo, pv, q_units, superunits, q_superunits):
    """
    Compute the purity values of the superunits described as a unification of
    the units

    Args:
       bo (matrix-like): Fragment Bond Orders of the Units
       pv (array-like): Purity values of the Units
       superunits (dict): lookup dictionary containing the list
         of units per superunit
       q_units (array-like): charges of the units
       q_superunits (array-like): charges of the superunits

    Returns:
       dict: purities of the superunits
    """
    purities = {}
    for fragid, atlist in superunits.items():
        purity = 0

        # Sum the bond orders
        for a1 in atlist:
            for a2 in atlist:
                if a1 != a2:
                    purity += bo[a1][a2]

        # Sum the un-normalized atomic purity values
        for a in atlist:
            charge = q_units[a]  # sum([x.nel for x in atsys["ATOM:"+str(a)]])
            purity += pv[a] * charge / 2.0

        # Normalize
        purity *= 2.0 / q_superunits[fragid]
        # (sum([x.nel for x in resys[fragid]]))
        purities[fragid] = purity
    return purities


def units_quadratic_traces(unitlist1, unitlist2, frag_indices, mat1, mat2,
                           a1=1.0, a2=1.0):
    from numpy import einsum
    # Double loop over the pairs to compute.
    interaction = {}
    for fragid1 in unitlist1:
        interaction[fragid1] = {}
        idx1 = frag_indices[fragid1]
        lmat1 = mat1[:, idx1].todense()  # smat1
        rmat1 = mat2[idx1, :].todense()  # smat2.T

        for fragid2 in unitlist2:
            idx2 = frag_indices[fragid2]

            # Optimization for the case where the indices are consecutive
            # values.
            sidx2 = sorted(idx2)
            mindx = min(sidx2)
            maxdx = max(sidx2)
            if len(sidx2) == maxdx - mindx + 1:
                lmat2 = lmat1[mindx:maxdx+1, :]  # smat12
                rmat2 = rmat1[:, mindx:maxdx+1]  # smat21.T
            else:
                lmat2 = lmat1[idx2, :]  # smat12
                rmat2 = rmat1[:, idx2]  # smat21.T

            interaction[fragid1][fragid2] = a1*a2*einsum("ij,ji->", lmat2,
                                                         rmat2)
    return interaction


def superunits_quadratic_quantities(bo, superunits):
    """
    Compute quantities that transforms like the bond orders of the superunits
    described as a unification of the units

    Args:
       bo (matrix-like): Quantities like Fragment Bond Orders of the Units
       superunits (dict): lookup dictionary containing the list
         of units per superunit

    Returns:
       dict: quantities of the superunits
    """
    bo_new = {}
    for f1, atlist1 in superunits.items():
        for f2, atlist2 in superunits.items():
            order = 0
            for a1 in atlist1:
                for a2 in atlist2:
                    order += bo[a1][a2]
            bo_new.setdefault(f1, {})[f2] = order
    return bo_new


def systems_heatmap(data, restrict_to=None, axs=None, columns=None, **kwargs):
    """
    Create a heatmap for a set of systems.

    Args:
        data (dict): a dictionary mapping system names to another dictionary
          which maps fragment ids to property values.
        restrict_to (list): a list of lists saying what fragments we should
          draw. It is a list of lists instead of just a list because this
          way we can put a separator between values (for example, when making
          a jump in the sequence).
        axs (matplotlib.axis): the axis to draw on.
        columns (list): the order of the columns on the x axis
        kwargs (dict): any extra argument you wish to pass to matplotlib's
          imshow command.

    Returns:
        a reference to the imshow value.
    """
    from numpy import zeros
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Handle optional axis
    if axs is None:
        fig, axs = plt.subplots(1, 1, figsize=(12, 8))

    # Handle the restrict_to option.
    if restrict_to is not None:
        keylist = [x for x in [y for y in restrict_to]]
        keylist = [x for y in restrict_to for x in y]
    else:
        key = list(data.keys())[0]
        keylist = list(data[key])

    cols = columns if columns is not None else data.keys()

    # Create the matrix of values and plot it.
    mat = zeros((len(keylist), len(cols)))
    for i, keyi in enumerate(cols):
        for j, keyj in enumerate(keylist):
            mat[j, i] = data[keyi][keyj]
    im = axs.imshow(mat, **kwargs)

    # Labels
    if restrict_to is not None:
        tickloc = []
        ticklab = []
        i = 0
        for x in restrict_to:
            for y in x:
                tickloc.append(i)
                ticklab.append(y)
                i += 1
            if x == restrict_to[-1]:
                break
            tickloc.append(i-0.5)
            ticklab.append("-------")
        axs.set_yticks(tickloc)
        axs.set_yticklabels(ticklab, family='monospace')
    else:
        axs.set_yticks(range(mat.shape[0]))
        axs.set_yticklabels(keylist, family='monospace')
    axs.set_xticks(range(mat.shape[1]))
    axs.set_xticklabels(list(cols))

    for tick in axs.get_xticklabels():
        tick.set_rotation(75)

    # Draw the colorbar
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="2%", pad=0.15)
    plt.colorbar(im, cax=cax)

    # Draw lines that separate out the boxes.
    for i in range(mat.shape[0]-1):
        axs.axhline(i+0.45, color='k', linestyle='-', linewidth=0.5)
    for i in range(mat.shape[1]-1):
        axs.axvline(i+0.48, color='k', linestyle='-', linewidth=0.5)

    # Not having a fixed aspect is good for larger plots.
    axs.set_aspect("auto")

    return im


def _example():
    """
    Postprocessing Example
    """
    from BigDFT.Systems import System, FragmentView
    from BigDFT.Fragments import Fragment
    from BigDFT.IO import XYZReader
    from BigDFT.Calculators import SystemCalculator
    from BigDFT.Inputfiles import Inputfile
    from scipy.linalg import eigh
    from copy import deepcopy

    # Create a system
    sys = System()
    sys["FRA:0"] = Fragment()
    with XYZReader("CH2") as ifile:
        for at in ifile:
            sys["FRA:0"].append(at)
    sys["FRA:1"] = Fragment()
    with XYZReader("CH3") as ifile:
        for at in ifile:
            sys["FRA:1"].append(at)
    sys["FRA:1"].translate([0, 0, -3])
    sys["FRA:2"] = deepcopy(sys["FRA:0"])
    sys["FRA:2"].translate([0, 0, 3])
    sys["FRA:2"].rotate(y=150, units="degrees")
    sys["FRA:3"] = deepcopy(sys["FRA:0"])
    sys["FRA:3"].translate([3, 0, 1.5])

    print(list(sys))

    # Run a calculation. `write_support_function_matrices` and `linear` are
    # key. You also should be sure to set the atom multipoles.
    inp = Inputfile()
    inp.set_xc("PBE")
    inp.set_hgrid(0.4)
    inp.write_support_function_matrices()
    inp["import"] = "linear"

    code = SystemCalculator()
    code.update_global_options(verbose=False)
    log = code.run(input=inp, posinp=sys.get_posinp(), run_dir="work")
    sys.set_atom_multipoles(log)

    # Create the post processing tool.
    from BigDFT.PostProcessing import BigDFTool
    btool = BigDFTool()

    # Purity
    purity = btool.run_compute_purity(sys, log)
    print(purity)

    # Charges
    charges = {fragid: sum(at.nel for at in frag)
               for fragid, frag in sys.items()}

    # Bond Orders
    bo = btool.fragment_bond_order(sys, sys.keys(), sys.keys(), log)

    # Population values.
    population = btool.fragment_population(sys, log)
    print(population)

    # These three things define a fragment view.
    view = FragmentView(purity, bo, charges)

    # Auto Fragmentation
    mapping = btool.auto_fragment(sys, view, 0.10)
    print(mapping)

    # This defines a new view.
    new_view = view.refragment(mapping)
    print(new_view.purities)

    # Eigenvalues.
    H = btool.get_matrix_h(log)
    S = btool.get_matrix_s(log)
    w = eigh(H.todense(), b=S.todense(), eigvals_only=True)
    print(w)


if __name__ == "__main__":
    _example()
