"""A module to define typical operations that can be done on biological systems

"""
from BigDFT.Systems import System


def load(archive, serialization_version=None, options={}):
    """
    Create a Biosystem from a serialization archive.
    It extracts the archive files in a temporary directory
    and loads the related information. Also it deletes the directory
    once the load is finished.

    Args:
       archive(path): the path of the serialization archive
       serialization_version (str): version of the archive.
           Latest if not specified.

    Returns:
       BioSystem: Instance of the class, loaded from the archive
    """
    from futile.Utils import unpack_tarball
    from shutil import rmtree
    from os.path import join  # ,basename
    # import tempfile
    # import tarfile
    # import yaml
    import json
    tmpdir, files = unpack_tarball(archive, tmpdir_prefix='tmp_BioQM_archive_')
    # # creates the directory
    # tmpdir = tempfile.mkdtemp(prefix='tmp_BioQM_archive_' +
    #                           basename(archive) + '_')
    # # extract the archive
    # arch = tarfile.open(archive)
    # arch.extractall(path=tmpdir)
    # arch.close()
    # load the data
    sys = BioSystem(join(tmpdir, 'posinp.pdb'), **options)
    sys.set_qm_run(join(tmpdir, 'log.yaml'))
    for version, attrs in BioSystemSerialization.cached_attributes.items():
        if serialization_version is not None:
            desired = serialization_version
        else:
            desired = BioSystemSerialization.version
        if not _version_is_compatible(desired, version):
            continue
        for attr in attrs:
            f = open(join(tmpdir, attr+'.json'), 'r')
            value = json.load(f)
            setattr(sys, '_'+attr, value)
            f.close()
    # remove the directory
    rmtree(tmpdir)

    return sys


def _version_is_compatible(desired_version, present_version):
    desired = list(map(int, desired_version.split('.')))
    obtained = list(map(int, present_version.split('.')))
    not_ok = obtained[0] != desired[0]
    if not not_ok:
        not_ok = any([ll < m for ll, m in zip(desired[1:], obtained[1:])])
    return not not_ok


def serialize_biosystem(archive, posinp, logfile, version=None):
    """
    Serialize the biosystem associated to the provided files

    Args:
       archive (path): name of the archive to be employed for dumping
       posinp (path): pdbfile of the structure
       logfile (path): logfile of the structure
       version(str): version of the serialization
    """
    system = BioSystem(posinp)
    system.set_qm_run(logfile)
    if version is not None and _version_is_compatible('2.0', version):
        system.refragment(0.05)  # to create the atomic purities
    # system.refragment(0.025)
    serialization = BioSystemSerialization(system, version=version)
    serialization.dump(archive)


class BioSystemSerialization():
    """
    Serialization of the BioSystem class that would prevent the need to
    copy the entire matrices of the QM calculations

    Args:
        biosys (BioSystem): a instance of the Biosystem class
        version (str): the version of the desired serialization.
           If absent, the default is considered.
    """

    version = '1.1'  #: version of the serialization
    cached_attributes = {'1.0': {'purities': 'purities',
                                 'fragindices': 'frag_indices',
                                 'pairwise_BO': 'bond_orders',
                                 'cached_refragment': '_cached_refragment'},
                         '1.1': {'interactions': 'interactions'},
                         '2.0': {'atomic_purities': '_atomic_purities',
                                 'atomic_BO': 'atomic_bond_orders',
                                 'atomic_interactions': 'atomic_interactions',
                                 'cached_refragment': '_cached_refragment'}}

    def __init__(self, biosys, version=None):
        self.files = {'posinp.pdb': biosys.structure_file_path,
                      'log.yaml': biosys.logfile_path}

        self.true_version = version if version is not None else self.version

        for version, attributes in self.cached_attributes.items():
            if not _version_is_compatible(self.true_version, version):
                continue
            for attr, trueattr in attributes.items():
                setattr(self, attr, getattr(biosys, trueattr))

    def dump(self, archive):
        """
        Create an archive with the entire set of information of the
        Serialization. Such a tarfile should be such that the same
        analysis of the BioSystem is possible
        """
        from futile.Utils import create_tarball, serialize_objects
        # import json
        from BigDFT import Systems

        def if_system(sys):
            return sys.get_posinp()

        # class NumpyEncoder(json.JSONEncoder):
        #     """ Special json encoder for numpy types and System"""
        #     def default(self, obj):
        #         import numpy as np
        #         import Systems
        #         if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
        #                             np.int16, np.int32, np.int64, np.uint8,
        #                             np.uint16, np.uint32, np.uint64)):
        #             return int(obj)
        #         elif isinstance(obj, (np.float_, np.float16, np.float32,
        #                               np.float64)):
        #             return float(obj)
        #         elif isinstance(obj, (np.ndarray,)):
        #             return obj.tolist()
        #         elif isinstance(obj, (set,)):
        #             return list(obj)
        #         elif isinstance(obj, (Systems.System,)):
        #             return obj.get_posinp()
        #         return json.JSONEncoder.default(self, obj)

        dictionaries = {}
        for version, attrs in self.cached_attributes.items():
            if not _version_is_compatible(self.true_version, version):
                continue
            dictionaries.update({att+'.json': getattr(self, att)
                                 for att in attrs})

        system_encoding = {'cls': Systems.System, 'func': if_system}

        objects = serialize_objects(
                    dictionaries, extra_encoder_functions=[system_encoding])

        create_tarball(archive, self.files, objects)


class Sequence():
    """
    Initialize and define a sequence, may be amminoacidic or nitrogenous.

    Args:
       seq(Bio.Seq.Sequence) : the starting sequence
    """

    def __init__(self, seq):
        self.sequence = seq

    def __add__(self, other):
        return Sequence(self.sequence + other.sequence)

    def __radd__(self, other):
        return self if other == 0 else self.__add__(other)

    def display(self, colors=None, labels=None, axs=None, boxcount=None,
                remove_letters=[]):
        """
        Display the sequence in a easy-to-read form

        Args:
            boxcount(int): the number of boxes per row
            remove_letters (list): list of letter positions that will be
                removed
        """
        def formatfun(x):
            return ("{: ^"+str(llen)+"}").format(x)
        # Handle optional axis
        if axs is None:
            import matplotlib.pyplot as plt
            fig, axs = plt.subplots()
        axs.axis('off')
        # bbox = axs.get_window_extent()

        # Handle the sizes based on the presence of labels
        if labels:
            yoffset = 24
            boxcount = 30 if boxcount is None else boxcount
            llen = max([len(x) for x in labels if x is not None])
            xoffset = 3*llen + 5
        else:
            xoffset = 7
            yoffset = 14
            boxcount = 40 if boxcount is None else boxcount

        # Handle the matplotlib options
        matargs = {"size": 18, "va": "center", "ha": "center",
                   "multialignment": "left"}
        matargs["fontdict"] = {"family": "monospace", "weight": "normal"}

        xcoord = 0
        ycoord = 100
        removed = False
        for i, letter in enumerate(str(self.sequence)):
            if i in remove_letters:
                if not removed:
                    matargs['bbox'] = {"linewidth": 0, 'fill': False}
                    axs.text(0.01 * xcoord, 0.01 * ycoord, "-",
                             **matargs)
                    # Update coordinates
                    xcoord += xoffset
                    removed = True
                continue
            removed = False
            label = None if labels is None or labels[i] is None else labels[i]
            fc = "none" if colors is None or colors[i] is None else colors[i]

            if isinstance(fc, tuple):
                lc = fc[1]
                fc = fc[0]
            else:
                lc = fc

            if label:
                matargs["bbox"] = {"fc": lc, "linewidth": 0}
                axs.text(0.01 * xcoord, 0.01 * ycoord + 0.05, formatfun(""),
                         **matargs)
                matargs["bbox"] = {"fc": fc, "linewidth": 0}
                axs.text(0.01 * xcoord, 0.01 * ycoord - 0.05, formatfun(""),
                         **matargs)

                matargs["bbox"] = {"fill": False, "linewidth": 1,
                                   "edgecolor": "k"}
                axs.text(0.01 * xcoord, 0.01 * ycoord, formatfun(label) +
                         "\n" + formatfun(letter), **matargs)
            else:
                matargs["bbox"] = {"fc": fc, "linewidth": 1, "edgecolor": "k"}
                axs.text(0.01 * xcoord, 0.01 * ycoord, letter, **matargs)

            # Update coordinates
            xcoord += xoffset
            if xcoord > boxcount * xoffset:
                xcoord = 0
                ycoord -= yoffset
        return axs


def order_chains(fragments):
    from Bio.PDB.Polypeptide import three_to_one as letter
    chains = {}
    for f in fragments:
        ch, res, num = construct_frag_tuple(f)
        resnum = ':'.join([res, str(num)])
        try:
            lett = letter(res)
        except KeyError:
            lett = 'X'
        chains.setdefault(ch, {})[resnum] = lett
    return chains


def residue_to_chains(chdict, residues):
    chlett = {}
    for ires, res in enumerate(residues):
        ch, resname, num = construct_res_tuple(res)
        resnum = ':'.join([resname, str(num)])
        lett = chdict[ch].pop(resnum)
        if lett != 'X':
            chlett.setdefault(ch, {})[num] = (lett, ires)
    assert all(len(ch) == 0 for ch in chdict.values())
    return chlett


def chains_to_sequences(chlett):
    from Bio.Seq import Seq
    from Bio.Alphabet import ProteinAlphabet
    seqs = []
    chains_to_residues = []
    for ch in sorted(chlett.keys()):
        chain = chlett[ch]
        strseq = ''.join([chain[i+1][0] for i in range(len(chain))])
        seqs.append(Sequence(Seq(strseq, alphabet=ProteinAlphabet)))
        chains_to_residues.append([chain[i+1][1] for i in range(len(chain))])
    return seqs, chains_to_residues


def pdb_sequences(filename):
    """
    Sequences of the chains of the proteins, in the order of the residue_list

    Args:
       filename (str): protein file, in pdb format

    Returns:
       list: strings of the protein sequence in the FASTA sequence order
    """
    import Bio.SeqIO as S
    return [Sequence(rec.seq) for rec in S.parse(filename, format='pdb-atom')]


def read_structure(filename, name='protein'):
    """
    Extract a biological structure from a PDB file. Use the package Bio.PDB

    Args:
       filename (str): path to the file of the structure
       name (str): the name of the structure

    Returns:
       (Bio.PDB.Structure) A structure instance
    """
    import Bio.PDB as pdb
    parser = pdb.PDBParser()
    return parser.get_structure(name, filename)


def residue_list(structure, system, strict=True):
    """
    Identify the residues of the sequence in the fragment of the system

    Args:
       structure (structure): A structure read from the PDB parser of Bio.PDB
       system (BigDFT.Systems.System): a system coming from BigDFT package
       strict (bool): define if the identification has to be performend via
          centroids control or just by the residue name

    Returns:
       tuple: tuple of list the fragments of the `system` in the order of the
           residues of the structure. Some of the fragments may be missing as
           they may have been merged of not associated to elements of the
           sequence. The second list is the list of the ordered residues.
    """
    from BigDFT.Atoms import AU_to_A
    import numpy as np
    reslist = []
    centroids_dict = {frag: system[frag].centroid * AU_to_A for frag in system}
    original_residues = []
    for residue in structure.get_residues():
        restuple = construct_res_tuple(residue)
        centroid = sum([a.coord for a in residue.get_atoms()]) / len(residue)
        found = False
        for frag in centroids_dict:
            fragtuple = construct_frag_tuple(frag)
            if not strict:
                found = restuple == fragtuple
            else:
                found = np.linalg.norm(centroids_dict[frag] - centroid) < 1.e-4
            if found:
                reslist.append(frag)
                centroids_dict.pop(frag)
                break
        if not found:
            reslist.append(None)
        original_residues.append(residue)
    return reslist, original_residues


def sequences_to_residues(reslist, chains):
    """
    Identify for each element of the chains which are the residues that may be
    associated to them.

    Args:
       reslist(list): the list of all the residues of the system
       chains (list): a list containing the sequences of the system

    Returns:
       list: a list of length of the chains indicating the lookup array for the
           sequence
    """
    lookup = [[None for k in c.sequence] for c in chains]
    chain_letters = []
    for ires, res in enumerate(reslist):
        chain_id = res.full_id[2]
        if chain_id not in chain_letters:
            chain_letters.append(chain_id)
        ch = chain_letters.index(chain_id)
        pos = res.full_id[3][1] - 1
        if ch < len(lookup) and ch >= 0 and pos < len(lookup[ch]) and pos >= 0:
            lookup[ch][pos] = ires
    return lookup


class Graph():
    """
    A Graph representation of a System or a subportion of it.

    Identify the connectivity of the system given a threshold, possibly
    restricted to a subset of fragments
    Args:
       cutoff (float) : the threshold of the cumulative fragment bond order
        that is applied to define the connectivity
       restrict_to (list): list of the fragment id to which calculate
        the connectivity
    """

    def __init__(self, fragkeys, bond_orders, frag_labels, cutoff=0.01,
                 restrict_to=None):
        import numpy as np
        self.connectivity_matrix = graph_bond(
            fragkeys, cutoff, bond_orders, restrict_to=restrict_to)

        self.nw, self.target_nodes, self.target_edges = get_BO_network(
                                                    self.connectivity_matrix.T)

        self.labels = {i: l for i, l in enumerate(frag_labels)
                       if i in self.target_nodes}

        if restrict_to is not None:
            self.restrict_nodes = [fragkeys.index(f) for f in restrict_to]
        else:
            self.restrict_nodes = None

        frags = fragkeys

        self.edge_weights = []
        for i, j in self.nw.edges:  # self.target_edges:
            wgt = bond_orders[frags[i]][frags[j]] if i != j else -1
            self.edge_weights.append(wgt)

        self.edge_weights = np.array(self.edge_weights)
        wh = np.where(self.edge_weights == -1)
        mx = np.max(self.edge_weights)
        self.edge_weights[wh] = mx

    def display(self, colors=None, edge_colorcode='Reds', axs=None,
                edge_labels=None, **kwargs):
        """
        Draw the graph in the Kamada Kawai representation.
        Employ colors if present.

        Args:
           colors (list): the list of the colors that have to be employed for
              all the graph nodes. If the target nodes are present only the
              selected nodes are employed
           edge_colors (list): the list of the colors that have to be employed
              for each of the edges of the node
           axs (matplotlib.axis): the axis to draw on (optional)
           edge_labels (dict): a dictionary from tuples of node labels to
             a given label.
           kwargs (dict): any other argument you wish to pass to networkx.

        """
        from networkx import kamada_kawai_layout, draw_networkx
        from networkx import draw_networkx_edge_labels, draw_networkx_nodes

        # Handle optional axis
        if axs is None:
            import matplotlib.pyplot as plt
            fig, axs = plt.subplots(1, 1, figsize=(18, 12))

        # Setup colors
        node_color = 'r' if colors is None else [
            colors[i] for i in self.target_nodes]

        # Handle kwargs
        args = {"labels": self.labels, "node_size": 800,
                "font_size": 10, "font_family": 'serif', "linewidths": 4,
                "font_weight": 'bold', "alpha": 0.95,
                "ax": axs, "edge_color": 'k'}
        args.update(kwargs)

        layout = kamada_kawai_layout(self.nw)

        if self.restrict_nodes is not None:
            node_color = 'r' if colors is None else [
                colors[i] for i in self.restrict_nodes]
            draw_networkx_nodes(self.nw, layout, node_shape='s',
                                node_color=node_color,
                                nodelist=self.restrict_nodes, **args)

        node_color = 'r' if colors is None else [
            colors[i] for i in self.target_nodes]
        draw_networkx(self.nw, layout, node_color=node_color,
                      nodelist=self.target_nodes, **args)

        if edge_labels is not None:
            draw_networkx_edge_labels(self.nw, layout,
                                      nodelist=self.target_nodes,
                                      node_color=node_color,
                                      edge_labels=edge_labels, **args)


def construct_frag_tuple(frag):
    """
    Build the fragment tuple given the fragment name.

    Args:
        frag (str): the fragment name

    Returns:
        tuple: the (chain,fragname,fragno) tuple name
    """
    ch, resnum = frag.split('-')
    res, num = resnum.split(':')

    return ch, res, int(num)


def construct_res_tuple(res):
    """
    Build the BigDFT fragment tuple given the residue of the structure

    Args:
       res(Residue): A residue Class on the Biopython package

    Returns:
       tuple: the (chain,fragname,fragno) tuple
    """
    chain = res.full_id[2]
    if len(chain.lstrip(' ')) == 0:
        chain = 'A'
    resname = res.resname
    position = res.full_id[3][1]
    return chain, resname, position


def residue_name(res):
    """
    Name of the residue in the fragment specification of PyBigDFT

    Args:
       res(Residue): A residue Class on the Biopython package
    Returns:
       str: name of the fragment
    """
    chain, resname, position = construct_res_tuple(res)
    return ':'.join(('-'.join((chain, resname)), str(position)))


def biosystem_fragmentation(sys, structure):
    """
    Refragment the system names according to the residues name of the
    structure from the biopython module.

    Args:
       sys (System): and instance of the BioSystem class, fragmented in a
       traditional way
       structure (Bio.Structure): a structure coming from the BioPython module

    Returns:
       System: a system with a fragmentation that is dependent of the
          chain of the molecules
    """
    from BigDFT.Fragments import Fragment
    from BigDFT.Systems import System
    from warnings import warn
    # from Atom import Atom
    # newsys = [at for at in sys.get_posinp()['positions']]
    # units = sys.get_posinp()['units']
    from BigDFT.Atoms import AU_to_A
    # from Fragments import system_from_dict_positions
    mproat = [{at.element: list(at.coord/AU_to_A), 'atom': at,
               'frag': residue_name(at.get_parent()),
               'chain':
               at.full_id[2] if len(at.full_id[2].lstrip(' ')) > 0 else 'A'}
              for at in structure.get_atoms()]
    try:
        lookup = sys.compute_matching(atlist=mproat)
    except ValueError as e:
        lookup = sys.compute_matching(atlist=mproat, check_matching=False)
        warn('Some atoms were not matched' + str(e))
    newsys = System()
    for frag in sys:
        for at, jat in zip(sys[frag], lookup[frag]):
            # reat = mproat[jat]['atom']
            # newname = residue_name(reat.get_parent())
            newname = mproat[jat]['frag']
            if frag not in newname:
                newname = '-'.join([mproat[jat]['chain'], frag])
            # if newname not in newsys:
            #     newsys[newname] = Fragment()
            newsys.setdefault(newname, Fragment()).append(at)
    return newsys


def sequences_to_fragments(chains_to_residues, frag_to_residues):
    """
    Defines the lookup array of the sequences to the fragments

    Args:
       chains_to_residues(list): lookup of the chains to the list of recogized
             residues
       frag_to_residues (list): name of the associated fragments in the order
             of residues list

    Returns:
       list: array of the name of the fragments in the same representation of
           chain_to_residues
    """
    from copy import deepcopy
    lookup = deepcopy(chains_to_residues)
    for c in lookup:
        for i, ires in enumerate(c):
            if ires is None:
                continue
            frag = frag_to_residues[ires]
            c[i] = frag
    return lookup


def _aggregate_values(fragments, data, view, op=sum):
    newdata = [d for d in data]
    for multiple_frag in view.values():
        datasum = op(data[fragments.index(f)] for f in multiple_frag)
        for f in multiple_frag:
            newdata[fragments.index(f)] = datasum
    return newdata


class BioSystem(System):
    """
    Initialize a System for a Biological analysis.
    Retrieve, from the initial file, the basic fragments, the system sequence,
    their mutual interactions and the informations about the purity and the
    observables

    Args:
       filename (str): a pdb file with fragment information inside
       sys (System): an initial system, useful for a custom refragmentation.
            If present, overrides the association written in the pdbfile.
            It is user's responibility to ensure coherency of the total
            number of atoms.
        sequence_from_fragments (bool): method to determine the sequences.
            If True, the sequence are determined by the fragment names.
            It is user's responsibility to provide a inputfile that has a
            set of residues compatible with the fragmentation.
        disable_warnings (bool): disable the BioPython warnings if true
        structure (Bio.Structure): a system structure that is meant to replace
            the one genrated by the file parser. Useful for corrupted pdb
            files. To be used with `sequence from_fragments`
    """
    def __init__(self, filename, sys=None, structure=None,
                 sequence_from_fragments=False, disable_warnings=False):
        from BigDFT.IO import read_pdb
        from os import path
        import warnings
        # disable warnings if asked for
        if disable_warnings:
            from Bio import BiopythonWarning
            warnings.simplefilter('ignore', BiopythonWarning)
        # path of the file
        self.structure_file_path = path.abspath(filename)
        if structure is None:
            # read structure from the Bio.PDB format
            # (can replace custom read_pdb)
            self.structure = read_structure(filename)
        else:
            self.structure = structure

        if sys is None:
            # read pdb from the custom routine
            with open(filename) as ifile:
                sys = read_pdb(ifile)
            # reorganize the fragmentation
            newsys = biosystem_fragmentation(sys, self.structure)
        else:
            newsys = sys

        System.__init__(self, **newsys.dict())

        self._residue_identification(sequence_from_fragments)

        self._sequence_identification(sequence_from_fragments)

        # once the system is identified proceed with the association
        self._fragment_identification()

    def _residue_identification(self, sequence_from_fragments):
        # match the system's fragment into the amminoacids of he structure
        self.frag_to_residues, self.residues = residue_list(
                                        self.structure, self,
                                        strict=not sequence_from_fragments)

    def _sequence_identification(self, sequence_from_fragments):
        # identify sequences
        if sequence_from_fragments:
            chdict = order_chains(self.keys())
            chlett = residue_to_chains(chdict, self.residues)
            self.chains, self.chains_to_residues = chains_to_sequences(chlett)
        else:
            self.chains = pdb_sequences(self.structure_file_path)  # filename)
            # convert the element of the sequences to residues (if possible)
            self.chains_to_residues = sequences_to_residues(
                self.residues, self.chains)
        # associate to each of the chain elements the corresponding fragment
        # names
        self.sequences_to_fragments = sequences_to_fragments(
                              self.chains_to_residues,
                              self.frag_to_residues)

    def _fragment_identification(self):
        # store the fragment keys in a given order
        self.fragment_names = list(self.keys())
        # obtain the list of the letters of each fragment (if possible)
        self.fragment_letters = fragment_letters(
                        self.fragment_names,
                        [seq.sequence for seq in self.chains],
                        self.frag_to_residues,  self.chains_to_residues)
        # identify the fragments that are unmatched either in the residues or
        # in the sequence
        self.unmatched_fragments = [frag for frag, lett in zip(
            self.fragment_names, self.fragment_letters) if lett == 'X' or
            lett == 'None' or lett is None]

    def set_qm_run(self, log):
        """
        Associate the system to a QM run that has been performed by BigDFT.
        Such run assumes that the positions of the system atoms are consistent
        with the PDB

        Args:
           log (Logfiles.Logfile): A Logfile of the BigDFT run
        """
        import numpy as np
        from os import path
        from BigDFT.Logfiles import Logfile
        from BigDFT.Systems import system_from_log
        # path of the logfile
        self.logfile_path = path.abspath(log)
        self.logfile = Logfile(log)
        # TODO: We should insert the compute_matching procedure which checks if
        # the systems are coherent (all the atoms are OK)
        self.set_logfile_info(self.logfile)
        self.fragment_charges = np.array(
            [self[frag].q0[0] for frag in self.fragment_names])
        self.fragment_dipoles = np.array([self[frag].d1()
                                          for frag in self.fragment_names])
        self.fragment_dipole_norms = np.array(
            [np.linalg.norm(d1) for d1 in self.fragment_dipoles])
        atomic_system = system_from_log(self.logfile,
                                        fragmentation='atomic')
        rename_mapping = {'ATOM:'+str(i): i for i in range(len(atomic_system))}
        self.atomic_system = atomic_system.rename_fragments(rename_mapping)
        lookup = self.compute_matching(
            [x[0] for x in self.atomic_system.values()])
        self.atomic_lookup = {x: [z for z in y] for x, y in
                              lookup.items()}
        self._cached_refragment = {}

    def refragment(self, cutoff=0.05, view=None):
        """
        Identify the system's fragmentation and represent it in terms of the
        systems' residues.
        Employs the user defined cutoff to find the correct fragmentation

        Args:
            cutoff (float): cutoff employed for the refragmentation
            view (dict): the view to be imposed as a starting
                fragmentation in the system
        """
        from BigDFT import Systems

        data = self._cached_refragment.get(str(cutoff))
        if data is not None:
            if isinstance(data, (list, tuple)):
                if not isinstance(data[0], Systems.System):
                    pos = data[0]['positions']
                    unt = data[0]['units']
                    data[0] = Systems.system_from_dict_positions(pos,
                                                                 units=unt)
                remapping = {frag: frag.split('+') for frag in data[0]}
            else:
                remapping = data
            return remapping
        fw = Systems.FragmentView(self.purities, self.bond_orders,
                                  charges={k: f.nel for k, f in self.items()})
        if view is None:
            new_view = fw
            sys = Systems.System(self)
        else:
            new_view = fw.refragment(view)
            sys = self.reform_superunits(view)
        remapping = self.bigdft_tool.auto_fragment(sys, new_view, cutoff,
                                                   verbose=True,
                                                   criteria="bondorder")
        # the mapping has to be reworked in terms of a view
        # to prevent renaming of the original fragments
        if view is not None:
            remap = {}
            for frag, lst in remapping.items():
                if frag in view and lst == view[frag]:
                    remap[frag] = lst
                    continue
                remap[frag] = []
                for orig_frag in view:
                    if all(f in frag.split('+') for f in orig_frag.split('+')):
                        remap[frag] += view[orig_frag]
            remapping = remap
        self._cached_refragment[str(cutoff)] = remapping
        return remapping

    def chessboards(self, cutoffs, initial_view=None):
        """
        Define a coherent group of fragmentations that are provided
        from the loosest to tightest cutoffs.
        The views of the fragmentations are defined such as the tighest
        cutoffs always contain regrouping of fragments from the loosest

        Args:
            cutoffs (list): list of floating point numbers of the desired
                cutoffs values
            initial_view (dict): mapping of the initial fragmentation that will
                be imposed strating from the loosest cutoff
        Returns:
            dict: dictionary of the view. the values of the dictionary are
                length-1 lists such that subsequent views (for instance found
                when instantiating a population) can be appended
        """
        coff = sorted(cutoffs, reverse=True)
        cb = {}
        current_view = None if initial_view is None else initial_view.copy()
        self._cached_refragment = {}  # remove the caching for refragmenting
        for c in coff:
            current_view = self.refragment(c, view=current_view)
            cb[c] = [current_view]
        return cb

    @property
    def fragment_ids(self):
        if not hasattr(self, '_fragment_ids'):
            self._fragment_ids = {}
            for frag in self.fragment_names:
                ifrag = 0
                found = False
                for ch in self.sequences_to_fragments:
                    if frag in ch:
                        jfrag = ch.index(frag)
                        found = True
                        break
                    ifrag += len(ch)
                if not found:
                    jfrag = self.unmatched_fragments.index(frag)
                self._fragment_ids[frag] = jfrag + ifrag
        return self._fragment_ids

    def _fragment_labels(self, view=None):
        label = {}
        if view is None:
            view = {f: [f] for f in self.fragment_names}
        for lf in view.values():
            lb = min([self.fragment_ids[f] for f in lf])
            for f in lf:
                if self.fragment_ids[f] == lb:
                    relb = lb
                    for ch in self.chains_to_residues:
                        if relb < len(ch):
                            break
                        # if relb >= len(ch):
                        relb -= len(ch)
                    label[f] = str(relb+1)
                else:
                    label[f] = '   '
        return label

    @property
    def purities(self):
        from BigDFT.PostProcessing import superunits_purities
        if not hasattr(self, '_atomic_purities') and \
           not hasattr(self, '_purities'):
            self._atomic_purities = self.bigdft_tool.run_compute_purity(
                             self.atomic_system, self.logfile,
                             kxs=self.ks_matrix,
                             frag_indices=self.atomic_frag_indices)
        if not hasattr(self, '_purities'):
            atomic_charges = {iat: at.nel for iat, at in
                              self.atomic_system.items()}
            fragment_charges = {frag: f.nel for frag, f in
                                self.items()}
            self._purities = superunits_purities(self.atomic_bond_orders,
                                                 self._atomic_purities,
                                                 atomic_charges,
                                                 self.atomic_lookup,
                                                 fragment_charges)
        return self._purities

    @property
    def fragment_purities(self):
        """
        Define the quality of the fragmentation that is provided in the system
        """
        return [abs(self.purities[frag]) for frag in self.fragment_names]

    @property
    def ks_matrix(self):
        if not hasattr(self, '_kxs'):
            self._kxs = self.bigdft_tool.get_matrix_kxs(self.logfile)
        return self._kxs

    @property
    def sinvh_matrix(self):
        if not hasattr(self, '_sinvh'):
            self._sinvh = self.bigdft_tool.get_matrix_sinvh(self.logfile)
        return self._sinvh

    @property
    def atomic_frag_indices(self):
        if not hasattr(self, '_atomic_fragindices'):
            self._atomic_fragindices = self.bigdft_tool.get_frag_indices(
                self.atomic_system, self.logfile)
        return self._atomic_fragindices

    @property
    def frag_indices(self):
        if not hasattr(self, '_fragindices'):
            self._fragindices = self.bigdft_tool.get_frag_indices(
                self, self.logfile)
        return self._fragindices

    @property
    def bigdft_tool(self):
        from BigDFT.PostProcessing import BigDFTool
        if not hasattr(self, '_btool'):
            from os import environ as env
            if 'BIGDFT_MPIRUN' not in env:
                # to be customized
                env['BIGDFT_MPIRUN'] = 'OMP_NUM_THREADS=2 mpirun -np 2'
            self._btool = BigDFTool()
        return self._btool

    def d3(self):
        from ase.calculators.dftd3 import DFTD3
        from os import environ
        from os.path import join
        from futile.Utils import unique_filename
        if 'ASE_DFTD3_COMMAND' not in environ:
            environ["ASE_DFTD3_COMMAND"] = join(environ["BIGDFT_ROOT"],
                                                "dftd3")
        return DFTD3(xc='pbe', grad=False, damping="bj",
                     label=unique_filename(prefix='ased3_'))

    @property
    def atomic_bond_orders(self):
        if not hasattr(self, '_atomic_BO'):
            atoms = [i for i in range(len(self.atomic_system))]
            self._atomic_BO = self.bigdft_tool.fragment_bond_order(
                       self.atomic_system, atoms, atoms,
                       self.logfile, kxs=self.ks_matrix,
                       frag_indices=self.atomic_frag_indices)
        return self._atomic_BO

    @property
    def bond_orders(self):
        from BigDFT.PostProcessing import superunits_quadratic_quantities
        if not hasattr(self, '_pairwise_BO'):
            self._pairwise_BO = superunits_quadratic_quantities(
                                self.atomic_bond_orders, self.atomic_lookup)
            # self._pairwise_BO = self.bigdft_tool.fragment_bond_order(
            #            self, self.fragment_names, self.fragment_names,
            #            self.logfile, kxs=self.ks_matrix,
            #            frag_indices=self.frag_indices)
        return self._pairwise_BO

    @property
    def atomic_interactions(self):
        if not hasattr(self, '_atomic_interactions'):
            atoms = [i for i in range(len(self.atomic_system))]
            self._atomic_interactions = \
                self.bigdft_tool.fragment_interaction_energy(
                         self.atomic_system, atoms, atoms,
                         self.logfile, kxs=self.ks_matrix,
                         frag_indices=self.atomic_frag_indices,
                         sinvh=self.sinvh_matrix)
        return self._atomic_interactions

    @property
    def interactions(self):
        from BigDFT.PostProcessing import superunits_quadratic_quantities
        if not hasattr(self, '_atomic_interactions') and \
           not hasattr(self, '_interactions'):
            self._interactions = self.bigdft_tool.fragment_interaction_energy(
                     self, self.fragment_names, self.fragment_names,
                     self.logfile, kxs=self.ks_matrix,
                     frag_indices=self.frag_indices,
                     sinvh=self.sinvh_matrix)
        elif not hasattr(self, '_interactions'):
            self._interactions = superunits_quadratic_quantities(
                                 self.atomic_interactions,
                                 self.atomic_lookup)
        return self._interactions

    def fragment_distances(self, target_fragments):
        """
        Calculate the distance between each of the fragment and
        a list of target fragments.

        Args:
            target_fragments (list): Id of the fragments to focus on

        Returns:
            array-like: list of the distances in the order of fragment_names
        """
        from BigDFT.Fragments import pairwise_distance
        from numpy import array

        distances = []
        for f in self.fragment_names:
            distances.append(
                min([pairwise_distance(self[f], self[g])
                     for g in target_fragments]))
        return array(distances)

    def fragment_interaction_strengths(self, target_fragments,
                                       criteria='hamiltonian'):
        """
        Fragment Bond Order between the different fragments of the system,
        possibly specified on target regions.

        Args:
           target_fragments (list): Id of the fragments to focus on
           criteria (str): may be 'hamiltonian', 'electrostatic',
               or 'bond_order'

        Returns:
           array-like: list of the strengths of the interactions of the
                fragments wrt to the target fragments.
                The interaction of the target fragments with themselves is put
                to the maximum value of the other interactions in order
                to avoid scale misbehaviour
        """
        ints = {'bond_order': 'bond_orders',
                'hamiltonian': 'interactions',
                'electrostatic': 'electrostatic_interactions'}
        return interaction_strengths(self.fragment_names, target_fragments,
                                     getattr(self, ints[criteria]))

    def d3PBE_energy(self):
        """
        Full dispersion energy of the system

        Returns:
            float: the d3 PBE energy for the PBE functional, in Hartree
        """
        energy = self.ase_potential_energy(self.d3())
        return energy

    def display_sequences(self, fig=None, labeldict=None, boxcount=None,
                          with_fragment_labels=False, remove_letters=None,
                          **kwargs):
        """
        Display the full sequencies of the system.
        Employ the given field to define the color dict if available

        Args:
           field_vals(list) : values of the field to decide
               the colors of the keys
           labeldict (dict): an optional dictionary mapping fragment ids
               to strings that will be printed on each sequence square.
           colordict (dict) :  the dictionary of the keys,
               and the corresponding colors. Overrides field_vals if present
           fig (matplotlib.pyplot.Figure): figure on which to put the sequence
               into
           with_fragment_labels (bool): a boolean useful to represent the
               fragment numbers. Can be combined with a view to visualize the
               QM fragments. Provides a default labeldict.
           remove_letters (list): list of the letters to be removed from
               the display. In elements of the chains.
           **kwargs: arguments to be passed to display of the colordict
        Returns:
           matplotlib.pyplot.Figure: the figure in which the sequences were
               displayed
        """
        # import matplotlib.pyplot as plt
        colord = self.colordict(**kwargs)
        if with_fragment_labels and labeldict is None:
            labeldict = self._fragment_labels(kwargs.get('view'))

        if labeldict is not None and 'colordict' not in kwargs:
            colord = self._colordict_for_labelling(colord, **kwargs)

        # if fig is None:
        #     fig = plt.figure()
        ifig = 1
        # nrows = len(self.chains)
        for i, (seq, lookup) in enumerate(zip(self.chains,
                                              self.sequences_to_fragments)):
            colorletters = [colord.get(frag) for frag in lookup]
            if remove_letters is None:
                remove_chain_letters = []
            else:
                remove_chain_letters = remove_letters[i]
            if labeldict is not None:
                labelletters = [labeldict.get(frag) for frag in lookup]
            else:
                labelletters = None
            # ax = fig.add_subplot(nrows, 1, ifig)
            # ax =
            seq.display(colors=colorletters, labels=labelletters,
                        boxcount=boxcount, remove_letters=remove_chain_letters)
            ifig += 1
        return fig

    def display_graph(self, restrict_to=None, bo_cutoff=0.01,
                      fragment_labels=None, **kwargs):
        """
        Display a graph view of the system

        Args:
           restrict_to (list): list of fragmens to which the graph
                 has to be restricted
           bo_cutoff (float): the limit for a non-negligible chemical link
           fragment_labels (list, dict): list (in fragment_names order) of the
               labels of the fragments or dictionary of the fragments which
               need relabelling
           **kwargs: colordict args
        """
        labels = self.fragment_letters
        if fragment_labels is not None:
            if isinstance(fragment_labels, dict):
                for ifrag, frag in enumerate(self.fragment_names):
                    if frag in fragment_labels:
                        labels[ifrag] = fragment_labels[frag]
            else:
                labels = fragment_labels
        G = Graph(self.fragment_names, self.bond_orders, labels,
                  restrict_to=restrict_to, cutoff=bo_cutoff)

        colord = self.colordict(**kwargs)

        G.display(colors=[colord[i] for i in self.fragment_names])

    def display(self, cartoon=False, **kwargs):

        System.display(self, colordict=self.colordict(**kwargs),
                       cartoon=cartoon)

    def _colordict_for_labelling(self, colord, **kwargs):
        from futile.Utils import kw_pop
        if 'color_by' not in kwargs and 'field_vals' not in kwargs:
            fragment_colors = colord
            cd = {name: ('white', fragment_colors.get(name, 'white'))
                  for name in colord}
        else:
            new_kw, tt = kw_pop('color_by', None, **kwargs)
            new_kw, tt = kw_pop('field_vals', None, **new_kw)
            fragment_colors = self.colordict(**new_kw)
            cd = {name: (fc, fragment_colors.get(name, 'white'))
                  for name, fc in colord.items()}
        return cd

    def refragmentation_colordict(self, refrag_keys):
        """
        Define a colordict that associates the fragment which go together
        with the same color and put to white all the previously pure fragments.

        Args:
           refrag_keys (list): Keys of the refragmented system

        Returns:
           dict: the dictionary of the colors (equal color means same fragment)
        """
        colordict = self.colordict()
        rekeys = {}
        for frag in refrag_keys:
            # if isinstance(frag, tuple):
            if '+' in frag:
                allfrags = frag.split('+')
                recolor = colordict[allfrags[0]]
                rekeys[allfrags[0]] = []
                for refrag in allfrags[1:]:
                    colordict[refrag] = recolor
                    rekeys[allfrags[0]].append(refrag)
            else:
                colordict[frag] = 'white'
        cc = self.colordict(keys=rekeys.keys())
        for leadfrag in rekeys:
            for refrag in rekeys[leadfrag]:
                cc[refrag] = cc[leadfrag]
        return cc

    def colordict(self, color_by=None, keys=None, field_vals=None,
                  colordict=None, colorcode='html', highlight=None,
                  view=None):
        """
        Define the colordict to quantify the fragments
        according to the method

        Args:
           keys (list): keys of the fragment names
           color_by (str): can be 'charge', 'dipole', 'purity'
           field_vals(list) : values of the field to decide
                 the colors of the keys
           colordict (dict) :  the dictionary of the keys,
                 and the corresponding colors. Overrides field_vals if present
           colorcode (str): the string of the colorcode. can be 'html' or 'rgb'
                if field_dict is absent, otherwise it represents
                the colormap of matplotlib.
                Default is 'seismic' for diverging values
                (i.e. field_dict has negative data), otherwise 'Reds'
            highlight (list) : color in red only the fragments which
                are indicated in the list
            view (dict): the fragmentation of the system to aggregate
                the field_vals to
        """
        from BigDFT import Visualization as V
        if keys is None:
            keys = self.fragment_names
        if colorcode == 'html':
            colorcode = self._colorcode(criteria=color_by)
        if field_vals is None:
            if color_by is not None:
                field_vals = self.fragment_values(criteria=color_by, view=view)
            elif view is not None and colordict is None:
                return self.refragmentation_colordict(view)
        elif view is not None:
            field_vals = _aggregate_values(keys, field_vals, view)
        if highlight is not None:
            colord = {
                frag: 'red' if frag in highlight else 'white' for frag in self}
        elif colordict is None:
            colord = V.get_colordict(keys=keys, field_vals=field_vals,
                                     colorcode=colorcode)
        else:
            colord = colordict
        return colord

    def _colorcode(self, criteria=None):
        """
        Retrieve the colorcode to be used for the numpy cmap
        This is a commodity function made to make more intentional
        some treatment.

        Args:
           criteria (str): can can be 'charge', 'dipole', 'purity'

        Returns:
           str: matplotlib cmap to be used
        """
        if criteria is None:
            return 'html'
        if criteria == 'charge':
            colorcode = 'seismic'
        elif criteria == 'dipole':
            colorcode = 'Greens'
        elif criteria == 'purity':
            colorcode = 'Greys'
        else:
            colorcode = 'html'
        return colorcode

    def fragment_values(self, criteria=None, view=None):
        """
        Retrieve the fragment values that are associated to a give criteria.
        This is a commodity function made to make more intentional some data.

        Args:
           criteria (str): can can be 'charge', 'dipole', 'purity'
           view (dict): the dictionary of the fragment view. If present,
               aggregate the values accordingly.

        Returns:
           list: in order of the fragment names.
        """
        from BigDFT.Systems import FragmentView
        from numpy.linalg import norm
        if criteria is None:
            return None
        field_vals = None
        if criteria == 'charge':
            if view is None:
                field_vals = self.fragment_charges
            else:
                field_vals = _aggregate_values(self.fragment_names,
                                               self.fragment_charges, view)
        elif criteria == 'dipole':
            if view is None:
                field_vals = self.fragment_dipole_norms
            else:
                aggregated_dipoles = _aggregate_values(self.fragment_names,
                                                       self.fragment_dipoles,
                                                       view)
                field_vals = map(norm, aggregated_dipoles)
        elif criteria == 'purity':
            field_vals = self.fragment_purities
            if view is not None:
                fw = FragmentView(self.purities, self.bond_orders,
                                  {f: ff.nel for f, ff in self.items()})
                new_fw = fw.refragment(view)
                for frag, pv in new_fw.purities.items():
                    for f in view[frag]:
                        field_vals[self.fragment_names.index(f)] = abs(pv)
        return field_vals

    def fragment_scatterplot(self, errors=None, **kwargs):
        """
        Defines the scatterplot to have a look to in order to interpret
        sequence data

        Args:
           errors (list): the errors of the field_vals in the order
               of fragment_names
           **kwargs: the same arguments of the `py:meth:colordict` method

        Returns:
           module: reference to the `matplotlib.pyplot` module
        """
        import matplotlib.pyplot as plt
        colors = self.colordict(**kwargs)
        field_vals = kwargs.get('field_vals')
        view = kwargs.get('view')
        if field_vals is None:
            field_vals = self.fragment_values(criteria=kwargs.get('color_by'),
                                              view=view)
        elif view is not None:
            field_vals = _aggregate_values(self.fragment_names, field_vals,
                                           view)
        for i, (f, p) in enumerate(zip(self.fragment_names, field_vals)):
            j = self.fragment_ids[f]
            plt.scatter(j, p, color=colors[f])
        if errors is not None:
            for iname, (y, e) in enumerate(zip(field_vals, errors)):
                x = self.fragment_ids[self.fragment_names[iname]]
                plt.errorbar(x, y, e, fmt='none', elinewidth=0.5, ecolor='k')
        plt.xlabel('Fragment ID')
        return plt

    def represent(self, errors=None, with3d=False, with_fragment_labels=False,
                  boxcount=None, remove_letters=None, **kwargs):
        """
        Represent all the data of the systems, both in with the sequence
        and with the scatterplot. Useful only when there is data to be
        represented.

        Args:
           with3d (bool): if True represent the 3d vision of the system
           errors (list): the errors of the field_vals in the order
               of fragment_names
           **kwargs: the same arguments of the `py:meth:colordict` method

        Returns:
           module: reference to the `matplotlib.pyplot` module
        """
        plt = self.fragment_scatterplot(errors=errors, **kwargs)
        self.display_sequences(with_fragment_labels=with_fragment_labels,
                               boxcount=boxcount,
                               remove_letters=remove_letters,
                               **kwargs)
        if with3d:
            self.display(**kwargs)
        return plt

    def unmatched_frag_info(self):
        """
        Retrieve the information about the fragments which have not been
        matched inside the sequence

        Returns:
           dict: a Dictionary of the fragment that have not been recognized
        """
        from BigDFT.Atoms import AU_to_A
        info = {'expected': {}, 'found': {}}
        expected = []
        for inum, (frag, residue) in enumerate(zip(self.frag_to_residues,
                                                   self.residues)):
            if frag is not None:
                letter = self.fragment_letters[self.fragment_names.index(frag)]
            if frag is not None and letter is not None:
                continue  # fragment has been matched
            centroid = sum([a.coord for a in residue.get_atoms()])/len(residue)
            info['expected'][residue_name(residue)] = {
                 'centroid': list(centroid), 'nat': len(residue),
                 'residue': residue}
            expected += [{at.element: list(at.coord/AU_to_A)}
                         for at in residue.get_atoms()]
        found = System()
        for frag in self.unmatched_fragments:
            info['found'][frag] = {'centroid':
                                   list(self[frag].centroid*AU_to_A),
                                   'nat': len(self[frag]),
                                   'fragment': self[frag]}
            found[frag] = self[frag]
        if len(expected) > 0:
            info['lookup'] = found.compute_matching(atlist=expected)
        return info

    def charge_at_pH(self, ph, closed_shell=True):
        """
        Identify the charge that the system should have for a given pH.
        Employ the
        `~py:meth:Bio.Bio.SeqUtils.ProtParam.ProteinAnalysis.charge_at_ph`
        method of BioPython

        Args:
            ph (float): the value of the ph
            closed_shell (bool): defines if the charge should be rounded
                to the nearest even integer considering the total number
                of electrons of the system.
        Returns:
            float : the value of the charge
        """
        from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
        protein = PA(str(sum(self.chains).sequence))
        charge = protein.charge_at_pH(ph)
        if closed_shell:
            nelec = float(sum(f.nel for f in self.values()))
            total_elec = round(nelec - charge)
            if int(total_elec/2)*2 != int(total_elec) and closed_shell:
                return int(charge+0.5)
        else:
            return round(charge) if closed_shell else charge

    def override_fragment_association(self, fragment_names, new_names=None):
        """
        Force the association of the fragment provided by the fragment_names
        list into the data represented by the new_names.
        New names should be provided as 'Ch-RES:num' where ``Ch`` is the chain,
        ``RES`` is the aminoacid residue at the position ``num``.

        Args:
           fragments_names(list): list of the fragments to be reassociated.
               Should be present in the unmatched_fragment list.
           new_names(list): list of the new names to be produced,
               in fragment_names element order.
               If absent, the fragment names are employed.
        """
        previously_unknown_fragments = {residue_name(res): ires
                                        for ires, (frag, res) in enumerate(
                                        zip(self.frag_to_residues,
                                            self.residues))
                                        if frag is None}
        if new_names is None:
            new_names = fragment_names
        for oldfrag, newfrag in zip(fragment_names, new_names):
            if oldfrag not in self.unmatched_fragments:
                raise ValueError(
                    'The fragment "' + oldfrag +
                    '" is not present in the list of unmatched fragments')
            iold = self.fragment_names.index(oldfrag)
            if self.fragment_letters[iold] != 'X':
                raise ValueError(
                    'The fragment "' + oldfrag +
                    '" has already a letter associated with it')
            if newfrag not in previously_unknown_fragments:
                raise ValueError(
                     'The fragment "' + newfrag +
                     '" cannot be associated to a valid sequence fragment')
            ires = previously_unknown_fragments.pop(newfrag)
            self.frag_to_residues[ires] = oldfrag
            for ic, chain in enumerate(self.chains_to_residues):
                if ires in chain:
                    iseq = chain.index(ires)
                    self.sequences_to_fragments[ic][iseq] = oldfrag
                    break
            self.fragment_letters[iold] = residue_letter(
                    ires, [seq.sequence for seq in self.chains],
                    self.chains_to_residues)
            self.unmatched_fragments.remove(oldfrag)
        if hasattr(self, '_fragment_ids'):
            delattr(self, '_fragment_ids')


def graph_bond(fraglist, threshold, pairwise_bo, restrict_to=None):
    """
    Defines the connectivity matrix of the fragments

    Args:
       fraglist (list): the list of the fragment names
       threshold (float) : the threshold of the cumulative fragment bond order
            that is applied to define the connectivity
       pairwise_bo (dict): list of the fragment BO connections
       restrict_to (list): list of fragments to which the graph has to be
            restricted

    Returns:
       (matrix) : connectivity matrix of the graph
    """
    from numpy import zeros
    mat = zeros((len(fraglist), len(fraglist)))
    lookup = {y: x for x, y in enumerate(fraglist)}

    for i, fragid1 in enumerate(fraglist):
        if restrict_to is not None and fragid1 not in restrict_to:
            continue
        bo = pairwise_bo[fragid1]
        sorted_bo = sorted(bo, key=bo.get, reverse=True)

        remainder = sum(bo.values())
        for key in sorted_bo:
            mat[i, lookup[key]] = 1
            remainder -= bo[key]
            if remainder < threshold:
                break
    return mat


def get_BO_network(connectivity_matrix):
    """
    Defines the connectivity matrix of the fragments
    Args:
       connectivity_matrix (matrix-like): matrix of the connnectivities of
           the system
       pairwise_bo (dict): list of the fragment BO connections

    Returns:
        tuple: graph,target_nodes,target_edges.
             Graph to be plotted and the list of target nodes and edges
             identified by the restriction if present.
    """
    import networkx as NX
    nw = NX.from_numpy_matrix(connectivity_matrix)
    target_nodes = []
    target_edges = []
    for cc in NX.connected_components(nw):
        subgraph = nw.subgraph(cc)
        if len(cc) > 1:
            target_nodes += list(subgraph.nodes)
            target_edges += list(subgraph.edges)
    return nw, target_nodes, target_edges


def fragment_letters(fragnames, sequences, residue_list,
                     sequences_to_residues):
    """
    Give a list of the fragments in terms of their letter

    Args:
        fragnames (list): list of the names of the fragments
        residue_list (list): The fragments of the `system` in the order
              of the residues of the structure
        sequences_to residues (list): strings of the protein sequence
              in the FASTA sequence order
        sequences(list): Sequences of the system
    Returns:
        list: the labels of the fragments as letters. the label is set to 'X'
           if the letter is not recognized
    """
    letters = ['X'] * len(fragnames)
    for ifrag, frag in enumerate(fragnames):
        if frag not in residue_list:
            continue
        ires = residue_list.index(frag)
        letters[ifrag] = residue_letter(ires, sequences, sequences_to_residues)
    return letters


def residue_letter(ires, sequences, sequences_to_residues):
    """
    Find the letter associated with the residue

    Args:
        ires (int): id of the residue
        sequences(list): Sequences of the system
        sequences_to residues (list): strings of the protein sequence
              in the FASTA sequence order

    Returns:
        str: letter of the residue
    """
    for iseq, seq in enumerate(sequences_to_residues):
        if ires in seq:
            ilett = seq.index(ires)
            return sequences[iseq][ilett]


def interaction_strengths(fragments, target_fragments, pairwise_bo):
    """
    Define the interaction strengths of the fragment list with the
    others with respect to a set of target fragments

    Args:
        fragments (list) : the name of the fragments of the
            system in pairwise_bo
        target_fragments (list) : the fragments constituting the target region
        pairwise_bo (dict): the bond order between fragments

    Returns:
        array: the  sum of the weigths of the interactions of the fragments
            in the target region
    """
    import numpy as np
    lookup = [fragments.index(frag) for frag in target_fragments]
    BOtot = np.zeros(len(fragments))
    for ipair in lookup:
        frag_i = fragments[ipair]
        BOtot += np.array(
            [0.5*(pairwise_bo[frag_i][frag] + pairwise_bo[frag][frag_i])
             for frag in fragments])
    if len(lookup) > 0:
        BOtot[np.array(lookup)] = np.nan
    return BOtot


def find_relevant_fragments(sys, data, criteria, columns=None):
    """
    Identify the fragments of the system that fulfill the criteria
    provided by the function in argument

    Args:
        sys (BigDFT.System): the System to identify
        data (array, Dataframe, dict): data to take care of. If ana array,
            it is considered in order of fragment_names.
            If it is a dataframe or a dict, the name of fragments are the keys.
        criteria (func): function that can be applied to the data items.
            The function is applied to each of the data elements and the
            criteria is assumed to be satisfied
            if any of the element fulfill the criteria
        columns (list): the fragment to focus the search to.
            Useful in the case of a dictionary
    Returns:
        list: fragments that fulfill the criteria
    """
    from BigDFT.IO import reorder_fragments
    accepted = []
    for ikey, key in enumerate(sys.fragment_names):
        if isinstance(data, dict):
            if key in data:
                if columns is None:
                    columns = data[key].keys()
                datarr = [data[key][val] for val in columns]
            else:
                datarr = []
        else:  # this assumes that the value is a scalar
            datarr = [data[ikey]]
        condition = any(criteria(d) for d in datarr)
        if condition:
            accepted.append(key)
    return reorder_fragments(accepted)


def _get_system(system=None, archive=None, filename=None,
                logfile=None, deepcopy=False, **kwargs):
    from copy import deepcopy as dc
    if system is not None:
        return dc(system) if deepcopy else system
    if archive is not None:
        sys = load(archive,
                   serialization_version=kwargs.get('serialization_version'),
                   options=kwargs.get('options', {}))
        return sys
    if filename is not None:
        sys = BioSystem(filename, **kwargs.get('options', {}))
    if logfile is not None:
        sys.set_qm_run(logfile)
    return sys


class BioSystemPopulation(BioSystem):
    """
    An ensemble of BioSystem. Useful to provide and analyze averaged quantities

    Args:
       systems(dict, list): BioSystems defining the populations.
           They can be provided as a dict of pdb files plus logfiles,
           or as a dictionary of systems,
           or as a dictionary of archives for serialization.
           The dictionary is composed as follows:
           'label': {'archive': <archive_filename>,
                     'system': <BioSystem instance>,
                     'filename': <pdbfile>,
                     'logfile': <logfile_path>,
                     'weight': weigth of the sample in the population,
                     'mapping': provides the expression of the fragment of
                         the representative in term of the fragments of the
                         item of the population}
            alternatively, a list can be employed, in which case the filling of
            the population is performed sequentially
       representative (str, int): label of the system representative of
           the population
       exclude_representative (bool): if True, the data of the population will
           be evaluated without considering the representative
       fragment_values (list): list of the internal fragment values that can
           be found inside each of the system. Can be left as an empty list
           to delegate  all the evaluation to the populations.
       to_evaluate(dict) : dictionary of functions to be employed to
           evaluate at the creation of a population.
           They must provide a quantity that is related to a system.
       chessboard_cutoffs (list): list of the cutoffs that will have
           to be employed for the fragmentation.
           Use an empty list to deactivate.
       initial_view (dict): view of the superunits fragmentation that will be
           employed to the representative in the lowest cutoff
       **kwargs: arguments to be passed at the loading of the system
    """
    def __new__(cls, systems, representative, **kwargs):
        from futile.Utils import merge_two_dicts
        sys = _get_system(deepcopy=True,
                          **merge_two_dicts(systems[representative], kwargs))
        sys.__class__ = BioSystemPopulation
        return sys

    def __init__(self, systems, representative,
                 fragment_values=['charge', 'dipole', 'purity'],
                 to_evaluate=None,
                 exclude_representative=False,
                 chessboard_cutoffs=[0.05, 0.045, 0.04, 0.035, 0.03, 0.025],
                 initial_view=None,
                 **kwargs):
        from BigDFT.Stats import Population
        self.populations_items = fragment_values
        self.dataframe_populations = ['bond_orders', 'interactions',
                                      'electrostatic_interactions']
        self.populations = {key: Population(labels=self.fragment_names)
                            for key in self.populations_items +
                            self.dataframe_populations}
        self.extra_populations = {}
        self.representative = representative
        self.exclude_representative = exclude_representative
        self._extra_population('system_dfs', None)
        if to_evaluate is not None:
            for label, func in to_evaluate.items():
                self._extra_population(label, func)
        self.systems = systems
        if len(chessboard_cutoffs) > 0:
            self.chessboard_dict = super(
                BioSystemPopulation, self).chessboards(
                cutoffs=chessboard_cutoffs, initial_view=initial_view)
        else:
            self.chessboard_dict = None
        self._fill_populations(**kwargs)

    def _fill_populations(self, **kwargs):
        from futile.Utils import fill_dictionary_in_parallel, kw_pop
        new_kwargs, nthreads = kw_pop('nthreads', 1, **kwargs)
        new_kwargs['biopop'] = self
        new_kwargs['cbs'] = self.chessboard_dict
        if nthreads == 1:
            dict_sys = {}
            if isinstance(self.systems, dict):
                keys = [k for k in self.systems.keys()
                        if not self.exclude_representative or
                        k != self.representative]
            else:
                keys = range(len(self.systems))
                if self.exclude_representative:
                    keys.remove(self.representative)
            for label in keys:
                dict_sys[label] = _data_of_one_system(label, **new_kwargs)
        else:
            dict_sys = fill_dictionary_in_parallel(nthreads,
                                                   self.systems.keys(),
                                                   _data_of_one_system,
                                                   **new_kwargs)
        for results in dict_sys.values():
            for pop_label, pop in self.populations.items():
                pop.append(**results[pop_label])
        self._assign_chessboard()

    def _extra_population(self, label, func):
        from BigDFT.Stats import Population
        self.populations[label] = Population()
        self.populations_items.append(label)
        self.extra_populations[label] = func

    def _assign_chessboard(self):
        if self.chessboard_dict is not None:
            self._cached_refragment = {str(c): clean_view(cb[-1], self)
                                       for c, cb in
                                       self.chessboard_dict.items()}

    def serialize(self, archive, chessboard_dict=None):
        """
        Dump the entire population into an archive, containing
        the representative archive, if it is provided as such
        (otherwise the serialization of the system)
        as well as the numpy files of all the population data.

        Args:
             archive(str): archive to dump the population to
             chessboard_dict (dict): optional dictionary of the fragmentation
                 of the population that has been previously found with
                 a chessboard algorithm. If absent, the internal chessboard is
                 employed
        """
        from BigDFT.Stats import dump_populations
        from futile.Utils import serialize_objects, create_tarball
        if 'archive' in self.systems[self.representative]:
            files = [self.systems[self.representative]['archive']]
        elif 'filename' in self.systems[self.representative]:
            files = [self.systems[self.representative]['filename']]
        else:
            self._cached_refragment = {str(c): cb[-1]
                                       for c, cb in chessboard_dict.items()}
            serialization = BioSystemSerialization(self, version='1.1')
            archive = 'representative.tar.bz2'
            serialization.dump(archive)
            files = [archive]
        files += list(dump_populations(self.populations))
        if chessboard_dict is None:
            chessboard_dict = self.chessboard_dict
        if chessboard_dict is not None and len(chessboard_dict) > 0:
            chb_obj = serialize_objects({'chessboard.json': chessboard_dict})
        else:
            chb_obj = {}
        create_tarball(archive, set(files), chb_obj)

    @classmethod
    def load(cls, archive, **kwargs):
        """
        Create a BioSystemPopulation instance from a serialized archive.

        Args:
            archive (str): the path of the archive
            **kwargs: other arguments that have to be passed to the
             BioSystemPopulation instantiation

        Returns:
            BioSystemPopulation: a biosystem population instance
        """
        from futile.Utils import unpack_tarball
        from shutil import rmtree
        from os.path import join
        from BigDFT.Stats import Population
        import json

        tmpdir, files = unpack_tarball(archive)

        for f in files:
            if '.tar.bz2' in f:
                # we found the system
                sys_dict = [{'archive': join(tmpdir, f)}]
                break
            if '.pdb' in f:
                # we found the filename
                sys_dict = [{'filename': join(tmpdir, f)}]
        loaded = cls(systems=sys_dict, representative=0,
                     fragment_values=[], chessboard_cutoffs=[], **kwargs)
        for f in files:
            if '.json' in f:
                ff = open(join(tmpdir, f), 'r')
                loaded.chessboard_dict = json.load(ff)
                loaded._assign_chessboard()
                break

        pop_data = {f.replace('.npy', ''): Population.load(join(tmpdir, f))
                    for f in files if '.npy' in f}
        loaded.populations = pop_data

        rmtree(tmpdir)
        return loaded

    def _present_as_population(self, method):
        return hasattr(self, 'populations') and method in self.populations \
            and len(self.populations[method].datas) > 0

    @property
    def purities(self):
        if self._present_as_population('purity'):
            pop = self.populations['purity']
            mean = {f: -pop.mean[i] for i, f in enumerate(pop.feature_labels)}
            return mean
        else:
            return super(BioSystemPopulation, self).purities

    @property
    def bond_orders(self):
        if self._present_as_population('bond_orders'):
            pop = self.populations['bond_orders']
            return pop.mean.to_dict()
        else:
            return super(BioSystemPopulation, self).bond_orders

    @property
    def interactions(self):
        if self._present_as_population('interactions'):
            pop = self.populations['interactions']
            return pop.mean.to_dict()
        else:
            return super(BioSystemPopulation, self).interactions

    def fragment_values(self, criteria, view=None):
        if criteria is None:
            return None
        field_vals = self.populations[criteria].mean
        errors = self.populations[criteria].std
        if view is not None:
            field_vals = _aggregate_values(self.fragment_names,
                                           field_vals, view)
        return {'field_vals': field_vals, 'errors': errors}

    def colordict(self, field_vals=None, **kwargs):
        if field_vals is None:
            field_vals = self.fragment_values(criteria=kwargs.get('color_by'))
            if field_vals is not None:
                field_vals = field_vals['field_vals']
        elif isinstance(field_vals, dict) and 'field_vals' in field_vals:
            field_vals = field_vals['field_vals']
        return super(BioSystemPopulation, self).colordict(
            field_vals=field_vals, **kwargs)

    def fragment_scatterplot(self, field_vals=None, errors=None, **kwargs):
        colors = self.colordict(field_vals=field_vals, **kwargs)
        if errors is None and field_vals is None:
            errors = self.fragment_values(
                                    criteria=kwargs.get('color_by'))['errors']
        if field_vals is None:
            field_vals = self.fragment_values(
                                 criteria=kwargs.get('color_by'))['field_vals']
        return super(BioSystemPopulation, self).fragment_scatterplot(
            field_vals=field_vals, errors=errors, colordict=colors, **kwargs)


def _data_of_one_system(label, **kwargs):
    from futile.Utils import merge_two_dicts, write
    biopop = kwargs['biopop']
    args = biopop.systems[label]
    wg = args.get('weight', 1.0)
    verbose = kwargs.get('verbose', True)
    if verbose:
        write(str(label) + ':  # weight='+str(wg))
    sys = _get_system(**merge_two_dicts(args, kwargs))
    results = {}
    # fragment the population
    cbs = kwargs.get('cbs')
    cutoffs = list(cbs.keys()) if cbs is not None else []
    if len(cutoffs) > 0:
        loosest_coff = sorted(cutoffs, reverse=True)[0]
        previous_view = cbs[loosest_coff][-1]
        new_view = clean_view(previous_view, sys)
        sys_cbs = sys.chessboards(cutoffs, initial_view=new_view)
    else:
        sys_cbs = []
    for cutoff in sys_cbs:
        view = cbs[cutoff][-1]
        review = sys_cbs[cutoff][0]
        # review = update_fragmentation(cutoff, view, sys)
        cbs[cutoff].append(review)
        if len(view) != len(review):
            write('The number of fragments has been modified', len(view),
                  len(review))
    for pop_label in biopop.populations:
        of = pop_label
        if verbose:
            write('  - ' + str(pop_label))
        if of in biopop.extra_populations and of != 'system_dfs':
            array = biopop.extra_populations[of](sys)
        elif of in biopop.dataframe_populations:
            array = getattr(sys, of, 0.0)
        else:
            array = sys.fragment_values(of)
        if pop_label == 'system_dfs':
            order = [frag for frag in biopop.fragment_names if frag in sys]
            order += [frag for frag in sys if frag not in order]
            data = sys.to_dataframe(order=order)
        else:
            data = reorder(array, sys, biopop, mapping=kwargs.get('mapping'))
        results[pop_label] = dict(data=data, weight=wg, label=label)
    return results


def update_fragmentation(cutoff, view_old, sys):
    """
    Update the fragmentation views starting from the last element of
    a list of views

    Args:
         cutoff (float): the cutoff to employ for the fragmentation
         view_old (dict): previous view of the system
         sys (BioSystem): the system employed to update the frament list
    Returns:
         dict: new view of the system
    """
    view = clean_view(view_old, sys)
    # erase any previously existing refragmentation
    sys._cached_refragment = {}
    return sys.refragment(cutoff, view=view)


def reorder(array, of, with_respect_to, switch_chains=None, mapping=None):
    """
    Reorder an array of data of systems in such a way that it can be applied
    to characterize another one.

    Args:
        array (array-like, dict): data in the order of `of.fragment_names`.
            in case of a dictionary, it is assumed that it contains the data
            coming from interactions or bond orders.
        of (BioSystem): the system associated to the original data
        with_respect_to (BioSystem): the system to which reorder the data
        switch_chains (dict): mapping of the chains of the old system with
           respect to the new.
        mapping (dict): dictionary of the view that decomposes the fragments of
           ``with_respect_to`` in terms of the fragment of ``of`` in case of
           multiple elements in the values, an aggregation is considered.
           If mapping is present, there should not be undefined fragments
           in their values

    Returns:
        numpy.array: data in the new order. It contains `numpy.nan` for the
            fragments which are not present in the new system.
    """
    from numpy import ones, nan
    from BigDFT.PostProcessing import superunits_quadratic_quantities
    if isinstance(array, dict):
        if mapping is not None:
            newdata = superunits_quadratic_quantities(array, mapping)
        else:
            newdata = array
            # newdata = {}
            # for frag in with_respect_to.fragment_names:
            #     newdata[frag] = {f: array[frag].get(f, 0)
            #                      for f in with_respect_to.fragment_names}
        return newdata
    data = ones(len(with_respect_to.fragment_names))*nan
    if mapping is not None:
        newarray = _aggregate_values(of.fragment_names, array, view=mapping)
    else:
        newarray = array
    for ifrag, frag_orig in enumerate(with_respect_to.fragment_names):
        if switch_chains is None:
            frag = frag_orig
        else:
            frag = None
            for ch_o, ch_n in switch_chains.items():
                if ch_o not in frag_orig:
                    continue
                frag = frag_orig.replace(ch_o+'-', ch_n+'-')
                break
            if frag is None:
                frag = frag_orig
        if frag not in of.fragment_names:
            continue
        jfrag = of.fragment_names.index(frag)
        data[ifrag] = newarray[jfrag]
    return data


def clean_view(view_orig, superunits):
    """
    Define a view of the system that comes from a view of another one.
    Superunits are cleaned in such a way that no inconsistencies are present
    in the new view.

    Args:
       view_orig (dict): the original view we would like to clean
       superunits (list): list of the new system's superunits
    Returns:
       dict: the new view of the system, where combined superunits are
           joined with the '+' symbol
    """
    if view_orig is None:
        return None
    view = view_orig.copy()
    allvals = []
    for k, val_orig in view_orig.items():
        val = [v for v in val_orig]
        for v in val_orig:
            if v not in superunits:
                previous_k = '+'.join(val)
                if k in view:
                    view.pop(k)
                else:
                    view.pop(previous_k)
                val.remove(v)
                rek = '+'.join(val)
                view[rek] = val
            else:
                allvals.append(v)
    for k in list(view.keys()):
        if len(view[k]) == 0:
            view.pop(k)
    for k in superunits:
        if k not in allvals:
            view[k] = [k]
    return view


def rename_residue(res):
    """
    Rename the residue of the system such that alphabetical order
    corresponds to sequence order

    Args:
        res (str): residue of the system, in form "chain-resname:pos"
    Returns:
        str: residue in the form "chain-pos:letter"
    Warning:
        Only put letters to residues which have letter associated
    """
    from Bio.PDB.Polypeptide import three_to_one
    ch, namepos = res.split('-')
    name, pos = namepos.split(':')
    try:
        letter = three_to_one(name)
    except Exception:
        letter = name
    rename = ch+'-'+str(pos).zfill(3)+':'+letter
    return rename


def rename_labels(labels, as_dict=False, from_dict=False,
                  mappable=rename_residue):
    """
    Rename the names of the residues of a BioSystem in order
    to be compatible with other naming schemes.

    Args:
        labels (iterable): contains the list of residue to be renamed
        as_dict (bool): if True, returns a dictionary as a relabel
        from_dict (bool): if True assumes that the remapping has to be
           performed also on the labels values.
        mappable (func): the function that will be applied for renaming
    Returns:
        dict, list: dictionary or list of the renamed labels, depending on
            the arguments.
    """
    relabel = [] if not from_dict else {}
    for label in labels:
        rename = '+'.join([mappable(res) for res in label.split('+')])
        if from_dict:
            relabel[rename] = map(mappable, labels[label])
        else:
            relabel.append(rename)
    return {k: v for k, v in zip(labels, relabel)} if as_dict else relabel
