"""
A module which contains the routines needed for computing the spillage values.
"""

#: Conversion between Atomic Units and Bohr
AU_to_A = 0.52917721092


class MatrixMetadata:
    """
    This class contains the information stored in the sparse matrix
    metadata file.

    Args:
      filename (str): the name of the metadata file to initialize from.

    Attributes:
      atoms (list): a list of atoms with their positions and associated
        basis functions.
      matdim (int): the dimension of the matrix.
    """
    def __init__(self, filename):
        from BigDFT.Atoms import Atom
        self.atoms = []
        symbol_lookup = []
        with open(filename, "r") as ifile:
            # Read the first line
            matinfo = next(ifile).split()
            self.matdim, natoms, ntypes = [int(x) for x in matinfo[:3]]
            # Units
            next(ifile)
            # skip geocode
            line = next(ifile)
            # skip shift?
            line = next(ifile)
            # get the symbol lookup information
            for i in range(0, ntypes):
                line = next(ifile).split()
                symbol_lookup.append(line[2])
            # Get the atom positions
            for i in range(0, natoms):
                line = next(ifile).split()
                adict = {}
                adict["sym"] = symbol_lookup[int(line[0]) - 1]
                adict["r"] = [float(x) for x in line[1:4]]
                # THIS IS BECAUSE THE METADATA FILE CURRENTLY PRINTS THE
                # WRONG UNITS!
                adict["units"] = "bohr"
                adict["indices"] = []
                self.atoms.append(Atom(adict))
            # Get the indices
            for i in range(0, self.matdim):
                line = next(ifile).split()
                self.atoms[int(line[0]) - 1]["indices"].append(i)

    def get_frag_indices(self, system):
        """
        Retrive the indices of the matrix associated with each fragment
        of a system.

        Args:
          system (Systems.System): the set of fragments to get the indices
            of.

        Return:
          (dict): a mapping from fragment to indices.
        """
        matching = system.compute_matching(self.atoms)

        frag_indices = {}
        for fragid in system:
            frag_indices[fragid] = []
            for at in matching[fragid]:
                frag_indices[fragid] += self.atoms[at]["indices"]

        return frag_indices


def _get_csc(filename):
    from scipy.sparse import csc_matrix
    from scipy.io import mmread

    if '.mtx' in filename:
        return csc_matrix(mmread(filename))
    else:
        return read_ccs(filename)


def serial_compute_puritybase(sfile, dfile):
    """
    This routine computes the matrix K*S using python.

    You will need this for computing the purity values.

    Args:
     sfile (str): the file name of the overlap matrix.
         Must be in ccs of mtx format.
     dfile (str): the file name of the density kernel.
         Must be in ccs of mtx format.

    Returns:
      (scipy.sparse.csc, scipy.sparse.csc): K*S
    """
    # Read from file
    smat = _get_csc(sfile)
    dmat = _get_csc(dfile)

    return dmat.dot(smat)


def serial_compute_spillagebase(sfile, hfile):
    """
    This routine computes the matrix S^{-1}*H using python.

    You will need this for computing the spillage values.

    Args:
     sfile (str): the file name of the overlap matrix.
         Must be in ccs of mtx format.
     hfile (str): the file name of the hamiltonian.
         Must be in ccs of mtx format.

    Returns:
      (scipy.sparse.csc, scipy.sparse.csc): Sinverse * H
    """
    # Read from file
    smat = _get_csc(sfile)
    hmat = _get_csc(hfile)

    # Invert
    sinvmat = invert_ccs(smat)

    return sinvmat.dot(hmat)


def invert_ccs(mat):
    """
    Invert a CSC matrix and restore a sparsity pattern such as it can be
    written without the need for large files.

    Args:
       mat (csc_matrix): the sparese matrix to invert

    Returns:
       csc_matrix: the inverse
    """
    from numpy import array, abs
    from scipy.sparse.linalg import inv
    # compute the inverse
    minv = inv(mat)
    # create the mask
    mask = array(abs(minv[minv.nonzero()]) < 1e-8)[0]
    r = minv.nonzero()[0][mask]
    c = minv.nonzero()[1][mask]
    # Filter the matrix.
    minv[r, c] = 0
    minv = 0.5 * (minv + minv.T)
    minv.eliminate_zeros()
    return minv


def read_ccs(fname):
    """
    Read a CCS matrix from file into a scipy csc matrix.

    This routine uses the CCS matrix files that CheSS generates.
    Args:
      fname (str): name of the file.

    Result:
      scipy.sparse.csc_matrix: the matrix in the file.
    """
    from scipy.sparse import csc_matrix

    data = []
    indices = []
    indptr = []

    with open(fname, "r") as ifile:
        # Read the meta data.
        line = next(ifile)
        split = line.split()
        matdim = int(split[0])
        nel = int(split[2])
        split = next(ifile).split()

        # Read in the index pointers.
        indptr = []
        indntx = [int(x) - 1 for x in split]
        indptr += indntx
        while len(indntx) == 10000:
            split = next(ifile).split()
            indntx = [int(x) - 1 for x in split]
            indptr += indntx

        # Read the indices.
        added = 0
        while(added < nel):
            split = next(ifile).split()
            values = [int(x) - 1 for x in split]
            indices.extend(values)
            added += len(values)

        # Read the data
        for line in ifile:
            data.append(float(line))

    matrix = csc_matrix((data, indices, indptr), shape=(matdim, matdim))
    return matrix
