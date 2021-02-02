"""
This module contains operations for the treatment of quantities related
to particles in space.
"""


def distance_matrix(X):
    """
    Extract the distance matrix of a array of positions with euclidean
    metric.

    Args:
        X (array): array of the positions of the system

    Returns:
        matrix-like: matrix of distances in squareform

    Warning:
        The distance provided here is meaningful only in isolated boundary
        conditions
    """
    from scipy.spatial import distance as d
    compactform = d.pdist(X)
    return d.squareform(compactform)


def inverse_distance_matrix(Rab):
    """
    Calculate the inverse distance matrix, with zeros on the diagonal

    Args:
         Rab (matrix-like): matrix of the atomic distances in squareform
    Returns:
         matrix-like: inverse distance matrix, with diagonal values set to 0
    """
    from numpy import seterr, fill_diagonal
    err = seterr(divide='ignore')
    rabm1 = 1.0/Rab
    seterr(**err)
    fill_diagonal(rabm1, 0.0)
    return rabm1


def electrostatic_energies_matrix(X, Rm1, Q, P):
    """
    Extract the matrix of the electrostatic energies of atoms described
    By their positions, charges and dipoles

    Args:
         X (array): the array of the atom positions
         Q (array): charge monopoles of the atoms (total)
         P (array): array of the dipoles, in Px,Py,Pz convention
         Rm1 (matrix): inverse distance matrix of the positions

    Returns:
         matrix: matrix of the electrostatic energies per atom pair

    Warning:
         This energy is meaningful only for isolated boundary conditions
    """
    from numpy import array, diag
    E1 = Q.dot(Q.T)
    E1 *= Rm1
    if P is not None:
        E3 = P.dot(P.T)
        E3 *= Rm1**3

        PX = P.dot(X.T)
        D = array(diag(PX))
        PX = (D-PX.T).T
        PX *= Rm1**2
        PXQ = PX*Q.T[0]
        PX *= PX.T

        PXQ += PXQ.T
        # PXQ -= 3.0*PX
        PXQ *= Rm1

        E3 = E3 + 3.0*PX*Rm1 - PXQ
    else:
        E3 = 0.0
    return 0.5 * (E1 + E3)


def long_range_energy_matrix(X, Rm1, Z, Q, P):
    """
    Extract the matrix of the electrostatic energies of atoms described
    By their positions, charges and dipoles

    Args:
         X (array): the array of the atom positions
         Q (array): charge monopoles of the electrons
         Z (array): charge monopole of the ions
         P (array): array of the dipoles, in Px,Py,Pz convention
         Rm1 (matrix): inverse distance matrix of the positions

    Returns:
         matrix: matrix of the electrostatic energies per atom pair

    Warning:
         This energy is meaningful only for isolated boundary conditions
    """
    from numpy import array, diag
    QQ = Q + Z
    E1 = QQ.dot(Q.T)
    E1 += E1.T
    E1 *= Rm1
    if P is not None:
        E3 = P.dot(P.T)
        E3 *= Rm1**3

        PX = P.dot(X.T)
        D = array(diag(PX))
        PX = (D-PX.T).T
        PX *= Rm1**2
        PXQ = PX*QQ.T[0]
        PX *= PX.T

        PXQ += PXQ.T
        # PXQ -= 3.0*PX
        PXQ *= Rm1

        E3 = 2.0 * E3 + 3.0*PX*Rm1 - PXQ
    else:
        E3 = 0.0
    return -0.5 * (E1 + E3)


def electrostatic_energy_dict(E, slicing):
    """
    Dictionary of the electrostatic energies starting from a slicing
    defining the blocks of atom to be considered

    Args:
        E (matrix): electrostatic_energies_matrix of the atoms
        slicing (dict): slicing of the fragments, given for example
            by `py:meth:~BigDFT.Systems.System.dataframe_slicing`
    Returns:
        dict: dictionary of the electrostatic interaction per blocks
    """
    eldict = {}
    for frag1 in slicing:
        for frag2 in slicing:
            if frag2 in eldict and frag1 in eldict[frag2]:
                e = eldict[frag2][frag1]
            else:
                s1, e1 = slicing[frag1]
                s2, e2 = slicing[frag2]
                e = E[s1:e1, s2:e2].sum()
                # E[group1][:, group2].sum() useful for non-contiguous
            eldict.setdefault(frag1, {})[frag2] = e
    return eldict


class PointParticles():
    """
    An object storing the point particle information

    Args:
        X (array): array of atomic positions
        Q (array): charge monopoles of the electrons
        Z (array): charge monopole of the ions
        P (array): electrostatic dipoles
    """
    def __init__(self, X, Z=None, Q=None, P=None):
        self.X = X
        if Z is not None:
            self.Z = Z
        if Q is not None:
            self.Q = Q
        self.P = P

    @property
    def R(self):
        """
        Distance matrix between point PointParticles
        """
        if not hasattr(self, '_R'):
            self._R = distance_matrix(self.X)
        return self._R

    @property
    def Rm1(self):
        """
        Inverse Distance matrix
        """
        if not hasattr(self, '_Rm1'):
            self._Rm1 = inverse_distance_matrix(self.R)
        return self._Rm1

    @property
    def Eel(self):
        """
        Matrix of the pair atomic electrostatic energies
        """
        if not hasattr(self, '_Eel'):
            self._Eel = electrostatic_energies_matrix(self.X, self.Rm1,
                                                      self.Q + self.Z, self.P)
        return self._Eel

    @property
    def Elr(self):
        """
        Matrix of the long range energy
        """
        if not hasattr(self, '_Elr'):
            self._Elr = long_range_energy_matrix(self.X, self.Rm1,
                                                 self.Z, self.Q, self.P)
        return self._Elr

    def Eel_dict(self, slicing):
        """
        Dictionary of the electrostatic energies given a slicing defining
        the groups of atoms to gather
        """
        return electrostatic_energy_dict(self.Eel, slicing)

    def Elr_dict(self, slicing):
        """
        Dictionary of the electrostatic energies given a slicing defining
        the groups of atoms to gather
        """
        return electrostatic_energy_dict(self.Elr, slicing)
