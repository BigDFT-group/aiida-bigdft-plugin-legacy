"""
A module for handling unit cells for periodic calculations.
"""


class UnitCell:
    """
    Defines a wrapper for unit cells.

    Args:
        cell (list): a list of unit cell vectors. This can either be a 3x1 list
          or a 3x3 list. Currently, only orthorhombic cells are supported.
          units (str): the units of the cell parameters. If cell is set to
          None, the free boundary condition is enforced.
        units (str): the unit of length.
    """
    def __init__(self, cell=None, units="bohr"):
        from BigDFT.Atoms import AU_to_A, IsAngstroem

        # Early exit for free boundary condition
        if cell is None:
            self.cell = [[float("inf"), 0, 0],
                         [0, float("inf"), 0],
                         [0, 0, float("inf")]]
            return

        # Copy cell based on the format
        if isinstance(cell[0], list):
            self.cell = cell
        else:
            self.cell = [[cell[0], 0, 0], [0, cell[1], 0], [0, 0, cell[2]]]

        # Check that the unit cell is valid
        if self[2, 2] == float("inf") and not (self[1, 1] == float("inf") and
           self[0, 0] == float("inf")):
            raise ValueError("Boundary condition: ",
                             "If z is infinity, x and y must also be.")
        elif self[0, 0] == float("inf") and self[1, 1] != float("inf"):
            raise ValueError("Boundary condition: ",
                             "If x is infinity, y must also be.")
        elif self[0, 1] != 0 or self[0, 2] != 0 or \
                self[1, 0] != 0 or self[1, 2] != 0 or \
                self[2, 0] != 0 or self[2, 1] != 0:
            raise ValueError("Non-orthorhombic cells not supported")

        for i in range(3):
            for j in range(3):
                if self[i, j] < 0:
                    raise ValueError("Unit cell must be non-negative values.")

        # Store internally in bohr
        if IsAngstroem(units):
            for i in range(3):
                for j in range(3):
                    self.cell[i][j] /= AU_to_A

    def __getitem__(self, idx):
        return self.cell[idx[0]][idx[1]]

    def get_boundary_condition(self, units="bohr"):
        """
        Get a string description of the boundary condition (i.e. free,
        surface, wire, periodic)
        """
        from BigDFT.Atoms import AU_to_A, IsAngstroem

        # Match units
        if IsAngstroem(units):
            conv = AU_to_A
        else:
            conv = 1

        if all([self[i, i] == float("inf") for i in range(3)]):
            return "free"
        elif self[0, 0] == float("inf") and self[1, 1] == float("inf"):
            return "wire 0.0 0.0 " + str(self[2, 2]*conv)
        elif self[1, 1] == float("inf"):
            return "surface " + str(self[0, 0]*conv) + " 0.0 " + \
                   str(self[2, 2]*conv)
        else:
            return "periodic " + " ".join([str(conv*self[i, i])
                                           for i in range(3)])

    def get_posinp(self, units="bohr"):
        """
        Create the dictionary representation of the cell that is passed to
        BigDFT.
        """
        from BigDFT.Atoms import AU_to_A, IsAngstroem
        if IsAngstroem(units):
            return [self[i, i] * AU_to_A for i in range(3)]
        else:
            return [self[i, i] for i in range(3)]

    def minimum_image(self, pos, units="bohr"):
        """
        Given a vector of three positions, this wraps those positions inside
        the cell using the minimum image convention.
        """
        from BigDFT.Atoms import AU_to_A, IsAngstroem

        # Match units
        if IsAngstroem(units):
            conversion = 1/AU_to_A
        else:
            conversion = 1
        bohrpos = [x*conversion for x in pos]

        # Adjust position
        for i in range(3):
            if self[i, i] == float("inf"):
                continue
            while(bohrpos[i] > self[i, i]):
                bohrpos[i] -= self[i, i]
            while(bohrpos[i] < 0):
                bohrpos[i] += self[i, i]

        # Convert back
        return [x/conversion for x in bohrpos]


def _example():
    # Create a basic unit cell
    cell = UnitCell([10, 8, 4], units="angstroem")

    # Print out the posinp representation
    print(cell.get_posinp())
    print(cell.get_posinp(units="angstroem"))

    # Right now we enforce the orthorhombic condition
    try:
        cell = UnitCell([[10, 0, 0], [0, 8, 0], [0, 4, 10]])
    except ValueError as e:
        print(e)
    cell = UnitCell([[10, 0, 0], [0, 8, 0], [0, 0, 4]])
    print(cell.get_boundary_condition("angstroem"))

    # Wire boundary condition
    wire = UnitCell([float("inf"), float("inf"), 4])
    print(wire.get_posinp())
    print(wire.get_boundary_condition())

    # Surface boundary condition
    surface = UnitCell([10, float("inf"), 4])
    print(surface.get_posinp())
    print(surface.get_boundary_condition())

    # Wrap positions to the minimum image convention.
    pos = [-5, -2, -3]
    print(cell.minimum_image(pos))
    print(wire.minimum_image(pos))
    print(surface.minimum_image(pos))

    pos = [15, 12, 13]
    print(cell.minimum_image(pos))
    print(wire.minimum_image(pos))
    print(surface.minimum_image(pos))


if __name__ == "__main__":
    _example()
