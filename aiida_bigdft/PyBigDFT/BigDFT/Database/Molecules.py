class Molecule(dict):
    """Class to load a molecule from the BigDFT database"""

    def __init__(self, name):
        import os
        from BigDFT.Fragments import System, Fragment
        from BigDFT.XYZ import XYZReader
        # import the positions of the molecules from the XYZ directory
        dirXYZ = os.path.join(os.path.dirname(__file__), 'XYZs')
        filename = os.path.abspath(os.path.join(dirXYZ, name+'.xyz'))
        if not os.path.isfile(filename):
            raise ValueError('Molecule not available')
        ff = XYZReader(filename)
        frag = Fragment(xyzfile=ff)
        system = System({'molecule:0': frag})
        self.update(system.get_posinp())
        # temporary change of the keys 'values' into 'positions'
        if 'values' in self:
            self['positions'] = self.pop('values')
        if 'positions'in self:
            for at in self['positions']:
                if 'r' in at:
                    at.pop('r')
