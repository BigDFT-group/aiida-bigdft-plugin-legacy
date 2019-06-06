class Molecule(dict):
    def __init__(self,name):
        import os
        from BigDFT.Fragments import System
        #import the positions of the molecules from the XYZ directory
        dirXYZ = os.path.join(os.path.dirname(__file__),'XYZs')
        filename = os.path.abspath(os.path.join(dirXYZ,name+'.xyz'))
        if not os.path.isfile(filename):
            raise ValueError('Molecule not available')
        system = System(xyz=filename)
        self.update(system.dict())
        #temporary change of the keys 'values' into 'positions'
        self['positions']=self.pop('values')
        for at in self['positions']:
            if 'r' in at: at.pop('r')
