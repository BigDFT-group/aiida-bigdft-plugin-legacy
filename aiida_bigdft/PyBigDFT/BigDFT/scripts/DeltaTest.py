from futile.Utils import write as safe_print
from . import InputGenerator

cpickle_file = "bigdft-results.cpickle"


def process_element(run, var, dico, force_run, Elements, strains):
    import math
    import cPickle
    name = dico["name"]
    hgrids = var['dft']['hgrids']
    ngrids = (math.ceil(strains[0] * dico["a"] / hgrids[0]),
              math.ceil(strains[0] * dico["b"] / hgrids[1]),
              math.ceil(strains[0] * dico["c"] / hgrids[2]))
    for strain in strains:
        if run.iproc == 0:
            safe_print("(%s)" % strain)
        if strain in dico["eKS"] and not force_run:
            # The calculation is already done
            if run.iproc == 0:
                safe_print((name, strain), ": ", dico["eKS"][strain])
            continue
        var.update(InputGenerator.set_strain(strain, ngrids, dico))

        run.update(var)

        # Run and store the results
        out = run.run()

        var.update(InputGenerator.set_restart())

        run.update(var)
        # And restart
        out = run.run()

        Elements[name] = dico
        dico["eKS"][strain] = out.eKS
        dico["FKS"][strain] = out.etot
        dico["pressure"][strain] = out.pressure
        if run.iproc == 0:
            cPickle.dump(Elements, open(cpickle_file, "w"))
        # Parse the generated log to get the number of grid points.
        for l in open("log-%s.yaml" % var["radical"], "r"):
            if "Grid Spacing Units" in l:
                grid = l.split()
        ngrids = (int(grid[5][:-1]), int(grid[6][:-1]), int(grid[7]))
        ngrids = [a + 1 for a in ngrids]
    # Freeing memory
    run = None
    out = None
    # Build the file
    if run.iproc == 0:
        fd = open("%s.dat" % name, 'w')
        for strain in strains:
            dico = Elements[name]
            # Volume in A^3/atom
            volume = dico["volume"]*strain**3/dico["nat"]
            HatoeV = 27.21138386
            eKS = float(dico["eKS"][strain])*HatoeV/dico["nat"]
            fd.write("%16.9f %16.9f\n" % (volume, eKS))
            safe_print(strain, dico["eKS"][strain], volume, eKS)
        fd.close()


# Build a dictionary


def elements_from_cif(files, name):
    # filename match as Unix shell
    import fnmatch
    import os
    Elements = dict()
    nonortho = list()
    ortho = list()
    for file in files:
        if not fnmatch.fnmatch(file, "*.cif"):
            continue
        dico = dict()
        dico["file"] = file
        # open("%s/%s" %(dirCIF,file)).readlines():
        for line in open(os.path.abspath(file)).readlines():
            if "_" in line:
                # Remove quote '
                items = line.replace("'", "").split()
                if len(items) == 2:
                    dico[items[0]] = items[1]
            elif "Biso" in line or name in line:
                items = line.split()
                # We know that the maximum of atoms is 8 in the cell
                n = int(items[0][-1])
                dico['nat'] = n
                # dico[n] = "%s  %s  %s " % (items[2], items[3], items[4])
                dico[n] = list(map(float, items[2:5]))
        # Convert name
        if name is not None:
            dico["name"] = name
        else:
            dico["name"] = dico["_pd_phase_name"]
        # We use bohrs
        atob = 1.0/0.5291772108
        dico["a"] = float(dico["_cell_length_a"])*atob
        dico["b"] = float(dico["_cell_length_b"])*atob
        dico["c"] = float(dico["_cell_length_c"])*atob
        dico["alpha"] = float(dico["_cell_angle_alpha"])
        dico["beta"] = float(dico["_cell_angle_beta"])
        dico["gamma"] = float(dico["_cell_angle_gamma"])
        # Only for non-orthorhombic and in angstroem^3
        dico["volume"] = dico["a"]*dico["b"]*dico["c"]/atob**3
        # Create a key results
        dico["eKS"] = dict()
        dico["FKS"] = dict()
        dico["pressure"] = dict()

        InputGenerator.transform_to_orthorombic(dico)

        name = dico['name']
        # Update volume after orthomrombic tranformation and in angstroem^3
        dico["volume"] = dico["a"]*dico["b"]*dico["c"]/atob**3
        if (dico["alpha"] != 90 or dico["beta"] != 90 or
                dico["gamma"] != 90):
            nonortho.append(name)
        else:
            ortho.append(name)
        # Add in the large dictionary
        Elements[name] = dico
    ortho.sort()
    nonortho.sort()
    return Elements, ortho, nonortho

def xyz_from_elements(dirXYZ, Elements, ortho, nonortho):
    import shutil
    import os
    import time
    format_xyz = """{0[nat]} reduced
 periodic {0[a]}   {0[b]}   {0[c]} 
"""
    # print "Remove the directory '%s'." % dirXYZ
    shutil.rmtree(dirXYZ, ignore_errors=True)
    # Create it
    os.mkdir(dirXYZ)
    safe_print("---")
    # Start and create the xyz files
    safe_print("Delta-test timestamp:", time.strftime('%X %x %Z'))
    safe_print("Delta-test code: BigDFT")
    safe_print("Number of elements: ", len(Elements))
    safe_print("List of elements:", list(Elements.keys()))
    safe_print("Number of orthorhombic elements: ", len(ortho))
    safe_print("Orthorhombic elements: ", ortho)
    safe_print("Number of non-orthorhombic elements: ", len(nonortho))
    safe_print("Non-orthorhombic elements: ", nonortho)
    for dico in Elements.values():
        name = dico['name']
        # We have all the specification
        fnew = os.path.join(dirXYZ, name+".xyz")  # "%s/%s.xyz" % (dirXYZ,name)
        fd = open(fnew, "w")
        fd.write(format_xyz.format(dico))
        for i in range(dico['nat']):
            fd.write("%s " % name)
            print(i, dico['nat'])
            fd.write("%f %f %f\n" % tuple(dico[i+1]))
            # fd.write("%s %s\n" % (name,dico[i+1]))
        fd.close()
        # print "#Creation of the file '{0:s}' from '{1:s}/{2:s}'".format(
        #       fnew,dirCIF,dico['file'])

def xyz_from_pymatgendict(dirXYZ, name, dico):
    import shutil
    import os
    import time
    format_xyz = """{1} reduced
 periodic {0[a]}   {0[b]}   {0[c]} 
"""
    #danger zone.. don't do ?
    if dirXYZ != ".":
        print ("Remove the directory '%s'." % dirXYZ)
        shutil.rmtree(dirXYZ, ignore_errors=True)
    # Create it
    os.mkdir(dirXYZ)
    safe_print("---")
    # Start and create the xyz file
    # We have all the specification
    fnew = os.path.join(dirXYZ, name+".xyz")  # "%s/%s.xyz" % (dirXYZ,name)
    fd = open(fnew, "w")
    nat = len(dico['sites'])
    fd.write(format_xyz.format(dico['lattice'],nat))
    for atom in dico['sites']:
        fd.write("%s " % atom['label'])
        fd.write("%f %f %f\n" % tuple(atom['xyz']))
    fd.close()
        # print "#Creation of the file '{0:s}' from '{1:s}/{2:s}'".format(
        #       fnew,dirCIF,dico['file'])

class Benchmark():
    def __init__(self, calculator):
        import os
        # import the positions of the periodic table from the CIFs directory
        dirCIF = os.path.join(os.path.dirname(__file__), 'CIFs')
        files = [os.path.abspath(os.path.join(dirCIF, f))
                 for f in os.listdir(dirCIF)]
        self.Elements, self.ortho, self.nonortho = elements_from_cif(files)
        # Strains (The original)
        self.strains = [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06]
        # From Francois Jollet
        # strains = [ 0.9795861087155615, 0.986484829732188,
        #             0.9932883883792687,  1.00, 1.006622709560113,
        #             1.0131594038201772, 1.0196128224222163 ]
        # ortho = [ "Ag"]
        self.strains.reverse()
        # strains = [ 1.00 ]
        self.calculator = calculator
        # Now build "*.xyz" for each elements
        if self.calculator.iproc == 0 and False:
            xyz_from_elements("XYZs", self.Elements, self.ortho, self.nonortho)

    def run(self, atoms=None, force=False):
        import os
        # To dump python object in a file
        import cPickle
        if atoms is None:
            atomlist = self.ortho
        else:
            atomlist = atoms
        # Loop over the elements
        # Test if the file bigdft-results.cpickle exists and use it instead
        if os.path.exists(cpickle_file) and not force:
            self.Elements = cPickle.load(open(cpickle_file, "r"))
        for name in atomlist:
            dico = self.Elements[name]
            if self.calculator.iproc == 0:
                safe_print("Start calculation of %s" % dico["name"])
            hgrid = 0.30
            var = InputGenerator.set_inputfile(hgrid, dico)
            self.calculator.set(var)
            process_element(self.calculator, var, dico, force,
                            self.Elements, self.strains)
