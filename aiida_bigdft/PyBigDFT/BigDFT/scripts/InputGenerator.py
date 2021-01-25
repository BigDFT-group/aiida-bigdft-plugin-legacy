from futile.Utils import write as safe_print
import io

def set_inputfile(hgrid, dico, init_input=None, psp=[], units="reduced"):
    basicinput = """
logfile: Yes
dft:
    ixc: PBE
    ncong: 2
    rmult: [10, 8]
    itermax: 3
    idsx: 0
    gnrm_cv: 1e-8
#Control the diagonalisation scheme
mix:
    iscf: 7
    itrpmax: 200
    rpnrm_cv: 1e-12
    tel: 1e-3
    alphamix: 0.5
    norbsempty: 1000
    alphadiis: 1.0

#perf:
#    accel: OCLGPU
#    ocl_devices: Tesla K40c
#    blas: Yes

"""
    import yaml
    if not init_input:
        var = yaml.load(basicinput)
    else:
        var = init_input
    # Spin parameters
    var["dft"].update(set_spin(dico["name"], dico["nat"]))
    # K point parameters
    var["kpt"] = set_kpoints(dico["nat"])
    var["dft"]["hgrids"] = (hgrid, hgrid, hgrid)
    # Control the diagonalisation scheme
    if dico["name"] in ("Cr", ):
        var["mix"]["iscf"] = 3
        var["mix"]["alphamix"] = 0.9
    if dico["name"] in ("Ba", "Ca"):
        var["mix"]["norbsempty"] = 8
    var["ig_occupation"] = {dico["name"]: {"empty_shells": ("s", "p", "d")}}
    set_psp(dico["name"], psp)
    if units == "reduced":
        positions = [{dico["name"]: dico[i + 1]} for i in range(
            dico["nat"])]
    else:
        positions = [{dico["name"]: [dico[i + 1][0]*dico["a"], dico[i + 1][1]*dico["b"],dico[i + 1][2]*dico["c"]]} for i in range(
            dico["nat"])]

    var["posinp"] = {"positions": positions, "units": units,
        "cell": (dico["a"], dico["b"], dico["c"])}
    # We round robin the igspins.
    if "mpol" in var["dft"]:
        mpol = 0
        while mpol < var["dft"]["mpol"]:
            for at in var["posinp"]["positions"]:
                if mpol < var["dft"]["mpol"]:
                    if "IGSpin" in at:
                        at["IGSpin"] += 1
                    else:
                        at["IGSpin"] = 1
                    mpol += 1
    elif "nspin" in var["dft"] and var["dft"]["nspin"] == 2:
        for (i, at) in enumerate(var["posinp"]["positions"]):
            at["IGSpin"] = 1 - 2 * (i % 2)
    return var

def set_psp(name, psp):
        # Atoms
    import os
    pspfile = "psppar."+name
    dirname = os.path.dirname(__file__)
    filename =  os.path.join(dirname, "psppar", pspfile)
    if not os.path.isfile(filename):
        safe_print("WARNING: Using default PSP for atom", filename, name)
    else:
        psp.append(filename)

def set_spin(name, nat):
    "Define the spin in function of the nature of the atoms"
    dspin = {}
    if name == 'O':
        dspin["nspin"] = 2
    elif name == 'Cr' or name == 'Mn':
        dspin["nspin"] = 2
    elif name == 'Fe' or name == 'Co' or name == 'Ni':
        dspin["nspin"] = 2
        mpol = {"Fe": 2.22, "Co": 1.72, "Ni": 0.60}
        dspin["mpol"] = int(mpol[name] * nat)
        if dspin["mpol"] % 2 == 1:
            dspin["mpol"] += 1
    else:
        dspin["nspin"] = 1
    return dspin


def set_kpoints(nat):
    "Define the k point mesh"
    dkpt = {}
    dkpt["method"] = "mpgrid"
    if nat == 1:
        dkpt["ngkpt"] = [19, 19, 19]
    elif nat == 2:
        dkpt["ngkpt"] = [15, 15, 15]
    elif nat == 3:
        dkpt["ngkpt"] = [14, 14, 14]
    elif nat == 4:
        dkpt["ngkpt"] = [12, 12, 12]
    elif nat == 6:
        dkpt["ngkpt"] = [11, 11, 11]
    elif nat == 8:
        dkpt["ngkpt"] = [10, 10, 10]
    elif nat == 12:
        dkpt["ngkpt"] = [9, 9, 9]
    return dkpt


def set_strain(strain, ngrids, dico):
    return {"radical": dico["name"] + "-%s" % strain,
            "posinp": {"cell": (strain * dico["a"],
                                strain * dico["b"],
                                strain * dico["c"])},
            "dft": {"hgrids": ((strain * dico["a"] + 1e-4) / ngrids[0],
                               (strain * dico["b"] + 1e-4) / ngrids[1],
                               (strain * dico["c"] + 1e-4) / ngrids[2])}}


def set_restart():
    return {"mix": {"tel": 1.e-5}, "dft": {"inputpsiid": 1}}


def transform_to_orthorombic(dico):
    import math
    btype = None
    # Look for possible orthorombic transformation :
    if (dico["alpha"] != 90 and dico["b"] == dico["c"] and
            dico["beta"] == 90 and dico["gamma"] == 90):
        btype = "hexagonal"
        la = "b"
        lb = "c"
        ang = "alpha"
        P = ((1.,  0.,  0.),
             (0.,  0.5, 0.5),
             (0., -0.5, 0.5))
        du = 0.
        dv = 0.
        dw = 1.
    elif (dico["beta"] != 90 and dico["a"] == dico["c"] and
          dico["alpha"] == 90 and dico["gamma"] == 90):
        la = "a"
        lb = "c"
        btype = "hexagonal"
        ang = "beta"
        P = ((0.5, 0., 0.5),
             (0.,  1., 0.),
             (-0.5, 0., 0.5))
        du = 0.
        dv = 0.
        dw = 1.
    elif (dico["gamma"] != 90 and dico["a"] == dico["b"] and
          dico["alpha"] == 90 and dico["beta"] == 90):
        la = "a"
        lb = "b"
        btype = "hexagonal"
        ang = "gamma"
        P = ((0.5, 0.5, 0.),
             (-0.5, 0.5, 0.),
             (0.,  0.,  1.))
        du = 0.
        dv = 1.
        dw = 0.
    elif (dico["gamma"] == 90 and dico["alpha"] == 90 and
          dico["beta"] == 90):
        btype = "orthorombic"
    elif (dico["a"] == dico["b"] and dico["a"] == dico["c"] and
          dico["alpha"] == dico["beta"]and dico["alpha"] == dico["gamma"]):
        btype = "rhombohedral"
    # Transform to orthorombic when possible.
    if btype == "hexagonal":
        a = dico[la]
        b = dico[lb]
        # - hexagonal case.
        alpha = float(dico[ang]) * math.pi / 180.
        dico[la] = math.sqrt((a+b*math.cos(alpha)) **
                             2 + (b*math.sin(alpha))**2)
        dico[lb] = math.sqrt((a-b*math.cos(alpha)) **
                             2 + (b*math.sin(alpha))**2)
        for i in range(dico["nat"]):
            u, v, w = dico[i+1]
            a = P[0][0] * u + P[0][1] * v + P[0][2] * w
            b = P[1][0] * u + P[1][1] * v + P[1][2] * w
            c = P[2][0] * u + P[2][1] * v + P[2][2] * w
            dico[i + 1] = (a, b, c)
            u += du
            v += dv
            w += dw
            a = P[0][0] * u + P[0][1] * v + P[0][2] * w
            b = P[1][0] * u + P[1][1] * v + P[1][2] * w
            c = P[2][0] * u + P[2][1] * v + P[2][2] * w
            dico[dico["nat"] + i + 1] = (a, b, c)
        dico["nat"] *= 2
        dico[ang] = "90"
    elif btype == "rhombohedral" and dico["alpha"] == 60:
        a = dico["a"]
        dico["a"] = a * math.sqrt(2.)
        dico["b"] = a * math.sqrt(2.)
        dico["c"] = a * math.sqrt(2.)
        P = ((0.,  0.5, 0.5),
             (0.5, 0.,  0.5),
             (0.5, 0.5, 0.))
        dd = ((0., 0., 0.),
              (1., 0., 0.),
              (0., 1., 0.),
              (0., 0., 1.))
        nat = dico["nat"]
        dico["nat"] *= len(dd)
        for i in range(nat):
            for (j, (du, dv, dw)) in enumerate(dd):
                u, v, w = dico[i+1]
                u += du
                v += dv
                w += dw
                a = P[0][0] * u + P[0][1] * v + P[0][2] * w
                b = P[1][0] * u + P[1][1] * v + P[1][2] * w
                c = P[2][0] * u + P[2][1] * v + P[2][2] * w
                dico[j * nat + i + 1] = (a, b, c)
                dico["alpha"] = 90.0
                dico["beta"] = 90.0
                dico["gamma"] = 90.0
    elif btype != "orthorombic":
        safe_print("to be treated", dico["name"], dico["alpha"],
                   dico["beta"], dico["gamma"], dico["a"], dico["b"],
                   dico["c"])
