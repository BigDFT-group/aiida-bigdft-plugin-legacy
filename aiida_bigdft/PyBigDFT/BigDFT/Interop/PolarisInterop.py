"""A module to define typical operations that can be done on biological systems
   starting from the API and the PDB standards employed in the Polaris code

"""


def read_polaris_sle(slefile):
    """
    Read the setup file of polaris with the atomic information.

    Args:
        slefile(str): the sle polaris file.

    Returns:
        list: a list of the atoms of the system with the same order of
            the pdbfile, ready to be passed to the atoms attribute
    """
    attributes = []
    iat = 0
    with open(slefile) as ifile:
        for line in ifile:
            sleline = line.split()
            try:
                zion, ch, pn, mass, chg,\
                    polchg, chi, name, chain, iat, ires, res = sleline
                attributes.append({'q0': [float(chg)], 'name': name,
                                   'resnum': ':'.join(map(str, (res, ires)))})
            except Exception:
                pass
    return attributes


def write_bigdft_charges_in_sle(charges, slefile, outfile):
    """
    Override polaris charges coming from BigDFT file in the sle

    Args:
        charges (dict): dictionary of the {iat: charges}, with iat from the pdb
        slefile (str): the sle polaris file.
        outfile (str): the output file
    """
    fout = open(outfile, 'w')
    with open(slefile) as ifile:
        for line in ifile:
            sleline = line.split()
            try:
                zion, ch, pn, mass, chg,\
                    polchg, chi, name, chain, iat, ires, res = sleline
                newchg = charges[int(iat)]
                fout.write(' '.join(map(str, (zion, ch, pn, mass, newchg,
                                              polchg, chi, name, chain, iat,
                                              ires, res, '\n'))))
            except Exception:
                fout.write(line)
    fout.close()


def read_polaris_pdb(pdbfile, chain_as_letter=False, slefile=None):
    """
    Read coordinates in the PDB format of POLARIS

    Args:
       pdbfile (str): path of the input file
       chain_as_letter (bool): If True, the fifth column
           is assumed to contain a letter
       slefile (str): path of the file ``.sle`` of Polaris from which
           to extract the system's attributes.

    Warning:
       Assumes Free Boundary conditions for the molecule.
       Only accepts atoms that have one letter in the symbol.
       Switch representation if there is a single letter in the fifth column

    Returns:
       System: A system class
    """
    from BigDFT.Fragments import Fragment
    from BigDFT.Systems import System, GetFragId
    from BigDFT.Atoms import Atom
    sys = System()
    units = 'angstroem'
    with open(pdbfile) as ifile:
            for line in ifile:
                if 'ATOM' not in line:
                    continue
                atomline = line.split()
                if chain_as_letter:
                    iat, name, frag, lett, ifrag, x, y, z, sn = atomline[1:10]
                    chain = lett
                    segname = sn
                else:
                    iat, name, frag, ifrag, x, y, z, chlett = atomline[1:9]
                    chain = chlett[2]
                    segname = chlett
                atdict = {str(name[:1]): map(float, [x, y, z]),
                          'frag': [chain+'-'+frag, int(ifrag)], 'name': name,
                          'iat': int(iat), 'segname': segname}
                fragid = GetFragId(atdict, iat)
                if fragid not in sys:
                    sys[fragid] = Fragment()
                sys[fragid].append(Atom(atdict, units=units))
    if slefile is None:
        return sys
    attributes = read_polaris_sle(slefile)
    from BigDFT import Systems as S, Fragments as F, Atoms as A
    system = S.System()
    for name, frag in sys.items():
        refrag = F.Fragment()
        for at in frag:
            atdict = at.dict()
            att = attributes[atdict['iat']-1]
            assert att['name'] == atdict['name']
            atdict.update(att)
            refrag.append(A.Atom(atdict))
        system[name] = refrag
    return system


def split_pdb_trajectory(directory, filename, prefix):
    """
    Split a trajectory file into various PDB files containing a single
    snapshot.

    Args:
        directory (str): the path of the directory in which the trajectory
            file exists
        filename (str): the trajectory file name
        prefix (str): the name to be provided to the resulting files.
           Files are numbered by prefix_0,_1, etc.

    Returns:
        list: list of the files created, including the directory
    """
    files_to_treat = []
    from os.path import join
    f = open(join(directory, filename))
    fileid = 0
    for line in f.readlines():
        if 'CRYST1' in line:
            newfilename = join(directory, prefix + str(fileid) + '.pdb')
            newf = open(newfilename, 'w')
        newf.write(line)
        if 'ENDMDL' in line:
            newf.close()
            files_to_treat.append(newfilename)
            fileid += 1
    f.close()
    return files_to_treat


def convert_pdb(filein, fileout, **kwargs):
    """
    Call the `func:read_polaris_pdb` function to convert the system, and verify
    that the converted filename is giving the same system.

    Args:
        filein (str): input file path
        fileout (str): output file path. Can be identical to filein
        **kwargs: arguments to be passed to `func:read_polaris_pdb`

    Returns:
        tuple: the systems read from filein and fileout, respectively
    """
    from BigDFT.IO import read_pdb, write_pdb
    sys1 = read_polaris_pdb(filein, **kwargs)
    write_pdb(sys1, open(fileout, 'w'))
    # Verify
    sys2 = read_pdb(open(fileout))
    return sys1, sys2
