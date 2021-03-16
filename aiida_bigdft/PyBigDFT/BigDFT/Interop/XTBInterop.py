"""
This module contains some wrappers for using XTB to perform calculations.

https://xtb-docs.readthedocs.io/en/latest/contents.html
"""
from BigDFT.Calculators import Runner
from futile.Utils import write as safe_print


class XTBCalculator(Runner):
    """
    A calculator which can drive simulations using the XTB program.

    XTB uses a mixed approach for control. First, you can pass command line
    arguments to the program. Second, you can modify the input file. We
    mirror this split here. Command line arguments can be passed as arguments
    to the run command. To modify the input file, the run command accepts
    an argument `inp` which is a dictionary. This dictionary is written
    to the input file (in the appropriate format). See the example for
    more details.

    https://xtb-docs.readthedocs.io/en/latest/contents.html
    """
    import os

    def __init__(self, omp=os.environ.get('OMP_NUM_THREADS', '1'),
                 dry_run=False, skip=False, verbose=True):
        # Use the initialization from the Runner class (so all options inside
        # __global_options)
        Runner.__init__(self, omp=str(omp), dry_run=dry_run, skip=skip,
                        verbose=verbose)
        # Build the command setting the number of omp threads
        self.command = 'xtb'
        safe_print(
            'Initialize a Calculator with OMP_NUM_THREADS=%s and command %s' %
            (self._global_options['omp'], self.command))

    def pre_processing(self):
        """
        Process local run dictionary to create the input directory and identify
        the command to be passed

        Returns:
            :py:class:`dict`: dictionary containing the command to be passed to
            :meth:`process_run`
        """
        self._ensure_run_directory()
        sys = self.run_options["sys"]
        inp = self.run_options.get("inp", {})

        # Write out the geometry file
        with open(self._get_inputfilename(), "w") as ofile:
            if sys.cell.get_boundary_condition() != "free":
                if any([sys.cell[x, x] == float("inf") for x in range(3)]):
                    raise ValueError("Only 3D periodic or free boundary " +
                                     "condition is supported")
                ofile.write("$periodic 3\n")
                ofile.write("$lattice bohr \n")
                for i in range(3):
                    ofile.write("  ")
                    ofile.write(" ".join([str(sys.cell[i, j])
                                          for j in range(3)]))
                    ofile.write("\n")

            ofile.write("$coord\n")
            for fragid, frag in sys.items():
                for at in frag:
                    pos = at.get_position("bohr")
                    ofile.write("  ")
                    ofile.write(str(pos[0]) + " ")
                    ofile.write(str(pos[1]) + " ")
                    ofile.write(str(pos[2]) + " ")
                    ofile.write(at.sym.lower() + "\n")
            ofile.write("$end\n")

            constr = self.run_options.get("constraints", None)
            if constr:
                self._set_constraints(constr, ofile)

            for k1, v1 in inp.items():
                ofile.write("$"+k1+"\n")
                for k2, v2 in v1.items():
                    ofile.write("  "+k2+"="+str(v2)+"\n")
                ofile.write("$end\n")

        return {'command': self._get_command()}

    def process_run(self, command):
        """
        Run the xtb executable.
        """
        from os import environ, system
        # Set the number of omp threads only if the variable is not present
        # in the environment
        if 'OMP_NUM_THREADS' not in environ:
            environ['OMP_NUM_THREADS'] = self.run_options['omp']
        if 'MKL_NUM_THREADS' not in environ:
            environ['OMP_NUM_THREADS'] = self.run_options['omp']

        if self.run_options['verbose']:
            if self.run_dir != '.':
                safe_print('Run directory', self.run_dir)
            safe_print('Executing command: ', command)

        # Run the command
        system("cd " + self.run_dir + "; " + command)

        return {'logname': self._get_logname(True)}

    def post_processing(self, logname, command):
        """
        Post processing the calculation.

        Returns:
            (BigDFT.Interop.XTBInterop.XTBLogfile): a representation of the
            logfile.
        """
        try:
            return XTBLogfile(sys=self.run_options["sys"], logname=logname)
        except IOError:
            raise ValueError("The logfile (", logname, ") does not exist.")

    def _ensure_run_directory(self):
        from futile.Utils import ensure_dir
        run_dir = self.run_options.get('run_dir', '.')
        # Create the run_dir if not exist
        if ensure_dir(run_dir) and self.run_options['verbose']:
            safe_print("Create the sub-directory '%s'" % run_dir)

        self.run_dir = run_dir

    def _get_command(self):
        from copy import deepcopy

        command = self.command
        command += " " + self.run_options.get('name', 'input') + ".tmol"

        # Check if it is a dry run
        if self.run_options['dry_run']:
            return command + " --define"

        if self._check_skip():
            return '''echo "skip"'''

        options = deepcopy(self.run_options)

        # Adjust the command line with options
        name = self.run_options.get("name", "xtbrun")
        command += " --namespace " + name

        kv_ops = ["opt", "gfnff", "md", "omd", "hess", "gfn"]
        kv_ops += ["alpb", "chrg", "gbsa", "metadyn", "path", "uhf"]

        for k in kv_ops:
            if k in options:
                if isinstance(options[k], bool) and options[k] is True:
                    command += " --"+str(k)
                elif not isinstance(options[k], bool):
                    command += " --"+str(k) + " " + str(options[k])
                options.pop(k)

        command += " > " + self._get_logname(False)

        return command

    def _check_skip(self):
        from os.path import exists, join
        if not self.run_options["skip"]:
            return False
        fpath = join(self.run_dir, self._get_logname(False))
        if not exists(fpath):
            return False
        with open(fpath) as ifile:
            for line in ifile:
                if "* finished run on" in line:
                    return True
        return False

    def _get_logname(self, full):
        from os.path import join
        name = self.run_options.get("name", "xtbrun") + ".log"
        if full:
            run_dir = self.run_options.get('run_dir', '.')
            return join(run_dir, name)
        else:
            return name

    def _get_inputfilename(self):
        from os.path import join
        name = self.run_options.get('name', '')
        run_dir = self.run_options.get('run_dir', '.')
        input_file = name + '.tmol' if name else 'input.tmol'
        return join(run_dir, input_file)

    def _set_constraints(self, constr, ofile):
        from itertools import groupby
        from operator import itemgetter
        clist = []
        i = 1
        for fragid, frag in constr.items():
            for at in frag:
                if at:
                    clist.append(i)
                i += 1

        clist.sort()

        if len(clist) > 0:
            ofile.write("$fix\n")
            ofile.write("  atoms: ")
            ranges = []
            # Create ranges from consecutive indices.
            for k, g in groupby(enumerate(clist), lambda x: x[0]-x[1]):
                trange = list(map(itemgetter(1), g))
                if len(trange) == 1:
                    ranges.append(str(trange[0]))
                elif len(trange) > 1:
                    ranges.append(str(trange[0]) + "-" + str(trange[-1]))
            ofile.write(",".join(ranges))
            ofile.write("\n$end\n")


class XTBLogfile():
    """
    This class stores the results of an XTB calculation which might be later
    post-processed.

    Attributes:
        sys (BigDFT.Systems.System): a handle to the original system used
            in the calculation.
        log (str): a string representation of the logfile from XTB.
        charges (list): charge values computed by XTB.
        opt_energy (list): a list of energy values from the geometry
            optimization procedure.
        opt_fnorm (list): a list of gradient norm values from the geometry
            optimization procedure.
        opt_traj (list): a list of BigDFT.Systems.System from the geometry
            optimization procedure.
        md_energy (list): a list of energy values from the MD simulation.
        md_fnorm (list): a list of gradient norm values from the MD simulation.
        md_traj (list): a list of BigDFT.Systems.System from the MD simulation.
    """
    def __init__(self, sys, logname, opt=False, md=False):
        self.sys = sys
        fpath = logname.replace(".log", "")
        with open(fpath+".log") as ifile:
            self.log = "".join(ifile)
        try:
            with open(fpath+".charges") as ifile:
                self.charges = [float(x) for x in ifile]
        except FileNotFoundError:
            self.charges = None
        self.opt_energy, self.opt_fnorm, self.opt_traj = \
            self._process_trajectory(fpath + ".xtbopt.log")
        self.md_energy, self.md_fnorm, self.md_traj = \
            self._process_trajectory(fpath + ".xtb.trj")
        self.optimized_geometry = \
            self._get_optimized_geometry(fpath + ".xtbopt.tmol")

    def _process_trajectory(self, fname):
        from BigDFT.Atoms import Atom
        try:
            with open(fname) as ifile:
                ene = []
                fnorm = []
                systems = []
                while(True):
                    try:
                        natoms = int(next(ifile))
                        eneline = next(ifile).split()
                        ene.append(float(eneline[1]))
                        fnorm.append(float(eneline[3]))

                        atlist = []
                        for i in range(natoms):
                            split = next(ifile).split()
                            pos = [float(x) for x in split[1:]]
                            sym = split[0].upper()
                            atlist.append(Atom(sym=sym, r=pos,
                                               units="angstroem"))
                        systems.append(self._match_atlist(atlist))
                    except StopIteration:
                        break
        except FileNotFoundError:
            return (None, None, None)

        return (ene, fnorm, systems)

    def _get_optimized_geometry(self, fname):
        from BigDFT.Atoms import Atom
        try:
            with open(fname) as ifile:
                atlist = []
                line = next(ifile)
                while("$coord" not in line):
                    line = next(ifile)
                line = next(ifile)
                while("$" not in line):
                    split = line.split()
                    pos = [float(x) for x in split[:3]]
                    sym = split[-1].upper()
                    atlist.append(Atom(sym=sym, r=pos,
                                       units="bohr"))
                    line = next(ifile)
                system = self._match_atlist(atlist)
        except FileNotFoundError:
            return None

        return system

    def _match_atlist(self, atlist):
        from BigDFT.Systems import System
        from BigDFT.Fragments import Fragment
        from copy import deepcopy

        optsys = System()
        i = 0
        for fragid, frag in self.sys.items():
            optsys[fragid] = Fragment()
            for at in frag:
                updated = deepcopy(at)
                updated.set_position(atlist[i].get_position())
                optsys[fragid] += Fragment([updated])
                i += 1

        return optsys

    @property
    def energy(self):
        """
        The total energy of the system.
        """
        for line in self.log.split("\n"):
            if "| TOTAL ENERGY " in line:
                return float(line.split()[3])
        else:
            raise ValueError("Total energy not available")

    @property
    def gradient_norm(self):
        """
        The total gradient of the system.
        """
        for line in self.log.split("\n"):
            if "| GRADIENT NORM " in line:
                return float(line.split()[3])
        else:
            raise ValueError("Gradient not available")

    @property
    def homo_lumo_gap(self):
        """
        The homo lumo gap of the system.
        """
        for line in self.log.split("\n"):
            if "| HOMO-LUMO GAP " in line:
                return float(line.split()[3])
        else:
            raise ValueError("Gap not available")

    @property
    def gfnff_qest(self):
        """
        When using the GFNFF option, no charges are computed. Instead, there
        are estimated charges when the topology is generated. This will
        extract those values.
        """
        qest = []
        for i, line in enumerate(self.log.split("\n")):
            if "metchar sp-hybrid imet pi  qest" in line:
                break
        else:
            raise ValueError("qest not found in log file")

        for line in self.log.split("\n")[i+1:]:
            split = line.split()
            if len(split) == 0:
                break
            qest += [float(line.split()[8])]

        return qest


def xtb_equilibrate(sys, calc, prod_time=20, **kwargs):
    """
    This function will automate the process of taking a system from zero
    temperature up to room temperature (298.15 kelvin).

    This routine will work by performing 10ps simulations at temperatures
    from 0 to 298.15 kelvin, with each simulation spaced by 10 degrees
    (310 ps total). Finally, the production run can be performed.

    Args:
        sys (BigDFT.Systems.System): the system to calculate.
        calc (BigDFT.Interop.XTBInterop.XTBCalculator): a calculator to use.
        prod_time (int): the number of production picoseconds to run at
            room temperature.
        **kwargs: any argument you would normally want to pass to the
            run command.
    Returns:
        (list, list): a list of positions along the trajectory and a list of
        energy values.
    """
    from copy import deepcopy
    positions = []
    energies = []
    param = kwargs.pop("inp", {})
    if "md" not in param:
        param["md"] = {}

    # Loop over temperature values
    usys = deepcopy(sys)
    for temp in list(range(0, 300, 10)) + [298.15]:
        param["md"]["temp"] = temp
        param["md"]["time"] = 10
        if temp == 0:
            param["md"]["restart"] = "false"
            log = calc.run(sys=usys, omd=True, inp=param, **kwargs)
        else:
            param["md"]["restart"] = "true"
            log = calc.run(sys=usys, md=True, inp=param, **kwargs)

        # Extract results
        positions += log.md_traj
        energies += log.md_energy
        usys = positions[-1]

    # Production Run
    param["md"]["time"] = prod_time
    param["md"]["restart"] = "true"
    param["md"]["temp"] = 298.15
    log = calc.run(sys=usys, md=True, inp=param, **kwargs)

    positions += log.md_traj
    energies += log.md_energy
    usys = positions[-1]

    return positions, energies


def _example():
    """Example of using XTB interoperability"""
    from BigDFT.IO import XYZReader
    from BigDFT.Systems import System
    from BigDFT.Fragments import Fragment
    from os.path import join
    from copy import deepcopy

    # Create a system.
    reader = XYZReader("Ar")
    fsys = System()
    fsys["FRA:1"] = Fragment(xyzfile=reader)
    fsys["FRA:2"] = deepcopy(fsys["FRA:1"])
    fsys["FRA:2"].translate([-2, 0, 0])

    # Create a calculator.
    code = XTBCalculator()
    log = code.run(sys=fsys, name="test", run_dir="work")

    # Print some values.
    print(log.energy)
    print(log.gradient_norm)
    print(log.homo_lumo_gap)
    print(log.charges)

    # Geometry optimization.
    log = code.run(sys=fsys, name="opt", opt=True, run_dir="work")
    sys2 = log.optimized_geometry
    # this also works
    sys2 = log.opt_traj[-1]
    print(log.opt_energy)
    print(log.opt_fnorm)

    # Recompute.
    log = code.run(sys=sys2, name="opt2", run_dir="work")
    print(log.energy)
    print(log.gradient_norm)
    print(log.homo_lumo_gap)

    # Molecular dynamics.
    inp = {}
    inp["md"] = {"temp": 300, "time": 2}
    logmd = code.run(sys=fsys, name="md", omd=True, inp=inp, run_dir="work")
    print(logmd.md_energy)
    print(logmd.md_fnorm)
    sys3 = logmd.md_traj[-1]


if __name__ == "__main__":
    _example()
