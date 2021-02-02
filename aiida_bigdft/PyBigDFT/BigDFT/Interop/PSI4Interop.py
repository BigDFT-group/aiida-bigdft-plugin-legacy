"""
This module contains some wrappers for using PSI4 to perform calculations.

https://psicode.org/
"""
from BigDFT.Calculators import Runner
from futile.Utils import write as safe_print


def bigdft_to_psi4(sysA, sysB=None, spinA=0, multiplicityA=1, spinB=0,
                   multiplicityB=1):
    """
    Create a PSI4 dimer molecule which can be used for SAPT.

    If multiple molecules are specified, we create a Dimer system.

    Args:
        sysA (BigDFT.Systems.System): the first molecule.
        sysB (BigDFT.Systems.System): the second molecule (optional).
        spinA (int): the spin of the first molecule.
        spinB (int): the spin of the second molecule.
        multiplicityA (int): the multiplicity for unpaired electrons.
        multiplicityB (int): the multiplicity for unpaired electrons.

    Returns:
        (psi4.core.Molecule): the psi4 geometry.
    """
    from BigDFT.IO import write_xyz
    from psi4 import geometry
    # py2 workaround
    from sys import version_info
    if version_info[0] < 3:
        from io import BytesIO as StringIO
    else:
        try:
            from io import StringIO
        except ImportError:
            from StringIO import StringIO

    # Create the right string representations
    strA = StringIO()
    write_xyz(sysA, strA)
    xyzA = str(spinA) + " " + str(multiplicityA) + "\n"
    xyzA += "\n".join(strA.getvalue().split("\n")[2:])

    if sysB is not None:
        strB = StringIO()
        write_xyz(sysB, strB)
        xyzB = str(spinB) + " " + str(multiplicityB) + "\n"
        xyzB += "\n".join(strB.getvalue().split("\n")[2:])

    if sysB is None:
        return geometry("\n"+xyzA+"\nunits angstrom\n")
    else:
        return geometry("\n"+xyzA+"--\n"+xyzB+"units angstrom\n")


class PSI4Calculator(Runner):
    """
    Perform a calculation on a given system using the PSI4 code.

    Note that if you intend to use an SAPT method, you need to pass a
    system which is composed of exactly two fragments.

    PSI4 has a number of different actions, ab initio methods, and basis
    sets. Be sure to specify each of these to the run command.
    """
    import os

    def __init__(self, omp=os.environ.get('OMP_NUM_THREADS', '1'),
                 skip=False, verbose=True):
        # Use the initialization from the Runner class (so all options inside
        # __global_options)
        self.mol = None
        Runner.__init__(self, omp=str(omp), dry_run=False, skip=skip,
                        verbose=verbose)

    def pre_processing(self):
        """
        Process local run dictionary to create the input directory and identify
        the command to be passed

        Returns:
            :py:class:`dict`: dictionary containing the command to be passed to
            :meth:`process_run`
        """
        from psi4 import set_options, energy, optimize, frequency
        from psi4.core import set_output_file
        from BigDFT.Systems import System
        from os.path import join

        # Check Arguments
        try:
            self.sys = self.run_options["sys"]
        except KeyError:
            raise ValueError("sys= must be provided as an argument")

        try:
            action = self.run_options["action"]
        except KeyError:
            raise ValueError("You must specify a valid action=(energy,"
                             " optimize, frequency")

        try:
            method = self.run_options["method"]
        except KeyError:
            raise ValueError("You must specify an ab initio method=")

        try:
            basis = self.run_options["basis"]
        except KeyError:
            raise ValueError("You must specify a basis=")

        # Run directory
        self._ensure_run_directory()

        # Check skip
        if self.run_options["skip"]:
            try:
                log = PSI4Logfile(self.sys, self._get_logname(True),
                                  action, method)
                return {"command": None}
            except ValueError:  # invalid logfile, so we can't skip.
                pass
            except IOError:  # no logfile, so we can't skip.
                pass

        # Convert system
        if "sapt" in method:
            keys = list(self.sys)
            if len(keys) != 2:
                raise ValueError("For SAPT method, your system must be"
                                 " composed of exactly two fragments.")
            sysA = System()
            sysA[keys[0]] = self.sys[keys[0]]
            sysB = System()
            sysB[keys[1]] = self.sys[keys[1]]
            mol = bigdft_to_psi4(sysA=sysA, sysB=sysB)
        else:
            mol = bigdft_to_psi4(sysA=self.sys)

        # Set extra options.
        if "psi4_options" in self.run_options:
            set_options(self.run_options["psi4_options"])
        set_output_file(self._get_logname(True), False)

        # Instead of a command, we return a closure.
        if action == "energy":
            def cmd(method, mol):
                energy(method, molecule=mol)
        elif action == "optimize":
            def cmd(method, mol):
                optimize(method, molecule=mol)
        elif action == "frequency":
            def cmd(method, mol):
                frequency(method, molecule=mol)

        return {"command": cmd(join(method, basis), mol)}

    def process_run(self, command):
        """
        Run the psi4 command.
        """
        from os import environ, system
        # Set the number of omp threads only if the variable is not present
        # in the environment
        if 'OMP_NUM_THREADS' not in environ:
            environ['OMP_NUM_THREADS'] = self.run_options['omp']

        if self.run_options['verbose']:
            if self.run_dir != '.':
                safe_print('Run directory', self.run_dir)

        # Run the command
        if command is not None:
            command()

        return {'logname': self._get_logname(True)}

    def post_processing(self, logname, command):
        """
        Post processing the calculation.

        Returns:
            (BigDFT.Interop.PSI4Interop.PSI4Logfile): a representation of the
            logfile.
        """
        try:
            return PSI4Logfile(self.sys, fname=logname,
                               action=self.run_options["action"],
                               method=self.run_options["method"])
        except IOError:
            raise ValueError("The logfile (", logname, ") does not exist.")

    def _ensure_run_directory(self):
        from futile.Utils import ensure_dir
        run_dir = self.run_options.get('run_dir', '.')
        # Create the run_dir if not exist
        if ensure_dir(run_dir) and self.run_options['verbose']:
            safe_print("Create the sub-directory '%s'" % run_dir)

        self.run_dir = run_dir

    def _get_logname(self, full):
        from os.path import join
        name = self.run_options.get("name", "sapt") + ".log"
        if full:
            run_dir = self.run_options.get('run_dir', '.')
            return join(run_dir, name)
        else:
            return name


class PSI4Logfile(dict):
    """
    A logfile wrapper for an SAPT calculation.

    This logfile inherits from a dictionary, as the properties we get out
    of a calculation depend on the type of calculation performed.

    Args:
        fname (str): the name of the file to process.
        action (str): the action that was used when running SAPT (energy,
          optimize, frequency).
        method (str): the scf method used.
    """
    def __init__(self, sys, fname, action, method):
        dict.__init__({})
        self["energy"] = {}

        self.sys = sys

        # The logfile, broken up into a list of lines.
        self.log = self._log_to_str(fname)

        if "sapt" in method:
            self._extract_sapt_energies(method)
        elif "energy" in action:
            self._extract_total_energy()
        elif "optimize" in action:
            self._extract_optimization_energies()
            self._extract_optimization_geometries()

    def _log_to_str(self, fname):
        with open(fname) as ifile:
            lstr = "".join([x for x in ifile])
        return lstr.split("\n")

    def _look_for(self, goal):
        for offset, line in enumerate(self.log):
            if goal in line:
                break
        if offset > len(self.log):
            raise ValueError("Couldn't find property in logfile: ", goal)
        return offset

    def _extract_optimization_energies(self):
        self["max force"] = []
        self["energy"]["optimize"] = []
        offset = self._look_for("==> Optimization Summary <==")
        for line in self.log[offset+6:]:
            if "---" in line:
                break
            split = line.split()
            self["max force"].append(float(split[3]))
            self["energy"]["optimize"].append(float(split[1]))
        self["energy"]["total"] = self["energy"]["optimize"][-1]

    def _extract_optimization_geometries(self):
        from BigDFT.Systems import System
        from BigDFT.Fragments import RotoTranslation
        from copy import deepcopy

        # Read out the posinps
        posinps = []
        geomcount = 0
        for i in range(len(self.log)):
            line = self.log[i]
            if "==> Geometry <==" in line:
                # It is printed twice, once for scf and one for gradient.
                geomcount += 1
                if geomcount % 2 == 0:
                    continue
                # Extract
                geom = {"positions": [], "units": "angstroem"}
                j = 9
                while(True):
                    split = self.log[i + j].split()
                    if len(split) == 0:
                        break
                    atom = {split[0].title(): [float(x) for x in split[1:4]]}
                    geom["positions"].append(atom)
                    j += 1
                posinps.append(geom)

        # Convert the posinps into first class system objects.
        self["positions"] = []
        for pos in posinps:
            self["positions"].append(deepcopy(self.sys))
            self["positions"][-1].update_positions_from_dict(pos)

        # PSI4 does rototranslations of the geometry, so we need to undo it.
        rt = {}
        for fragid in self.sys:
            rt[fragid] = RotoTranslation(self["positions"][0][fragid],
                                         self.sys[fragid])
        for osys in self["positions"]:
            for fragid in self.sys:
                osys[fragid] = rt[fragid].dot(osys[fragid])

    def _extract_total_energy(self):
        for line in self.log:
            if "Total Energy =" in line:
                self["energy"]["total"] = float(line.split()[-1])
        if "total" not in self["energy"]:
            raise ValueError("Could you find total energy in log file.")

    def _extract_sapt_energies(self, method):
        def _check_vals():
            if any([x not in vals for x in ["total", "exchange", "induction",
                                            "dispersion", "electrostatics"]]):
                raise ValueError("Could not find sapt energies in log file")

        # Get main computation
        vals = {}
        offset = self._look_for("SAPT Results")
        for line in self.log[offset:]:
            if "Electrostatics" in line:
                vals["electrostatics"] = float(line.split()[1]) / 1000.0
            elif "Exchange" in line:
                vals["exchange"] = float(line.split()[1]) / 1000.0
            elif "Induction" in line:
                vals["induction"] = float(line.split()[1]) / 1000.0
            elif "Dispersion" in line:
                vals["dispersion"] = float(line.split()[1]) / 1000.0
            elif "Total" in line:
                vals["total"] = float(line.split()[2]) / 1000.0
                break

        _check_vals()
        self["energy"][method] = vals

        # Extract special recipe version
        vals = {}
        for line in self.log:
            if "Total sSAPT0" in line:
                vals["total"] = float(line.split()[2]) / 1000.0
            if "Exchange sSAPT0" in line:
                vals["exchange"] = float(line.split()[2]) / 1000.0
            if "Induction sSAPT0" in line:
                vals["induction"] = float(line.split()[2]) / 1000.0
            if "Dispersion sSAPT0" in line:
                vals["dispersion"] = float(line.split()[2]) / 1000.0
            if "Electrostatics sSAPT0" in line:
                vals["electrostatics"] = float(line.split()[2]) / 1000.0

        _check_vals()
        self["energy"]["sSAPT0"] = vals


def _example():
    """Example of using PSI4 interoperability"""
    from BigDFT.IO import XYZReader
    from BigDFT.Systems import System
    from BigDFT.Fragments import Fragment
    from os.path import join
    from copy import deepcopy

    # Create a system.
    reader = XYZReader(join("Database", "XYZs", "He.xyz"))
    fsys = System()
    fsys["FRA:1"] = Fragment(xyzfile=reader)
    fsys["FRA:2"] = deepcopy(fsys["FRA:1"])
    fsys["FRA:2"].translate([-2, 0, 0])

    # Create a calculator.
    code = PSI4Calculator()
    log = code.run(sys=fsys, action="energy", basis="jun-cc-pvdz",
                   psi4_options={"scf_type": "direct"}, method="scf",
                   name="test")
    print(log)

    # The raw logfile is also available.
    print(log.log[4])

    # Geometry optimization.
    log = code.run(sys=fsys, action="optimize", basis="jun-cc-pvdz",
                   method="scf", name="test-opt")
    print(log)

    # SAPT
    log = code.run(sys=fsys, action="energy", basis="jun-cc-pvdz",
                   method="sapt0", name="test-sapt")
    print(log)


if __name__ == "__main__":
    _example()
