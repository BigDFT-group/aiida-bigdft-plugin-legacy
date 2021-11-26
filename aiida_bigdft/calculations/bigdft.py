# -*- coding: utf-8 -*-

"""
Calculations provided by aiida_bigdft.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from __future__ import absolute_import

import six
import os

from aiida import orm
from aiida.common import datastructures, exceptions
from aiida.engine import CalcJob
from aiida.orm.nodes.data import List, SinglefileData, Float, Dict
from aiida.plugins import DataFactory

from BigDFT import Calculators as BigDFT_calc
from BigDFT import InputActions as BigDFT_input
from BigDFT import Inputfiles as BigDFT_files
BigDFTParameters = DataFactory('bigdft')
BigDFTLogfile = DataFactory('bigdft_logfile')

# just override SystemCalculator constructor to avoid checking BIGDFT_ROOT variable


class PluginSystemCalculator(BigDFT_calc.SystemCalculator):
    def __init__(self,
                 omp=os.environ.get('OMP_NUM_THREADS', '1'),
                 mpi_run=os.environ.get('BIGDFT_MPIRUN', ''),
                 dry_run=False, skip=False, verbose=True):
        # Use the initialization from the Runner class (so all options inside
        # __global_options)
        BigDFT_calc.Runner.__init__(self, omp=str(omp), mpi_run=mpi_run,
                                    dry_run=dry_run, skip=skip, verbose=verbose)
        self.command = ""

    # run_dir location in aiida case does not really matter, as it's temporary only.
    def _ensure_run_directory(self):
        from futile.Utils import ensure_dir
        run_dir = self.run_options.get('run_dir', '.')
        # Create the run_dir if not exist
        if ensure_dir(run_dir) and self.run_options['verbose']:
            print("Create the sub-directory '%s'" % run_dir)
        self.run_dir = run_dir


class BigDFTCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the BigDFT python interface.
    """
    _POSINP_FILE_NAME = 'posinp.xyz'
    _INPUT_FILE_NAME = 'input.yaml'
    _OUTPUT_FILE_NAME = 'log.yaml'
    _TIMING_FILE_NAME = 'time.yaml'

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super(BigDFTCalculation, cls).define(spec)
        spec.input('metadata.options.resources', valid_type=dict, default={
                   'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        spec.input('metadata.options.parser_name',
                   valid_type=six.string_types, default='bigdft')
        spec.input('metadata.options.output_filename',
                   valid_type=six.string_types, default=cls._OUTPUT_FILE_NAME)
        spec.input('metadata.options.jobname',
                   valid_type=six.string_types, required=False)
        spec.input('parameters', valid_type=BigDFTParameters,
                   help='Command line parameters for BigDFT')
        spec.input('structure', valid_type=orm.StructureData,
                   help='StructureData struct')
        spec.input('structurefile', valid_type=orm.Str, help='xyz file',
                   default=lambda: orm.Str(cls._POSINP_FILE_NAME))
        spec.input('pseudos', valid_type=List, help='',
                   default=lambda: List(), required=False)
        spec.input('kpoints', valid_type=Dict, help='kpoint mesh or kpoint path',
                   default=lambda: Dict(dict={}), required=False)
        spec.input('extra_retrieved_files', valid_type=List,
                   help='', default=lambda: List())
        spec.output('bigdft_logfile', valid_type=BigDFTLogfile,
                    help='BigDFT log file as a dict')
        spec.exit_code(100, 'ERROR_MISSING_OUTPUT_FILES',
                       message='Calculation did not produce all expected output files.')
        spec.exit_code(101, 'ERROR_PARSING_FAILED',
                       message='Calculation did not produce all expected output files.')
        spec.exit_code(400, 'ERROR_OUT_OF_WALLTIME',
                       message='Calculation did not finish because of a walltime issue.')
        spec.exit_code(401, 'ERROR_OUT_OF_MEMORY',
                       message='Calculation did not finish because of memory limit')

    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files needed by
            the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        # check that either structure or structurefile are set
        # if both are set, fail.
#        if self.inputs.structurefile is None and self.inputs.structure is None:
#            raise exceptions.InputValidationError("either structure or structurefile must be set")
#        if self.inputs.structurefile is not None and self.inputs.structure is not None:
#            raise exceptions.InputValidationError("Only one of structure or structurefile must be set")

        # somehow aiida sends back unicode strings here
        dico = BigDFT_files.Inputfile()
        dico.update(self.inputs.parameters.dict)

        bigdft_calc = PluginSystemCalculator()
        local_copy_list = []
        # if the structure is not already inside input dict
        if 'posinp' not in dico:
            posinp_filename = self.inputs.structurefile.value
            if self.inputs.structure is not None:
                if self.inputs.structure.cell_angles != [90.0, 90.0, 90.0]:
                    raise ValueError('non orthorhombic cells are not supported')
                print("writing input posinp file")
                posinp_string = self.inputs.structure._prepare_xyz()[0]
                # set bcs at the correct format (periodic only?)
                if self.inputs.structure.pbc == (True, True, True):
                    filestring = posinp_string.split(b'\n')
                    line = "periodic " + str(self.inputs.structure.cell_lengths[0]) + " " + str(self.inputs.structure.cell_lengths[1]) + " " + str(self.inputs.structure.cell_lengths[2])
                    filestring[1] = line.encode()
                    posinp_string = b'\n'.join(filestring)

                if "jobname" not in self.inputs.metadata.options:
                    posinp_filename = self._POSINP_FILE_NAME
                else:
                    posinp_filename = self.inputs.metadata.options.jobname + ".xyz"
    #            bigdft_calc.update_global_options(units="angstroem")
                posinp_file = open(posinp_filename, "wb")

                posinp_file.write(posinp_string)
                posinp_file.close()
            posinp_filedata = SinglefileData(
                file=os.path.abspath(posinp_filename)).store()

    #        BigDFT_input.set_atomic_positions(dico, posinp_filename)
    #        bigdft_calc._run_options(posinp={})

            local_copy_list = [
                (posinp_filedata.uuid, posinp_filedata.filename, posinp_filedata.filename)]

        # setup pseudopotentials if needed
        if self.inputs.pseudos is not None:
            for filename in self.inputs.pseudos:
                pseudo_filedata = SinglefileData(
                    file=os.path.abspath(filename)).store()
                local_copy_list.append(
                    (pseudo_filedata.uuid, pseudo_filedata.filename, pseudo_filedata.filename))
        # generate yaml input file from dict and whatever

        if "jobname" in self.inputs.metadata.options:
            bigdft_calc.update_global_options(
                name=self.inputs.metadata.options.jobname)
        bigdft_calc._run_options(input=dico, run_dir=folder.abspath)
        bigdft_calc.pre_processing()
        if "jobname" not in self.inputs.metadata.options:
            input_filename = self._INPUT_FILE_NAME
        else:
            input_filename = self.inputs.metadata.options.jobname + ".yaml"
        input_filedata = SinglefileData(
            file=folder.get_abs_path(input_filename)).store()
        local_copy_list.append(
            (input_filedata.uuid, input_filedata.filename, input_filename))

        codeinfo = datastructures.CodeInfo()
        codeinfo.code_uuid = self.inputs.code.uuid
        outfile = self.inputs.metadata.options.output_filename
        timefile = self._TIMING_FILE_NAME
        if "jobname" in self.inputs.metadata.options:
            outfile = "log-" + self.inputs.metadata.options.jobname + ".yaml"
            timefile = "time-" + self.inputs.metadata.options.jobname + ".yaml"

#        codeinfo.stdout_name = outfile
        codeinfo.withmpi = self.inputs.metadata.options.withmpi
        if "jobname" in self.inputs.metadata.options:
            codeinfo.cmdline_params = ["--name=" +
                                       self.inputs.metadata.options.jobname]
        #local_copy_list = []
        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.retrieve_list = [outfile,
                                  timefile,
                                  "forces_posinp.yaml",
                                  "forces_posinp.xyz",
                                  "final_posinp.yaml",
                                  "final_posinp.xyz",
                                  ["./debug/bigdft-err*", ".", 2]]
        calcinfo.retrieve_list.extend(self.inputs.extra_retrieved_files)
        return calcinfo
