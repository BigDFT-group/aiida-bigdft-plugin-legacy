# -*- coding: utf-8 -*-

"""
Calculations provided by aiida_bigdft.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from __future__ import absolute_import

from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida.orm.nodes.data import List
from aiida.plugins import DataFactory

RemoteData = DataFactory('remote')


class ScriptCalculation(CalcJob):
    """
    AiiDA calculation to add post treatments to a computation workcahin.
    post treatment scripts are to be registered as codes in aiida.
    They are python scripts accepting one argument : a remotefolder where data is stored
    Output files are not specified and can be added to the extra_retrieved_files list
    """

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super(ScriptCalculation, cls).define(spec)
        spec.input('metadata.options.resources', valid_type=dict, default={'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        spec.input('bigdft_data_folder', valid_type=RemoteData, help='Folder to the BigDFT data folder')
        spec.input('retrieved_files', valid_type=List, help='', default=lambda: List())
        spec.exit_code(101, 'ERROR_SCRIPT',
                       'Script execution failed')


    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files needed by
            the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        codeinfo = datastructures.CodeInfo()
        codeinfo.code_uuid = self.inputs.code.uuid
        if not self.inputs.code.can_run_on(
                self.inputs.bigdft_data_folder.computer):
            self.report("This post-processing script {}\
                         can't run on  {} where data resides",
                        format(self.inputs.code, self.inputs.bigdft_data_folder.get_computer_name()))
            return self.exit_codes.ERROR_SCRIPT

        codeinfo.withmpi = False
        codeinfo.cmdline_params = [self.inputs.bigdft_data_folder.get_remote_path()]
        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()

        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = []
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.extend(self.inputs.retrieved_files)

        return calcinfo
