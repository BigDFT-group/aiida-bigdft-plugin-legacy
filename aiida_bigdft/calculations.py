"""
Calculations provided by aiida_bigdft.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from __future__ import absolute_import

import six

from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida.orm import SinglefileData
from aiida.plugins import DataFactory

from BigDFT import Calculators as BigDFT_calc
from BigDFT import InputActions as BigDFT_input

#BigDFTParameters = DataFactory('bigdft')


class BigDFTCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the BigDFT python interface.
    """
    _POSINP_FILE_NAME = 'posinp.xyz'
    
    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super(BigDFTCalculation, cls).define(spec)
        spec.input('metadata.options.resources', valid_type=dict, default={'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default='bigdft')
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default='log.yaml')
        spec.input('parameters', valid_type=orm.Dict, help='Command line parameters for BigDFT')
        spec.input('structure', valid_type=orm.StructureData, required=False, help='StructureData struct')
        spec.input('structurefile', valid_type=six.string_types, required=False, help='xyz file', default=_POSINP_FILE_NAME)
        spec.input('pseudos', valid_type=List, help='')
        
        spec.exit_code(100, 'ERROR_MISSING_OUTPUT_FILES', message='Calculation did not produce all expected output files.')


    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files needed by
            the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        #check that either structure or structurefile are set
        #if both are set, fail.
        if self.inputs.structurefile is None and self.inputs.structure is None:
            raise exceptions.InputValidationError("either structure or structurefile must be set")
        if self.inputs.structurefile is not None and self.inputs.structure is not None:
            raise exceptions.InputValidationError("Only one of structure or structurefile must be set")
        dico = self.inputs.parameters
        posinp_filename=self.inputs.structurefile
        if self.inputs.structure is not None
            posinp_string = self.inputs.structure._prepare_xyz()[0]
            posinp_file = open(_POSINP_FILE_NAME,"w") 
            posinp_filename = _POSINP_FILE_NAME
        BigDFT_input.set_atomic_positions(dico, posinp_filename)
        

        local_copy_list=[posinp_filename]
        #setup pseudopotentials if needed
        if self.inputs.pseudos is not None
            for filename in raw_pseudos:
                local_copy_list += os.path.abspath(filename)
        
        #generate yaml input file from dict and whatever
        bigdft_calc=BigDFT_calc.SystemCalculator()
        bigdft_calc.update_global_options(dico)
        bigdft_calc.update_global_options(units=angstroem)
        bigdft_calc.pre_processing()
        
        
        codeinfo = datastructures.CodeInfo()
#        codeinfo.cmdline_params = self.inputs.parameters.cmdline_params()
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.withmpi = self.inputs.metadata.options.withmpi
        local_copy_list = []
        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.retrieve_list = [self.metadata.options.output_filename]

        return calcinfo
