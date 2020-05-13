# -*- coding: utf-8 -*-
"""
Parsers provided by aiida_bigdft.

Register parsers via the "aiida.parsers" entry point in setup.json.
"""
from __future__ import absolute_import

from aiida.engine import ExitCode
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory
from aiida.plugins import DataFactory
from aiida.common.exceptions import ValidationError

BigDFTCalculation = CalculationFactory('bigdft')
BigDFTLogfile = DataFactory('bigdft_logfile')


class BigDFTParser(Parser):
    """
    Parser class for parsing output of calculation.
    """

    def __init__(self, node):
        """
        Initialize Parser instance

        Checks that the ProcessNode being passed was produced
        by a BigDFTCalculation.

        :param node: ProcessNode of calculation
        :param type node: :class:`aiida.orm.ProcessNode`
        """
        from aiida.common import exceptions
        super(BigDFTParser, self).__init__(node)
        if not issubclass(node.process_class, BigDFTCalculation):
            raise exceptions.ParsingError("Can only parse BigDFTCalculation")

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.

        :returns: an exit code, if parsing fails
        (or nothing if parsing succeeds)
        """
        output_filename = self.node.get_option('output_filename')
        jobname = self.node.get_option('jobname')
        if jobname is not None:
            output_filename = "log-" + jobname + ".yaml"
        # Check that folder content is as expected
        files_retrieved = self.retrieved.list_object_names()
        files_expected = [output_filename]
        # Note: set(A) <= set(B) checks whether A is a subset of B
        if not set(files_expected) <= set(files_retrieved):
            self.logger.error("Found files '{}', expected to find '{}'".format(
                files_retrieved, files_expected))
            return self.exit_codes.ERROR_MISSING_OUTPUT_FILES

        # add output file
        self.logger.info("Parsing '{}'".format(output_filename))
#        print(self.retrieved._repository._get_base_folder().get_abs_path(output_filename))
        output = BigDFTLogfile(self.retrieved._repository._get_base_folder().
                               get_abs_path(output_filename))
        try:
            output.store()
        except ValidationError:
            self.logger.info("Impossible to store LogFile - ignoring '{}'".
                             format(output_filename))

#        with self.retrieved.open(output_filename, 'rb') as handle:
#            output_node = SinglefileData(file=handle)
#        output_dict_aiida=orm.Dict(dict=output_dict)
#        output_dict_aiida.store()
#        output_log_aiida=BigDFTLogfile(output)
        self.out('bigdft_logfile', output)

        return ExitCode(0)
