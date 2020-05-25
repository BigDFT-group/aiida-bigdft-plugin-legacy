"""
Basic wrapping workchain on a BigDFT computation.

"""
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import BaseRestartWorkChain, process_handler, while_, ExitCode
from aiida.plugins import CalculationFactory, DataFactory

from aiida.engine.processes.workchains.utils import process_handler, ProcessHandlerReport

from futile import YamlIO

Dict = DataFactory('dict') 

RemoteData = DataFactory('remote') 
BigDFTCalculation = CalculationFactory('bigdft')

class BigDFTBaseWorkChain(BaseRestartWorkChain):

    _process_class = BigDFTCalculation

    @classmethod
    def define(cls, spec):
        super(BigDFTBaseWorkChain, cls).define(spec)
        spec.input('show_warnings', valid_type=orm.Bool,
                   default=lambda: orm.Bool(True),
                   help='turn the warnings on/off.')
        spec.input('run_opts', valid_type=Dict,
                    required=False,
                    help='metadata')
        spec.expose_inputs(BigDFTCalculation, exclude=('metadata',))

        spec.outline(
            cls.setup,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )
        spec.expose_outputs(BigDFTCalculation)
        #this one needs to be optional to avoid being checked wrongly by the restartworkchain
        spec.exit_code(100, 'ERROR_INPUT',
                       message='BigDFT input error')
        spec.exit_code(200, 'ERROR_RUNTIME',
                       message='BigDFT runtime error')

    @process_handler(priority=600)
    def check_debug_output(self, calculation):
        repo = calculation.outputs.retrieved._repository._get_base_folder()
        try:
            repo.get_abs_path('debug', check_existence=True)
        except OSError:
            return
        debug_folder = repo.get_subfolder('debug')
        # debug folder exists, error probably happened.
        if "jobname" in self.ctx.inputs.metadata.options:
            jobname = self.ctx.inputs.metadata.options.jobname
        else:
            jobname = 'BigDFT job'
        posout_list = debug_folder.get_content_list(pattern="bigdft-err*")
        for filename in posout_list:
            log = YamlIO.load(debug_folder.get_abs_path(filename))
            err = log[0].get('BIGDFT_INPUT_VARIABLES_ERROR')
            info = log[0].get('Additional Info')
            if err is not None:
                self.report('{}<{}> input error : {} Id: {}.\n\
                            Additional Information : {}'.
                            format(jobname, calculation.pk,
                                   err['Message'], err['Id']))
                if info is not None:
                    self.report('Additional Info :', info)
                return ProcessHandlerReport(True, ExitCode(100))
            err = log[0].get('BIGDFT_RUNTIME_ERROR')
            if err is not None:
                self.report('{}<{}> runtime error : {} Id: {}'.
                            format(jobname, calculation.pk,
                                   err['Message'], err['Id']))
                if info is not None:
                    self.report('Additional Info :', info)
                return ProcessHandlerReport(True, ExitCode(200))

    @process_handler(priority=590)
    def check_warnings(self, calculation):
        if calculation.is_finished_ok and self.inputs.show_warnings:
            logfile=calculation.outputs.bigdft_logfile.logfile
            if(isinstance(logfile, list)):
                warnings = logfile[-1].get('WARNINGS')
            else:
                warnings = logfile.get('WARNINGS')
            if warnings is not None:
                self.report('Warnings were found :')
                for warn in warnings:
                    self.report(warn)

    @process_handler(priority=500)
    def finish(self, calculation):
        if calculation.is_finished_ok:
            if "jobname" in self.ctx.inputs.metadata.options:
                jobname = self.ctx.inputs.metadata.options.jobname
            else:
                jobname = 'BigDFT job'
            self.report('{}<{}> completed successfully'.
                        format(jobname, calculation.pk))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

    def setup(self):
        super().setup()
        self.ctx.inputs = AttributeDict(self.exposed_inputs(BigDFTCalculation))
        if self.inputs.get('run_opts') is not None:
            self.ctx.inputs.metadata = AttributeDict(self.inputs.run_opts.get_dict())
        else:
            self.ctx.inputs.metadata = {}


