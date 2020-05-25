"""
Relax workchain.

"""
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, append_, ToContext
from aiida.orm.nodes.data import List
from aiida.plugins import WorkflowFactory, DataFactory
from BigDFT import InputActions
import numpy
BigDFTBaseWorkChain = WorkflowFactory('bigdft')
BigDFTParameters = DataFactory('bigdft')
StructureData = DataFactory('structure')
ArrayData = DataFactory('array')


class BigDFTRelaxWorkChain(WorkChain):
    """Structure relaxation workchain."""
    _verbose = False

    @classmethod
    def define(cls, spec):
        super(BigDFTRelaxWorkChain, cls).define(spec)
        spec.expose_inputs(BigDFTBaseWorkChain,
                           exclude=['parameters',
                                    'extra_retrieved_files'])
        spec.input('parameters', valid_type=BigDFTParameters, required=False,
                   default=lambda: orm.Dict(), help='param dictionary')
        spec.input('extra_retrieved_files', valid_type=List, required=False,
                   help='', default=lambda: List())
        spec.input('relax.perform', valid_type=orm.Bool, required=False,
                   default=lambda: orm.Bool(True), help='perform relaxation')
        spec.input('relax.algo', valid_type=orm.Str,
                   default=lambda: orm.Str('FIRE'),
                   help='algorithm to use during relaxation')
        spec.input('relax.threshold_forces', valid_type=orm.Float,
                   required=False, help='energy cutoff value, in ev/Ang')
        spec.input('relax.steps', valid_type=orm.Int, required=False,
                   default=lambda: orm.Int(50),
                   help='number of relaxation steps to perform.')
        spec.outline(
            cls.relax,
            cls.results,
        )
        spec.expose_outputs(BigDFTBaseWorkChain)
        spec.output('relaxed_structure', valid_type=StructureData,
                    required=False)
        spec.output('forces', valid_type=ArrayData, required=False)
        spec.output('total_energy', valid_type=orm.Float, required=False)
        spec.exit_code(101, 'ERROR_FAILED_RELAX',
                       'Subprocess failed for relaxation')

    def relax(self):

        self.ctx.inputs = AttributeDict(
            self.exposed_inputs(BigDFTBaseWorkChain))
        inputdict = self.inputs.parameters.dict
        InputActions.optimize_geometry(inputdict,
                                       self.inputs.relax.algo.value,
                                       self.inputs.relax.steps.value)
        if self.inputs.relax.threshold_forces is not None:
            InputActions.dict_set(inputdict, 'geopt', 'forcemax',
                                  self.inputs.relax.threshold_forces.value/0.52917721067121)

        self.ctx.inputs.parameters = BigDFTParameters(dict=inputdict)

        # gather outputs
        extra_retrieved_files = List()
        if self.inputs.extra_retrieved_files is not None:
            extra_retrieved_files.set_list(
                self.inputs.extra_retrieved_files.get_list())
        extra_retrieved_files.extend([["./data*/*.xyz", ".", 2],
                                     ["./data*/geopt.mon", ".", 2],
                                     ["./final_*xyz", ".", 2]])
        self.ctx.inputs.extra_retrieved_files = extra_retrieved_files

        node = self.submit(BigDFTBaseWorkChain, **self.ctx.inputs)
        return ToContext(work=append_(node))

    def results(self):
        workchain = self.ctx.work[-1]

        if not workchain.is_finished_ok:
            self.report(
                'Relaxation failed with exit status {}'.
                format(workchain.exit_status))
            return self.exit_codes.ERROR_FAILED_RELAX

        # final_*xyz file not found, it did not finish.
        # We can output last *xyz file in data folder, restart, or fail
        outstruct = "final_posinp.xyz"
        repo = workchain.outputs.retrieved._repository._get_base_folder()
        if "jobname" in self.inputs.run_opts.get_dict()['options']:
            outstruct = "final_" +\
                self.inputs.run_opts.get_dict()['options']['jobname'] + ".xyz"
        try:
            sf = repo.get_abs_path(outstruct, check_existence=True)
        except OSError:
            # do we have posout files ?
            subname = "data"
            if "jobname" in self.inputs.run_opts.get_dict()['options']:
                subname = subname + "-" +\
                    self.inputs.run_opts.get_dict()['options']['jobname']

            data_folder = repo.get_subfolder(subname)
            posout_list = data_folder.get_content_list(pattern="posout*")
            if posout_list.empty():
                # not even, we failed. Should have been caught before.
                self.report(
                    'Relaxation failed - no output found')
                return self.exit_codes.ERROR_FAILED_RELAX
            sf = repo.get_abs_path(posout_list.sort()[-1][0])

        s = StructureData()
        # BigDFT xyz files have more data on the first line which confuse aiida
        # Get rid of them (and get energy and forces in the process)
        try:
            with open(sf) as f:
                content = f.readlines()
        except FileNotFoundError:
            self.report(
                'Relaxation failed - no output position file found')
            return self.exit_codes.ERROR_FAILED_RELAX

        firstline = content[0].split()
        content[0] = firstline[0] + '\n'
        self.out('total_energy', orm.Float(firstline[2]).store())

        forces_index = content.index(" forces\n")
        forces = content[forces_index + 1:]
        for i, f in enumerate(forces):
            forces[i] = f.split()[1:4]

        array_forces = ArrayData()
        array_forces.set_array('forces', numpy.array(forces))
        self.out('forces', array_forces.store())

        positions = content[0:forces_index]
        s._parse_xyz("".join(positions))
        self.out('relaxed_structure', s.store())
        self.out_many(self.exposed_outputs(workchain, BigDFTBaseWorkChain))
