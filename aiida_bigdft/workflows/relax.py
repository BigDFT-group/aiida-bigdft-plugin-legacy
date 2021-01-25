"""
Relax workchain.

"""
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, append_, ToContext
from aiida.orm.nodes.data import List
from aiida.plugins import WorkflowFactory, DataFactory
from BigDFT import InputActions
from futile import YamlIO

import numpy
BigDFTBaseWorkChain = WorkflowFactory('bigdft')
BigDFTParameters = DataFactory('bigdft')
StructureData = DataFactory('structure')
ArrayData = DataFactory('array')

HARTREE_TO_EV = 27.211386024367243
HARTREE_BOHR_TO_EV_AG = 51.42208619083232


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
        spec.input('relax.threshold_forces', valid_type=orm.Float, required=False,
                   default=lambda: orm.Float(0.0), help='energy cutoff value, in ev/Ang')
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
        if self.inputs.relax.perform:
            InputActions.optimize_geometry(inputdict,
                                           self.inputs.relax.algo.value,
                                           self.inputs.relax.steps.value)
        if self.inputs.relax.threshold_forces is not None:
            InputActions.dict_set(inputdict, 'geopt', 'forcemax',
                                  self.inputs.relax.threshold_forces.value / HARTREE_BOHR_TO_EV_AG)

        self.ctx.inputs.parameters = BigDFTParameters(dict=inputdict)

        # gather outputs
        extra_retrieved_files = List()
        if self.inputs.extra_retrieved_files is not None:
            extra_retrieved_files.set_list(
                self.inputs.extra_retrieved_files.get_list())
        extension = "xyz"
        posinp = self.inputs.parameters.dict.get('posinp')
        if posinp is not None:
            try:
                extension = posinp['properties']['format']
            except KeyError:
                extension = "yaml"
        extra_retrieved_files.extend([["./data*/*." + extension, ".", 2],
                                      ["./data*/geopt.mon", ".", 2],
                                      ["./final_*", ".", 2]])
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
        extension = "xyz"
        posinp = self.inputs.parameters.dict.get('posinp')
        if posinp is not None:
            try:
                extension = posinp['properties']['format']
            except KeyError:
                extension = "yaml"
        if self.inputs.relax.perform:
            outstruct = "final_posinp." + extension
            repo = workchain.outputs.retrieved._repository._get_base_folder()
            if "jobname" in self.inputs.run_opts.get_dict()['options']:
                outstruct = "final_" +\
                    self.inputs.run_opts.get_dict(
                    )['options']['jobname'] + ".xyz"
            try:
                sf = repo.get_abs_path(outstruct, check_existence=True)
            except OSError:
                # final_*xyz file not found, it did not finish.
                # We can output last *xyz file in data folder, restart, or fail
                # do we have posout files ?
                subname = "data"
                if "jobname" in self.inputs.run_opts.get_dict()['options']:
                    subname = subname + "-" +\
                        self.inputs.run_opts.get_dict()['options']['jobname']

                data_folder = repo.get_subfolder(subname)
                posout_list = data_folder.get_content_list(pattern="posout*")
                if not posout_list:
                    # not even, we failed. Should have been caught before.
                    self.report(
                        'Relaxation failed - no output found')
                    return self.exit_codes.ERROR_FAILED_RELAX
                sf = repo.get_abs_path(posout_list.sort()[-1])
        else:
            # no relaxation performed, file is named forces_posinp.xyz .. or yaml

            outstruct = "forces_posinp." + extension
            repo = workchain.outputs.retrieved._repository._get_base_folder()
            if "jobname" in self.inputs.run_opts.get_dict()['options']:
                outstruct = "forces_" +\
                    self.inputs.run_opts.get_dict(
                    )['options']['jobname'] + ".xyz"
            try:
                sf = repo.get_abs_path(outstruct, check_existence=True)
            except OSError:
                self.report(
                    'Relaxation failed - no output found')
                return self.exit_codes.ERROR_FAILED_RELAX

        s = StructureData()
        # BigDFT xyz files have more data on the first line which confuse aiida
        # Get rid of them (and get energy and forces in the process)

        if extension == "xyz":
            try:
                with open(sf) as f:
                    content = f.readlines()
            except FileNotFoundError:
                self.report(
                    'Relaxation failed - no output position file found')
                return self.exit_codes.ERROR_FAILED_RELAX
            firstline = content[0].split()
            content[0] = firstline[0] + '\n'
            self.out('total_energy', orm.Float(float(firstline[2]) * HARTREE_TO_EV).store())

            forces_index = content.index(" forces\n")
            forces = content[forces_index + 1:]
            for i, f in enumerate(forces):
                forces[i] = f.split()[1:4]

            array_forces = ArrayData()
            array_forces.set_array('forces', numpy.array(forces))
            self.out('forces', array_forces.store())

            positions = content[0:forces_index]
            s._parse_xyz("".join(positions))
            try:
                s._adjust_default_cell(vacuum_addition=0.0,
                                       pbc=self.inputs.structure.pbc)
            except ValueError:
                # we are probably on a plane, not a volume
                pass
        elif extension == "yaml":
            logfile = workchain.outputs.bigdft_logfile.logfile
            if(isinstance(logfile, list)):
                energy = logfile[-1].get('Energy (Hartree)')
            else:
                energy = logfile.get('Energy (Hartree)')
            self.out('total_energy', orm.Float(energy * HARTREE_TO_EV).store())
            if(isinstance(logfile, list)):
                forces = logfile[-1].get('Atomic Forces (Ha/Bohr)')
            else:
                forces = logfile.get('Atomic Forces (Ha/Bohr)')
            for i, f in enumerate(forces):
                forces[i] = list(f.values())[0]
            array_forces = ArrayData()
            array_forces.set_array('forces', numpy.array(forces))
            self.out('forces', array_forces.store())
            logl = YamlIO.load(sf)
            if logl is None:
                self.report(
                    'Relaxation failed - no output position file found')
                return self.exit_codes.ERROR_FAILED_RELAX
            log = logl[0]
            s.cell = log['abc']
            for i in log['positions']:
                s.append_atom(symbols=list(i.keys())[
                              0], position=list(i.values())[0])
        else:
            self.report(
                'Relaxation failed - no parsing available for output type ' + extension)
            return self.exit_codes.ERROR_FAILED_RELAX
        self.out('relaxed_structure', s.store())
        self.out_many(self.exposed_outputs(workchain, BigDFTBaseWorkChain))
