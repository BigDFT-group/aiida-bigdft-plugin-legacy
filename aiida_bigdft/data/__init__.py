# -*- coding: utf-8 -*-
"""
Data types provided by plugin

Register data types via the "aiida.data" entry point in setup.json.
"""

# You can directly use or subclass aiida.orm.data.Data
# or any other data type listed under 'verdi data'
from __future__ import absolute_import
from aiida.orm.nodes import Data
from voluptuous import Optional

from BigDFT import Logfiles
from BigDFT.scripts import InputGenerator

# A subset of BigDFT's command line options
cmdline_options = {
    Optional('ignore-case'): bool,
    Optional('ignore-file-name-case'): bool,
    Optional('ignore-tab-expansion'): bool,
    Optional('ignore-space-change'): bool,
    Optional('ignore-all-space'): bool,
}


class BigDFTParameters(Data):
    """
    Command line options for BigDFT.

    This class represents a python dictionary used to
    pass command line options to the executable.
    """

    # "voluptuous" schema  to add automatic validation
    # pylint: disable=redefined-builtin
    def __init__(self, dict=None, **kwargs):
        """
        Constructor for the data class

        Usage: ``BigDFTParameters(dict{'ignore-case': True})``

        :param parameters_dict: dictionary with commandline parameters
        :param type parameters_dict: dict

        """
#        dict = self.validate(dict)
        super(BigDFTParameters, self).__init__(**kwargs)
        self.set_attribute('dict', dict)

    def set_dict(self, dict):
        self.set_attribute('dict', dict)

    @property
    def dict(self):
        """
        Return the Logfile
        """
        return self.get_attribute('dict')

    def cmdline_params(self, file1_name, file2_name):
        """Synthesize command line parameters.

        e.g. [ '--ignore-case', 'filename1', 'filename2']

        :param file_name1: Name of first file
        :param type file_name1: str
        :param file_name2: Name of second file
        :param type file_name2: str

        """

    def __str__(self):
        """String representation of node.

        Append values of dictionary to usual representation. E.g.::

            uuid: b416cbee-24e8-47a8-8c11-6d668770158b (pk: 590)
            {'ignore-case': True}

        """
        string = super(BigDFTParameters, self).__str__()
        string += "\n" + str(self.get_attribute('dict'))
        return string

    def set_inputfile(hgrid, dico, init_input=None, psp=[], units="reduced"):
        return InputGenerator.set_inputfile(hgrid, dico, init_input, psp, units)

    def set_spin(name, nat):
        return InputGenerator.set_spin(name, nat)

    def set_kpoints(nat):
        return InputGenerator.set_kpoints(nat)

    def set_strain(strain, ngrids, dico):
        return InputGenerator.set_strain(strain, ngrids, dico)

    def set_restart():
        return InputGenerator.set_restart()

    def transform_to_orthorombic(dico):
        return InputGenerator.transform_to_orthorombic(dico)

    def set_psp(name, psp):
        return InputGenerator.set_psp(name, psp)


class BigDFTLogfile(Data):
    """
    Wrapper around a Logfile object

    This class represents a BigDFT Logfile object.
    """
    def __init__(self, path=None, **kwargs):
        super(BigDFTLogfile, self).__init__(**kwargs)
        self.logfile = path

    @property
    def logfile(self):
        """
        Return the Logfile
        """
        return self.get_attribute('logfile')

    @logfile.setter
    def logfile(self, path):
        """
        Set the logfile

        :raise ValueError:
        """
        self.bigdftlogfile = Logfiles.Logfile(path)
        # replace forbidden chars in aiida dicts.
        if len(self.bigdftlogfile) > 0:
            logs = []
            for log in self.bigdftlogfile:
                try:
                    log.log['Timestamp of this run'] = \
                        log.log['Timestamp of this run'].strftime("%Y-%m-%d %H:%M:%S.%f")
                except KeyError:
                    pass
                logs.append(log.log)
            self.set_attribute('logfile', logs)
        else:
            self.bigdftlogfile.log['Timestamp of this run'] = \
                self.bigdftlogfile.log['Timestamp of this run'].strftime("%Y-%m-%d %H:%M:%S.%f")
            self.set_attribute('logfile', self.bigdftlogfile.log)

    def __str__(self):
        """String representation of node.

        Append values of dictionary to usual representation. E.g.::

            uuid: b416cbee-24e8-47a8-8c11-6d668770158b (pk: 590)
            {'ignore-case': True}

        """
        string = super(BigDFTLogfile, self).__str__()
        string += "\n" + str(self.get_attribute('logfile'))
        return string
