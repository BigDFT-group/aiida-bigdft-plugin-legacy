# -*- coding: utf-8 -*-
"""Submit a test calculation on localhost.

Usage: verdi run submit.py
"""
from __future__ import absolute_import
from __future__ import print_function
import os
from aiida_bigdft import tests, helpers
from aiida.plugins import DataFactory, CalculationFactory
from aiida.engine import run

# get code
computer = helpers.get_computer()
code = helpers.get_code(entry_point='bigdft', computer=computer)

# Prepare input parameters
BigDFTParameters = DataFactory('bigdft')
StructureData = DataFactory('structure')


# TiO2 cubic structure
alat = 4. # angstrom
cell = [[alat, 0., 0.,],
[0., alat, 0.,],
[0., 0., alat,],
   ]
s = StructureData(cell=cell)
s.append_atom(position=(alat/2.,alat/2.,alat/2.),symbols='Ti')
s.append_atom(position=(alat/2.,alat/2.,0.),symbols='O')
s.append_atom(position=(alat/2.,0.,alat/2.),symbols='O')

# set up calculation
inputs = {
    'code': code,
    'structure': s,
    'metadata': {
        'description': "Test job submission with the aiida_bigdft plugin",
        'options' : {
            'jobname': 'TiO2',
            'max_wallclock_seconds': 30 * 60
        }
    },
}

bigdft_parameters = {}
bigdft_parameters["dft"] = { "ixc": "LDA", "itermax": "5" }
bigdft_parameters["output"] = { 'orbitals': 'binary' } 
inputs['parameters'] = BigDFTParameters(dict=bigdft_parameters)

inputs['extra_retrieved_files'] = List()
inputs['extra_retrieved_files'].set_list([["./data*/*", ".", 2]])
# Note: in order to submit your calculation to the aiida daemon, do:
# from aiida.engine import submit
# future = submit(CalculationFactory('bigdft'), **inputs)
result = run(CalculationFactory('bigdft'), **inputs)

#get a dict from the Yaml outputfile, which was stored as a result.
BigDFT_logfile = result['bigdft_logfile'].logfile

#for extra retrieved_files
data_folder = result['retrieved']

print (BigDFT_logfile['Energy (Hartree)'])
print (BigDFT_logfile['Walltime since initialization'])

