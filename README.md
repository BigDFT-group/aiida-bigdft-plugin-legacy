[![Build Status](https://travis-ci.org/BigDFT-group/aiida-bigdft-plugin.svg?branch=master)](https://travis-ci.org/BigDFT-group/aiida-bigdft-plugin) 
[![Coverage Status](https://coveralls.io/repos/github/BigDFT-group/aiida-bigdft-plugin/badge.svg?branch=master)](https://coveralls.io/github/BigDFT-group/aiida-bigdft-plugin?branch=master) 
[![Docs status](https://readthedocs.org/projects/aiida-bigdft/badge)](http://aiida-bigdft.readthedocs.io/) 
[![PyPI version](https://badge.fury.io/py/aiida-bigdft.svg)](https://badge.fury.io/py/aiida-bigdft)

# aiida-bigdft

Aiida plugin for BigDFT code.
This is a simple plugin to integrate bigdft computation in an AiiDA workflow. Input file is generated using PyBigDFT tools, and output is returned as a dict, using the LogFile feature of PyBigDFT. Both futile and PyBigDFT are included in this plugin, BigDFT itself is not needed.

This plugin is the default output of the 
[AiiDA plugin cutter](https://github.com/aiidateam/aiida-plugin-cutter),
intended to help developers get started with their AiiDA plugins.

## Features

 * Create input files and specify command line options via a python dictionary and `BigDFTParameters`:
   ```python
   d = { 'ignore-case': True }
   BigDFTParameters = DataFactory('bigdft')
   bigdft_parameters = {}
   bigdft_parameters["dft"] = { "ixc": "LDA", "itermax": "5" }
   bigdft_parameters["output"] = { 'orbitals': 'binary' } 
   inputs['parameters'] = BigDFTParameters(dict=bigdft_parameters)
   ```

 * Run computation and retrieve output files (by default : logfile, time.yaml, forces, can be extended to retrieve any geenrated file):
   ```python
   inputs['extra_retrieved_files'] = List()
   inputs['extra_retrieved_files'].set_list([["./data*/*", ".", 2]])
   result = run(CalculationFactory('bigdft'), **inputs)
   #or asynchronously
   future = submit(CalculationFactory('bigdft'), **inputs)
   ```

 * load back YAML logfile and turn it into a python dict (through PyBigDFT utilities) to analyze results
   ```python
   #only if run asynchronously, load results from database first after completion
   result=load_node(future.pk).outputs
   
   logfile = result['bigdft_logfile'].logfile
   print (logfile['Energy (Hartree)'])
   ```
## Installation

```shell
pip install aiida-bigdft
verdi quicksetup  # better to set up a new profile, or run reentry scan
verdi plugin list aiida.calculations  # should now show your calclulation plugins
```

## Usage

Here goes a complete example of how to submit a test calculation using this plugin.

A quick demo of how to submit a calculation:
```shell
verdi daemon start         # make sure the daemon is running
cd examples
verdi run submit.py        # submit test calculation
verdi process list -a  # check status of calculation
```

The plugin also includes verdi commands to inspect its data types:
```shell
verdi data bigdft list
verdi data bigdft export <PK>
```

## Development

```shell
git clone https://github.com/BigDFT-group/aiida-bigdft-plugin .
cd aiida-bigdft-plugin
pip install -e .[pre-commit,testing]  # install extra dependencies
pre-commit install  # install pre-commit hooks
pytest -v  # discover and run all tests
```

See the [developer guide](http://aiida-bigdft.readthedocs.io/en/latest/developer_guide/index.html) for more information.

## License

MIT


## Contact

bigdft-developers@lists.launchpad.net

