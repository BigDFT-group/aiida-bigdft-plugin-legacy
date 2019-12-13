""" Tests for calculations

"""
from __future__ import print_function
from __future__ import absolute_import

import os
from aiida_bigdft import tests


def test_process(aiida_code):
    """Test running a calculation
    note this does not test that the expected outputs are created of output parsing"""
    from aiida.plugins import DataFactory, CalculationFactory
    from aiida.engine import run
    #TODO
