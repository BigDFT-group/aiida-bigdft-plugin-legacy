#!/usr/bin/env python

from __future__ import absolute_import
from setuptools import setup, find_packages
import os
import json
import sys
if __name__ == '__main__':
    # Provide static information in setup.json
    # such that it can be discovered automatically
    # create symlinks to subdirectories for easier use
    try:
      os.symlink('aiida_bigdft/PyBigDFT/BigDFT/', 'BigDFT')
    except:
      pass
    try:
      os.symlink('aiida_bigdft/futile/src/python/futile/','futile')
    except:
      pass
    with open('setup.json', 'r') as info:
        kwargs = json.load(info)
    setup(
        packages=find_packages(exclude=["aiida", "*.pyc", "aiida_bigdft/futile", "aiida_bigdft/PyBigDFT"]),
        # this doesn't work when placed in setup.json (something to do with str type)
        package_data={
            "": ["*"],
            # TODO be more specific with package data (but the line below isn't working)
            # "tests.input_files": ["*"],
        },
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        **kwargs)
