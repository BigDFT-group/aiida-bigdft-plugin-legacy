"""
This module will generate a set of python unit tests that will automatically
run each of the `_example` subroutines in the PyBigDFT package.

python -m unittest discover -s $PYTHONPATH/BigDFT/ -p "check_examples.py"
"""
from pkgutil import walk_packages
from futile.Utils import write as safe_print
from future.utils import with_metaclass
from unittest import TestCase
import BigDFT
from contextlib import redirect_stdout


class MetaTest(type):
    """
    A meta test generator that creates test cases from the test spec.
    """
    def __new__(mcs, name, bases, dict):
        # Generate a closure based on a given module.
        def gen_test(name, module):
            def test(self):
                with open(name + ".txt", "w") as ofile:
                    with redirect_stdout(ofile):
                        safe_print("Running Test:", name)
                        module._example()
            return test

        # This dictionary will be a map from test names to test methods.
        dict = {}

        # Walk over all packages
        for _, modname, _ in walk_packages(path=BigDFT.__path__,
                                           prefix=BigDFT.__name__+'.',
                                           onerror=lambda x: None):
            try:
                # Import the module we found
                module = __import__(modname, fromlist="dummy")

                # See if it has an example routine.
                if hasattr(module, "_example"):
                    name = modname.replace(".", "_").replace("/", "")
                    if "Interop" in name:  # for now we disable interop
                        continue
                    dict["test_" + name] = gen_test(name, module)
            except ModuleNotFoundError:
                continue

        return type.__new__(mcs, name, bases, dict)


class ExampleTest(with_metaclass(MetaTest, TestCase)):
    """
    The actual test class. The MetaTest provides the test cases, and inheriting
    from TestCase lets it be automatically be run using ``unittest``.
    """
    pass
