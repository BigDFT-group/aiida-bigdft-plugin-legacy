"""Handling of the input options

This module contains the useful quantities to deal with the preparation and
the usage of inputfiles for BigDFT code. The main object is the
:class:`Inputfile` class, which inherits from a python dictionary. Such
inheritance is made possible by the internal representation of the BigDFT
inputfile, which employs the YAML syntax. This means that there is a one-to-one
correspondence between a python dictionary and a BigDFT inputfile.

"""


class Inputfile(dict):
    """ The BigDFT inputfile.

    Principal object needed to run a  BigDFT calculation.
    Might be initialized either from a dictionary, a
    :py:class:`~BigDFT.Logfiles.Logfile` instance
    or a (yaml-compliant) filename path

    Note:

       Each of the actions of the module :py:mod:`~BigDFT.InputActions`, which
       is defined on a generic dictionary, also corresponds to a method of of
       the `Inputfile` class, and it is applied to the class instance.
       Therefore also the first argument of the corresponding action is
       implicitly the class instance. For the
       :py:func:`~BigDFT.InputActions.remove` method, the action has to be
       invoked should come from the :py:mod:`~BigDFT.InputActions` module.


    .. _input_action_example:
    Example:

       >>> import InputActions as A, Inputfiles as I
       >>> inp=I.Inputfile()
       >>> inp.set_hgrids(0.3) # equivalent to A.set_hgrids(inp,0.3)
       >>> inp
       {'dft': {'hgrids': 0.3}}
       >>> inp.optimize_geometry() # equivalent to A.optimize_geometry(inp)
       >>> inp
       {'dft': {'hgrids': 0.3},'geopt': {'method': 'FIRE',
                                         'ncount_cluster_x': 50} }
       >>> # equivalent to A.remove(inp,A.optimize_geometry)
       >>> inp.remove(A.optimize_geometry)
       >>> inp
       {'dft': {'hgrids': 0.3}}

     .. todo ::

         Consider the possiblity of initializing an `Inputfile` instance
         from a ``yaml`` file. And also from a
         :py:class:`~BigDFT.Logfiles.Logfile` class


    """

    def __init__(self, *args, **kwargs):
        import BigDFT.InputActions as A
        dict.__init__(self, *args, **kwargs)
        functions = dir(A)
        for action in functions:
            if "__" in action:
                continue
            from functools import partial
            func = getattr(A, action)
            setattr(self, action, partial(func, self))
