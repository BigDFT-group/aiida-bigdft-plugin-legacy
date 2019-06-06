How to write documentation in this package
==========================================

The syntax in the text pages is ``ReStructured Text`` format.
see `this page`_ for examples of the syntax, and `this one`__ for a more comprehensive explanation of
the syntax.

.. __: http://docutils.sourceforge.net/rst.html

.. _this page: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

The syntax chosen for the auto-documentation of the Python module is ``Google Docstring`` syntax.
see `the google Docstring example`__ in the `Sphinx-napoleon`__ package.

.. __: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html

.. __: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/index.html

How to generate documentation
-----------------------------

 * Install ``sphinx`` package and the requirements written in the file  ``source/requirements.txt`` of
   provided with the sources of ``PyBigDFT``.

 * Copy the quick build file ``Makefile-sphinxbuild`` of the ``PyBigDFT`` source directory into ``Makefile``

 * Run the command ``make html``

 * The documentation of the package can be explored starting from the file ``build/html/index.html``.

Automodule generation
---------------------

The directive to be used is mainly the ``automodule`` one, see `this`__ page:

.. __: http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
