.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _siconos_doc:

=======================================================
 Writing and building the Siconos documentation
=======================================================

.. _about_doc:

Documentation overview
======================

The whole Siconos documentation is gathered on the Siconos web page, https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/.
It consists of

* several **"textbooks"**: details on algorithms, how to install and use the platform... everything users and developers should know about the Siconos software:

  * Getting and installing Siconos software
  * Users' guide
  * Developers' guide

  All these guides are generated with Sphinx (http://www.sphinx-doc.org/en/master/index.html) from files written in reStructuredText (http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

  The source files and the Sphinx configuration file are all in the docs/sphinx directory.

* the documentation of the **Siconos APIs** (Python and C++), generated automatically from inline comments in the source files (.h and .hpp).

  The rules and best practices for writing comments that give proper HTML documentation are detailed in :ref:`documenting_source_code`.



* in addition to these pages, "inline" documentation of the Python API based on docstrings, as usual for Python packages, e.g.:

  .. code-block:: python

   >>> import siconos.modeling as sm
   >>> help(sm.LagrangianDS)
   Help on ...
    ...



The complete documentation can be generated in one go with the doc target:

.. code-block:: bash

   cd build-dir
   cmake -DWITH_DOCUMENTATION=ON path_to_siconos_sources
   make
   make doc


This generates HTML pages in build-dir/docs/build/html.
The Siconos web site (http://siconos.org) is the online version of these pages.

The documentation build process is described in :ref:`build_doc`.



.. _doc_rules:

How to write Siconos documentation
==================================

Writing textbooks
-----------------

* Textbooks are generated from the rst files of the source-dir/docs/sphinx directory.
* Everything you need to know to write the Siconos docs is in the reStructuredText Primer: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html.
* All rst files are copied from the source directory to the binary directory by cmake, so only edit
  the files in the source directory.
* All figures and images must be in the sphinx/figures directory.

.. _documenting_source_code:

Document source code
--------------------

Since source code comments (in C/C++ header files) are used to generate the Doxygen, Sphinx and docstring documentation,
these comments must be written carefully, following the rules below.


.. rubric:: General rules

* Document all header files with Doxygen comments, as defined in http://www.stack.nl/~dimitri/doxygen/manual/index.html
* Try to follow the numpydoc (https://numpydoc.readthedocs.io/en/latest/) requirements.




.. rubric:: File description

Each header file must contain something like

.. code-block:: c++

     /*! @file SimpleMatrix.hpp
       Brief (no more than one line) description of the content of the file
     */


The name and description are used in the API content listings.
If this block is missing, the file (and all the objects or functions it contains) will not appear in the documentation.

.. rubric:: Classes and structs

Document each class like this

.. code-block:: c++

     /**
        @brief Short description of the class

        Detailed description
        equations (see the details about LaTeX below), references to textbook chapters and so on
     */
     class SiconosVector
     ...

.. rubric:: Class methods or functions

.. code-block:: c++

     /** @brief brief description
      *
      *  @param name_of_param1 description of the param
      *  @param name_of_param2 description of the param
      *  @return description of what is returned
     */
     double some_function(int p, int v)

There is no need to repeat the parameter types in the comments (param or return): they are extracted from the function prototype.
Something like

.. code-block:: c++

   /** @brief get size of A
   *   @param A double A
   *   @return unsigned int
   */
   unsigned int size(double * A) const;

is useless.


.. warning::

   Spaces, tabs and newlines matter, especially when generating docstrings and documentation from docstrings!

   - always add a newline before the first @param or @return;
   - align everything with the first line of the Doxygen comment.


Correct:

.. code-block:: c++

    /** @brief brief description
     *
     *  @param name_of_param1 description of the param
     *  @param name_of_param2 description of the param
     *  @return description of what is returned
     */

Wrong (missing newline and wrong indent):

.. code-block:: c++

    /** @brief brief description
     * @param name_of_param1 description of the param
     * @param name_of_param2 description of the param
     * @return description of what is returned
     */



.. rubric:: rst inside Doxygen comments

Use the "\\rst" / "\\endrst" tags to write reStructuredText (reST) specific comments (i.e. that Doxygen cannot handle).
See the details below for references and math formulas.

For comments with a leading asterisk, use the "\\rststar" / "\\endrststar" tags.

Note: most of the time, this breaks the rendering of the Python API.

.. rubric:: Enums, unions...

Since they will probably appear as global variables in the Python API,
each member of an enum must have an explicit comment, e.g.:

.. code-block:: cpp

    /** Global description of the enum */
    enum UBLAS_TYPE
    {
     /** id for dense matrix or vector */
     DENSE = 1,
     /** id for triangular matrix */
     TRIANGULAR,
    }


References to Sphinx documents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To refer to any other Sphinx document (reminder about Sphinx cross-references: http://www.sphinx-doc.org/en/stable/markup/inline.html),
use the "\\rst" / "\\endrst" tags:

.. code :: rst

  /** Class used to define friction-contact problems

  This class deals with blabla

  \rst

   See :ref:`global_fc_problem`

  \endrst

  */

or, with a leading asterisk:

.. code :: rst

  /** Class used to define friction-contact problems
   *
   * This class deals with blabla
   *
   * \rststar
   *
   *   See :ref:`global_fc_problem`
   *
   * \endrststar
   *
   *
   */



Math and LaTeX
~~~~~~~~~~~~~~

* inline math

  .. code:: rst

     use this \f$ \alpha \f$ to write inline math

     Important: always put a space before and after \f$

* displayed math

  - Wrap your formula between the "\f[" and "\f]" tags and write the math as you would in LaTeX.

.. _build_doc:

Building process
================

One target generates the whole documentation:

.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON ...
   make doc # The whole doc

The resulting files are in docs/build/html of the build directory.

Below are some details about the documentation generation process, useful only if you want to generate part of the doc or change the configuration and the process.


Tools, config and description
-----------------------------

.. rubric:: Tools

* `Doxygen`_: generates XML files from the comments in the source code (C/C++ headers);
* `Breathe`_: an extension to reStructuredText and Sphinx that reads and renders the Doxygen XML output.
  It converts the Doxygen XML output to rst files;
* `SphinxApidoc`_: generates rst files from the docstrings of the Python code (pure Python or from pybind11);
* `Sphinx`_: a powerful documentation generator (mostly for Python). It generates the HTML output from the rst files.


See the file `docs/requirements.txt <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/blob/main/docs/requirements.txt?ref_type=heads>`_ in the Siconos repository for the list of required Python packages.

.. rubric:: Config and sources:

* docs/CMakeLists.txt: main driver
* cmake/doc_tools.cmake: cmake macros and functions calling Doxygen, Sphinx and the other documentation tools
* cmake/doxygen_warning.cmake: included in LibraryProjectSetup, rules to build the "\*.warnings" files
* docs/sphinx/conf.py.in: main Sphinx configuration file
* docs/sphinx/index.rst.in: source of the main documentation page
* docs/sphinx/\*/\*.rst: inputs of the Sphinx doc (textbooks)
* docs/sphinx/figures/\*: all the figures used in the Sphinx doc
* docs/gendoctools/\*: Python tools used to generate the docs. This Python package is installed in <CMAKE_BINARY_DIR>/share
  at build time.

.. rubric:: Process

The whole process is run by "make doc" and is summarized in the figure below.

.. figure:: /figures/doc_process/build_breathe.*
   :figclass: align-center

   Generation of the documentation for the C++ and Python API


Remark: during the generation process, the Siconos Python packages are imported, and only
objects with non-empty docstrings are documented.

.. rubric:: Other (exotic) configuration options

* Full Doxygen HTML documentation (i.e. extract everything from the sources).

  Target: developers only. Everything is included in the Sphinx doc except things like call graphs or collaboration graphs.
  Use this when you need these graphs.

  Same as above, with an extra cmake option:

.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON -DUSE_DEVEL_DOXYGEN=ON
   make doc # The whole doc

The Doxygen output is "quiet", without warnings.
If you need the warnings (development), use:


* Generate Doxygen warnings

  Use the WITH_DOXYGEN_WARNINGS option (in USER_OPTIONS_FILE or on the command line), e.g.:

  .. code-block:: bash

     cmake -DWITH_DOCUMENTATION=ON -DWITH_DOXYGEN_WARNINGS=ON
     make filter_warnings

  This generates (during compilation) and prints the Doxygen warnings in files
  saved in CMAKE_BINARY_DIR/doxygen_warnings. A warnings file is generated for
  each input source file. The final call to 'make filter_warnings' concatenates all the relevant
  warnings into one file, doxygen_warnings/SUMMARY.warnings.

  The configuration of the Doxygen warnings is defined in docs/config/doxy_warnings.config.in and set up in
  cmake/doxygen_warnings.cmake.

.. code-block:: bash

   python3 -m pip install -rU ./docs/requirements.txt

See also the file `ci_gitlab/make_siconos_doc.sh <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/devel/siconos/-/blob/main/ci_gitlab/make_siconos_doc.sh?ref_type=heads>`_ in the Siconos repository.
It is used by the GitLab continuous integration to install all the dependencies needed to build the doc on Ubuntu, so it can help you set up your own doc build.



.. _Doxygen : http://www.stack.nl/~dimitri/doxygen/

.. _Sphinx : http://www.sphinx-doc.org/en/master/

.. _Breathe : https://github.com/michaeljones/breathe

.. _SphinxApidoc : https://www.sphinx-doc.org/en/master/man/sphinx-apidoc.html
