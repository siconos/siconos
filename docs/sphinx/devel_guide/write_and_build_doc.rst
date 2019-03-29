.. _siconos_doc:

=======================================================
 Writing and building documention for Siconos Software
=======================================================


.. contents::
   :local:


.. _about_doc:

Documentation overview
======================

The whole Siconos documentation is gathered on Siconos webpage, https://siconos.gforge.inria.fr.
It consists in

* Different **"textbooks"** : details on algorithms, how to install, use the platform, ... everything user and developers should know concerning Siconos software.

  * Getting and installing Siconos software
  * Users' guide
  * Developers' guide
  
  All these guides are generated with sphinx (http://www.sphinx-doc.org/en/master/index.html) from files written in reStructuredText (http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)

  Source files and sphinx configuration file are all in docs/sphinx directory.

* **Siconos APIs** (Python and C++) documentation, automatically generated from inline comments in source files ( *.h/ *.hpp).
  Rules and best practices to write comments leading to a proper html documentation are detailed in :ref:`documenting_source_code`.



* In addition to these pages, "inline" documentation for the python API, as usual for python packages, based on docstrings, e.g. :

  .. code-block:: python

   >>> import siconos.kernel as sk
   >>> help(sk.SimpleMatrix)
   Help on class SimpleMatrix in module siconos.kernel:

   class SimpleMatrix(SiconosMatrix)
   |  Matrix (embedded various types of Boost matrices of double)
   |
   |  SimpleMatrix is used in the platform to store matrices (mathematical object) of
   |  double.
    ...

    
    
The complete documentation can be generated in one shot using doc target

.. code-block:: bash

   cd build-dir
   cmake -DWITH_DOCUMENTATION=ON path_to_siconos_sources
   make doc


It results in html pages, generated in build-dir/docs/build/html.
Siconos web site (https://siconos.gforge.inria.fr) is the online version of those pages.
   
The building process for documentation is described in :ref:`build_doc`.



.. _doc_rules:

How to write Siconos documentation
==================================

Writing textbooks
-----------------

* Textbooks are generated from rst files from source-dir/docs/sphinx directory.
* Everything you need to know to write siconos docs in reStructuredTextPrimer is available here : http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html.
* All rst files are copied from source dir to binary dir by cmake, so edit only
  files in source dir.
* All figures and images must be in sphinx/figures directory.

.. _documenting_source_code:

Document source code
--------------------

Since source code comments (in c/c++ header files) are used to generate doxygen, sphinx and docstrings documentation,
one has to be very careful when writing those comments and follow the rules below.


.. rubric:: General rules

* Document all header files using doxygen comments, as defined in http://www.stack.nl/~dimitri/doxygen/manual/index.html
* Do not comment doxygen comments --> breaks doxy2swig outputs.
    e.g. commenting a function and its doc will append the doc to the next function in the file
    and so break doxy2swig outputs
* Try to follow numpydoc (https://numpydoc.readthedocs.io/en/latest/) requirements.
    
  


.. rubric:: files description

Each header file must contain something like 

.. code-block:: c++

     /*! \file SimpleMatrix.hpp
       Brief (no more than one line) description of the content of the file
     */


The name and description will be used in the API contents listings.
If this block is not present, the file (and all the objects or functions it contents) won't appear in the documentation.

.. rubric:: Classes and structs

Document each class like this

.. code-block:: c++

     /** Short description of the class
     *
     * Detailed description
     * equations (see details about latex below), reference to textbooks chapter and so on
     */
     class SiconosVector
     ...

.. rubric:: Class methods or functions

.. code-block:: c++

     /** brief description
      * \param name_of_param1 description of the param
      * \param name_of_param2 description of the param
      * \return description of what is returned
     */
     double some_function(int p, int v)

No need to repeat parameters types in comments (param or return)! They will be extracted from function prototype.
Something like

.. code-block:: c++

   /** get size of A
   * \param A double A
   *  \return unsigned int
   */
   unsigned int size(double * A) const;

is totally useless ...

.. rubric:: rst inside doxygen commments

Use  "\\rst" / "\\endrst" tags to write reStructuredText (reST) specific (i.e. that doxygen can not tackle) comments.
See details below for references and math formula.

In the case of comments with leading asterisk, use "\\rststar" / "\\endrststar" tags

.. rubric:: Enums, union ...

Since they will probably appear as global variables in python API,
it's important that each component of the enum has an explicit comment, e.g:

.. code-block:: cpp

    /** Global description of the enum */ 
    enum UBLAS_TYPE
    {
     /** id for dense matrix or vector */
     DENSE = 1, 
     /** id for triangular matrix */
     TRIANGULAR,
    }


References to sphinx documents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To refer to any other sphinx document (reminder about sphinx cross-ref : http://www.sphinx-doc.org/en/stable/markup/inline.html)
use "\\rst" / "\\endrst" tags :

.. code :: rst

  /** Class used to defined friction-contact problems
  
  This class deals with blabla

  \rst
  
   See :ref:`global_fc_problem`

  \endrst

  */

or with leading asterisk

.. code :: rst

  /** Class used to defined friction-contact problems
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

  

Math and latex
~~~~~~~~~~~~~~

* inline math

  .. code:: rst

     use this \f$\alpha\f$ to write inline math

* displayed math

  - Wrap your formula between "\rst" and "\endrst" tags and write math as you would with sphinx (see http://www.sphinx-doc.org/en/master/ext/math.html).
  - Between rst tags, replace all occurences of :math:'\'dot (one backlash)  with :math:'\\'dot (two backlashes), else doxygen will fail to produce documentation.
  
  A simple example :

  .. code:: rst

     \rst
     
     .. math::
      
        y &=& h(X,t,\lambda,Z)\\
        R &=& g(X,t,\lambda,Z)

     \endrst

  * New line after math keyword is required.
  * Indentation for formula (related to math keyword) is required.
  
  For more complicated maths, use nowrap keyword :
  
  .. code:: rst

     \rst
     
     .. math::
        :nowrap:
      
         \left\{\begin{array}{l}
         y \geq 0, \lambda \geq 0, y^{T} \lambda=0\\
         if y \leq 0 \quad \mbox{then} \quad \\dot y(t^{+}) - e \\dot y(t^{-}) \geq 0, \quad  \lambda \geq 0, (\\dot y(t^{+}) - e \\dot y(t^{-}))^{T} \lambda=0
         \end{array}\right.

     \endrst

     
If you need comments with leading asterisk, use "\rststar" / "\endrststar" tags :

.. code:: rst
   
 * \rststar
 *
 * .. math::
 *    :nowrap:
 *
 *    \begin{eqnarray}
 *    \begin{cases}
 *     M v =  q +  H r \\
 *     u = H^\top v + b \\
 *     \hat u = u +\left[
 *       \left[\begin{array}{c}
 *           \mu^\alpha \|u^\alpha_{T}\|\\
 *           0 \\
 *           0
 *         \end{array}\right]^T, \alpha = 1 \ldots n_c
 *      \right]^T \\ \\
 *      C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu}
 *     \end{cases}
 *    \end{eqnarray}
 *
 * \endrststar

  

  
.. _build_doc:

Building process
================

One target to generate the whole documentation : 
  
.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON ...
   make doc # The whole doc

Resulting files will be in docs/build/html of the current building path.

Below are some details about the documentation generation process, useful only if you want to generate a subpart of the doc or change the configuration and the process.


Tools, config and description
-----------------------------

**Tools:**

* `Doxygen`_ : tool able to generate documentation from annotated C++ sources, in html, xml ...
* `Sphinx`_ : powerful generator of documentation (mostly for Python)
* `Breathe`_ : an extension to reStructuredText and Sphinx to be able to read and render the Doxygen xml outputs.
* `Doxy2swig`_ : converter from doxygen XML to SWIG docstring.


Images are sometimes better than words : the different operations are  detailed on figures below

.. rubric:: Generation of rst files for C++ API


How does it work?

`Doxygen`_ is used to generate xml files from comments in headers of each Siconos component. Python scripts are used to postprocess those xml files and produce rst files fitting with `Breathe`_ requirements and ready for `Sphinx_`.


*Config and sources:*

* cmake/doc_tools.cmake : cmake macros and functions calling doxygen, sphinx or other tool related to documentation.
* docs/gendoctools/* : python tools used to generate docs. This python package will be installed in <CMAKE_BINARY_DIR>/share
  at build time.
* docs/config/doxyxml2sphinx.config.in : doxygen (xml output) for breathe/sphinx

.. figure:: /figures/doc_process/build_breathe.*
   :figclass: align-center

   Generation of rst files for C++ API


.. rubric:: Generation of rst files for Python API

How does it work?

Doxygen generates xml from comments in headers. Some python scripts are
used to postprocess those xml files and produce .i files (swig), ending in
docstrings in generated swig python modules.
We have written a wrapper to doxy2swig (https://github.com/m7thon/doxy2swig) to fit with our needs.
Finally, rst files are generated, based on those docstrings, in autodoc format, for sphinx.

*Config and sources:*

* cmake/swig_python_tools.cmake : python functions used to drive docstrings
  generation
* docs/gendoctools/* : python tools used to generate docs. This python package will be installed in <CMAKE_BINARY_DIR>/share
  at build time.
* docs/config/doxy2swig.config.in : doxygen (xml output) config, for swig and docstrings


.. figure:: /figures/doc_process/build_doxy2swig.*
   :figclass: align-center

   Generation of rst files for Python API

Remark : during generation process, siconos python packages are imported and only
objects with non-empty docstrings are documented. 
 

.. rubric:: html pages generation

How does it work?

All rst files (from source dir and generated for Python and C++ API) and processed by `Sphinx`_ to produce html documentation.
          
*Config and sources:*

* docs/CMakeLists.txt : main driver
* cmake/doc_tools.cmake : cmake macros and functions calling doxygen, sphinx or other tool related to documentation.
* cmake/doxygen_warning.cmake : included in LibraryProjectSetup, rules to build "\*.warnings" files.
* docs/sphinx/conf.py.in : main sphinx configuration file
* docs/sphinx/index.rst.in : source for documentation main page
* docs/sphinx/\*/\*.rst : inputs for sphinx doc (textbooks)
* docs/sphinx/figures/\* : all figures used in sphinx doc
* docs/gendoctools/* : python tools used to generate docs.
* docs/config/doxy.config.in : doxygen (html output) config
* docs/config/doxy_warnings.config.in : doxygen (log output) config

            
.. figure:: /figures/doc_process/build_html_process.*
   :figclass: align-center
              
   make doc toolchain
            
.. figure:: /figures/doc_process/targets_dep.*
   :figclass: align-center

   make targets (related to doc) dependencies

           
.. rubric:: Other (exotic) configuration options

* Full doxygen (i.e. extract all from sources!)  html documentation.

  Target : developpers only. Indeed, everything is included in sphinx doc except things like callgraph or collaboration graphs.
  Use this when those graphs are needed. 

  Same as above, with a new cmake option :

.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON -DUSE_DEVEL_DOXYGEN=ON
   make doc # The whole doc

Doxygen ouput is set to be "quiet" and without warnings.
But, if required (devel), use:


* Generate doxygen warnings

  Use WITH_DOXYGEN_WARNINGS option (in USER_OPTIONS_FILE or in command line), e.g. :

  .. code-block:: bash

     cmake -DWITH_DOCUMENTATION=ON -DWITH_DOXYGEN_WARNINGS=ON
     make filter_warnings
     
  It will generate (during compilation process) and print doxygen warnings in files
  saved in CMAKE_BINARY_DIR/doxygen_warnings. A warnings file is generated for
  each input source file. The final call to 'make filter_warnings' will concatenate all interesting
  warnings into one file, doxygen_warnings/SUMMARY.warnings

  doxygen warnings conf is defined in docs/config/doxy_warnings.config.in and setup in
  cmake/doxygen_warnings.cmake.


* Generate docstrings

  .. code-block:: bash

     cmake -DWITH_DOCUMENTATION=ON -DWITH_DOXY2SWIG=ON
     make numerics_docstrings

  This option is set to ON by default.
     
  To produce documentation in python interface, xml outputs from doxygen are used to create swig files containing 'docstrings' for python.

  Comments written in C++ (doxygen) will be available in python interface, e.g. :

  .. code-block:: python

     import siconos.kernel as sk
     help(sk.DynamicalSystem)
   
     Help on class LCP in module siconos.kernel:

     class LCP(LinearOSNS)
     |  Non Smooth Problem Formalization and Simulation.
     |
     |  author: SICONOS Development Team - copyright INRIA
     |
     |  This is an abstract class, that provides an interface to define a non smooth
     |  problem:
     |
     |  *   a formulation (ie the way the problem is written)
     |  *   a solver (algorithm and solving formulation, that can be different from
     |      problem formulation)
     |  *   routines to compute the problem solution.


 
Dependencies
------------

* `Doxygen`_
* `Sphinx`_, sphinxcontrib.bibtex, sphinxcontrib-youtube, sphinxcontrib-napoleon
* sphinx-bootstrap-theme
* `Breathe`_

See docs/requirements.txt for a list of required python packages, and try for example

.. code-block:: bash

   pip install -r ./docs/requirements.txt
   pip install git+https://github.com/sphinx-contrib/youtube.git

See also the file CI/make_siconos_doc.sh that may be helpful to install siconos docs, since it is used by continuous integration on gitlab to provide all dependencies required to build doc on ubuntu. 



More about Doxygen to sphinx rst
--------------------------------

Some other tools to generate rst from doxygen have been tested : Exhale and doxyrest. We choose breathe, that seems more appropriate to our case. Exhale and doxyrest configs are kept for the records in siconos-junk/sandbox project.

Existing tools (as far as we know ...):

* Sphinx/Exhale(breathe) : https://github.com/svenevs/exhale`Sphinx/Exhale
* doxyrest https://github.com/vovkos/doxyrest
* https://bitbucket.org/trlandet/doxygen_to_sphinx_and_swig

.. _Doxygen : http://www.stack.nl/~dimitri/doxygen/

.. _Sphinx : http://www.sphinx-doc.org/en/master/

.. _Breathe : https://github.com/michaeljones/breathe

.. _Doxy2swig : https://github.com/m7thon/doxy2swig
