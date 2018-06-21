.. _write_and_build_doc:

#####################################################
Writing and building documention for Siconos Software
#####################################################

* Rules for developpers when it comes to write proper documentation
* Details on the building process for documentation


.. _build_doc:

**********************
Building documentation
**********************

Which documentation for who?

* The whole web site, including user manual, tutorials, siconos api and so on.

Target : users and developpers.
  
  
.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON ...
   make doc # The whole doc

Resulting files will be in docs/build/html of the current building path.

Some other targets are available to build only sub-parts of the doc. Check :ref:`make_targets` for details.





* Full doxygen (i.e. extract all from sources!)  html documentation.

  Target : developpers only. Indeed, everything is included in sphinx doc except things like callgraph or collaboration graphs.
  Use this when those graphs are needed. 

  Same as above, with a new cmake option :

.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON -DUSE_DEVEL_DOXYGEN=ON
   make doc # The whole doc
  
* Other useful (devel) options:

  * USE_DOXYREST=ON, to generate rst files from xml outputs. Test purpose. Useful to produce rst files from "related pages" of doxygen.
  * USE_EXHALE=ON to generate rst files from xml outputs. On test. In the future, we should set this to ON all the time. Warning : combining this option with EXTRACT_ALL ON    (USE_DEVEL_DOXYGEN=ON) may result in a very long time to build documentation.

*****************************
Building process (doc) review
*****************************
    
Building process for docs is defined in :

* docs/CMakeLists.txt : main driver
* cmake/doxygen_tools.cmake
* cmake/doxygen_warning.cmake : included in LibraryProjectSetup, rules to build "\*.warnings" files.
* cmake/filedoxy.config.in

  
Dependencies
============

* doxygen
* sphinx, sphinxcontrib.bibtex
* doxyrest (optional)
* breathe (optional)

  
  
Doxygen docs
============

  
API doc from doxygen
--------------------

Automatic generation of documentation from source code (C/C++).

* **Sources** : headers
* **Config** : docs/config/doxy.config.in. This file is used to generate CMAKE_BINARY_DIR/docs/config/doxy.config that will be used by doxygen to build doc in DOXYGEN_OUTPUT.
* **Results** (all in in CMAKE_BINARY_DIR/doc/build/html/doxygen)
  * html files
  * png (class diagrams ...)

**Usage :**

.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON
   make doxygen
   # check : open in web browser CMAKE_BINARY_DIR/doc/build/html/doxygen/index.html


**Remark :**

Doxygen ouput is set to be "quiet" and without warnings.
But, if required (devel), use:

.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON -DWITH_DOXYGEN_WARNINGS=ON
   make filter_warnings
   # if WITH_DOXYGEN_WARNINGS_INFILE=ON, create doxygen_warnings/SUMMARY.warnings

It will generate (during compilation process) and print doxygen warnings, either on screen or in files
saved in CMAKE_BINARY_DIR/doxygen_warnings (if WITH_DOXYGEN_WARNINGS_INFILE=ON ...).

doxygen warnings conf is defined in docs/config/doxy_warnings.config.in and setup in
cmake/doxygen_warnings.cmake.


Docstrings (for swig)
---------------------

To produce documentation in python interface, xml outputs from doxygen are used to create swig
files containing 'docstrings' for python.

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


Usage 

.. code-block:: bash

   cmake -DWITH_DOCUMENTATION=ON -DWITH_DOXY2SWIG=ON
   make numerics_docstrings

   
Process :

#. Generates xml files for each component (doxygen).\n
   Config from docs/config/doxy2swig.config.in\n
   Results in CMAKE_BINARY_DIR/docs/build/doxygen/doxy2swig-xml/component_name\n
   target : make component_name_xml4swig
   
#. Generates swig files (.i) from xml for one component and concatenate into
   CMAKE_BINARY_DIR/wrap/siconos/component_name-docstrings.i. \n
   Tool = doxy2swig (https://github.com/m7thon/doxy2swig) saved in externals/swig.\n
   target : make component_name_docstrings

These \*-docstrings.i files are included into component.i (e.g. kernel.i) to produce doc during swig process.

Todo : test this tool (https://bitbucket.org/trlandet/doxygen_to_sphinx_and_swig/) which produces
both docstrings for swig and rst for sphinx from doxygen outputs, in one shot.



Sphinx doc
==========

* conf defined in docs/sphinx/conf.py.in, used to generate (cmake) CMAKE_BINARY_DIR/docs/sphinx/conf.py
  as sphinx main configuration file.

Process :
* generate (cmake) CMAKE_BINARY_DIR/docs/sphinx/conf.py from docs/sphinx/conf.py.in
==> sphinx main configuration file
* generate (cmake) CMAKE_BINARY_DIR/docs/sphinx/index.rst from docs/sphinx/index.rst.in
==> defines main page for the resulting doc (i.e. website home page)
* generate (copy only) CMAKE_BINARY_DIR/docs/sphinx/*.rst (recursive) from docs/sphinx/*.rst



Remarks
-------

* Great care is required when writing doxygen comments, doxy2swig, sphinx, exhale are very sensitive and fail easily ...

  * Do not comment doxygen comments --> breaks doxy2swig outputs.
    e.g. commenting a function and its doc will append the doc to the next function in the file
    and so break doxy2swig outputs
  * Write proper latex formulas in doxygen comments, else doxy2swig or sphinx output will fail.

  * Try to follow numpydoc (https://numpydoc.readthedocs.io/en/latest/) requirements.
    
  


.. _doc_rules:

***********
Writing doc
***********


Document source code
====================

General rules:

* Document all header files using doxygen comments, as defined in http://www.stack.nl/~dimitri/doxygen/manual/index.html
* Use  "\rst" / "\endrst" tags to write reStructuredText (reST) specific (i.e. that doxygen can not tackle) comments.
  See details below for references and math formula.
* In the case of comments with leading asterisk, use "\rststar" / "\endrststar" tags



References to sphinx documents
------------------------------

To refer to any other sphinx document (reminder about sphinx cross-ref : http://www.sphinx-doc.org/en/stable/markup/inline.html)
use "\rst" / "\endrst" tags :

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
--------------

* inline math

  .. code:: rst

     use this \f$\alpha\f$ to write inline math

* displayed math

  - Wrap your formula between "\rst" and "\endrst" tags and write math as you would with sphinx (see http://www.sphinx-doc.org/en/master/ext/math.html).
  - Between rst tags, replace all occurences of :code:'\dot' with :code:'\\dot', else doxygen will fail to produce documentation.
  
  
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

  
Enums, union ...
----------------

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

  


Doxygen to sphinx
=================
  
Existing tools (as far as we know ...)

* Sphinx/Exhale(breathe) : https://github.com/svenevs/exhale`Sphinx/Exhale
* doxyrest https://github.com/vovkos/doxyrest
* https://bitbucket.org/trlandet/doxygen_to_sphinx_and_swig

Both exhale and doxyrest are available in siconos, (use -DUSE_EXHALE=ON or -DUSE_DOXYREST=ON). We prefer exhale.


Exhale/sphinx pipeline :
* generates rst files from xml files (doxygen outputs) in CMAKE_BINARY_DIR/docs/sphinx/api
* build html (as usual) from rst files, in CMAKE_BINARY_DIR/docs/build/html/api

Exhale conf is defined in conf.py.in (sphinx) and may also handle doxygen run (xml outputs + rst generations from those outputs).

  
Doxyrest works the same way but is not as convenient as exhale. Outputs are in CMAKE_BINARY_DIR/docs/sphinx/from_doxygen.

Both (exhale and doxyrest) are quite slow and doc generation may take long time ...
It seems that it strongly depends on the chosen theme for sphinx (avoid bootswatch).
