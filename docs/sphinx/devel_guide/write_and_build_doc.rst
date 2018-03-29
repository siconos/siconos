.. _write_and_build_doc:

#####################################################
Writing and building documention for Siconos Software
#####################################################

* Rules for developpers when it comes to write proper documentation
* Details on the building process for documentation


.. _build_doc:

Building documentation
**********************

Which documentation for who?

* The whole web site, including user manual, tutorials, siconos api and so on.

Target : users and developpers.
  
  
.. code:: bash

   cmake -DWITH_DOCUMENTATION=ON ...
   make doc # The whole doc
   # or
   make doxygen (doxygen only ...)
   # or
   make html (html only)

Resulting files will be in docs/build/html of the current building path.


* Full doxygen (i.e. extract all from sources!)  html documentation.
  Target : developpers only. Indeed, everything is included in sphinx doc except things like callgraph or collaboration graphs.
  Use this when those graphs are needed. 

  Same as above, with a new cmake option :

.. code:: bash

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

3 differents parts :

* API doc generation : standard doxygen doc, from C/C++ headers.
* Doxygen warnings
* docstring generation : create docstrings for swig outputs, based on doxygen.

API doc from doxygen
--------------------
  
Defined in file docs/config/doxy.config.in.
During cmake process this file is used to generate CMAKE_BINARY_DIR/docs/config/doxy.config that will
be used by doxygen to build doc in DOXYGEN_OUTPUT.

html, latex, xml and man outputs

man? Does it works?


Doxygen warnings
----------------
  What for?

  Process:

* If WITH_DOXYGEN_WARNINGS is ON (default=OFF)
  
  * Conf generated from docs/config/doxy_warnings.config.in
  * Work done following include(doxygen_warnings) for each component (numerics, kernel ...)
    called in LibraryProjectSetup.
  * Result : .config and .warnings files for each header of the current component.
    All files generated in CMAKE_BINARY_DIR/doxygen_warnings.

Doxygen to swig
---------------

Describe this ...


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





.. _doc_rules:

Writing doc
***********


* Document all sources files (headers) using doxygen, as defined in http://www.stack.nl/~dimitri/doxygen/manual/index.html


  
  

