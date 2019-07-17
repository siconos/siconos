# ============================================================
# Prepare  config to generate siconos documentation :
#
# - automatic doxygen documentation, from source code
# - sphinx documentation for manuals
# 
# ============================================================
include(FindPythonModule)
include(FindSphinxModule)

# Devel option : generate pure doxygen doc (no sphinx), extracting all and with graphs.
option(USE_DEVEL_DOXYGEN "Build doxygen (html, latex ...) doc, extracting all. Devel purpose. Default = OFF" OFF)


# -- Doc config --
# Name and create build path for documentation
set(DOC_ROOT_DIR ${CMAKE_BINARY_DIR}/docs/build CACHE INTERNAL "Root path to generated documentation.")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/html)
# -- Gendoctools install --
# gendoctools is installed in binary dir as a python package
# used to generate docstrings, rst files and so on.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gendoctools/__init__.py.in ${CMAKE_BINARY_DIR}/share/gendoctools/__init__.py @ONLY)
file(GLOB files_py RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/gendoctools/*.py)
foreach(fname IN LISTS files_py)
  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/${fname}"
    "${CMAKE_BINARY_DIR}/share/${fname}" COPYONLY)
endforeach()
if(WITH_DOXY2SWIG)
  configure_file(${CMAKE_SOURCE_DIR}/externals/swig/doxy2swig.py ${CMAKE_BINARY_DIR}/share/gendoctools/doxy2swig.py COPYONLY)
endif()

# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/replace_py_escape.sh ${CMAKE_BINARY_DIR}/share/replace_py_escape.sh COPYONLY)

#  ====== doxygen documentation setup ======
# This setup is required if one at least of
# the following is needed:
# - doxygen warnings (see WITH_DOXYGEN_WARNINGS)
# - python docstrings from doxygen (xml) (see WITH_DOXY2SWIG)
# - doxygen documentation generated from source files (see WITH_DOCUMENTATION)
if(USE_DOXYGEN)

  # Search doxygen and set config
  find_package(Doxygen REQUIRED)
  # Output path for doxygen documentation
  if(NOT DOXYGEN_OUTPUT)
    set(DOXYGEN_OUTPUT ${DOC_ROOT_DIR}/doxygen CACHE INTERNAL "Path to doxygen outputs.")
  endif()
  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT}/man)
  
  # --- Config files and their options. --
  # --> a lot for doxygen warnings (from doxy_warnings.config.in)
  # --> one for documentation (html and xml to sphinx) (from doxy.config.in)

  # == set configuration file for doxygen warnings ==
  #  - Results in binary_dir/doxygen_warnings
  #  - input config from config/doxy_warnings.config
  #  - only xml output.
  if(WITH_DOXYGEN_WARNINGS)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings/man)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings/xml)
    set(DOXY_WARNINGS_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/config/doxy_warnings.config.in CACHE INTERNAL "config file for doxygen warnings output.")
    # All the configure process will be done for each component
    # thanks to 'include(doxygen_warnings)'. See LibraryProjectSetup.

    # A target to filter out doxygen warnings : remove empty files and concatenate the others.
    # call python function from gendoctools (filter_doxygen_warnings_files).
    add_custom_target(filter_warnings
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
      "from gendoctools.filters import filter_doxygen_warnings_files as f; f('${CMAKE_BINARY_DIR}/doxygen_warnings', 'SUMMARY.warnings')"
      VERBATIM
      COMMENT "Filter doxygen warnings (result : SUMMARY.warnings)."
      DEPENDS ${last_component} # No warnings if siconos has not been compiled. Should we remove this dep ? 
      )
  endif()
endif()

  
#  ====== Create a Target to generate documentation ======
if(WITH_DOCUMENTATION)

  # --- Search dependencies required to build documentation ---
  message("--------------------------------------------------")
  message("------ Start global setup for documentation ------\n")
  find_python_module(breathe REQUIRED)
  find_python_module(sphinx REQUIRED)
  find_sphinx_module(sphinxcontrib bibtex REQUIRED)
  find_sphinx_module(sphinxcontrib youtube)

  if(NOT ${youtube_FOUND})
    # This python package can not be found using standard ways such as conda or pip,
    # so we force download and install.
    message("Can not found sphinxcontrib-youtube. It will be downloaded and installed ...")
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -m pip install git+https://github.com/sphinx-contrib/youtube.git)
    find_sphinx_module(sphinxcontrib youtube REQUIRED)
  endif()
  # Sphinx theme for documentation
  find_python_module(sphinx_bootstrap_theme REQUIRED)

  # --- copy/conf sphinx files into build directory ---

  # main page
  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/index.rst.in"
    "${CMAKE_CURRENT_BINARY_DIR}/sphinx/index.rst" @ONLY)

  # siconos biblio
  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/siconos.bib"
    "${CMAKE_CURRENT_BINARY_DIR}/sphinx/siconos.bib" COPYONLY)

  # List and copy all rst files from source dir to binary dir
  file(GLOB_RECURSE files_rst RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx
    ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/*.rst)
  foreach(filepath IN LISTS files_rst)
    #get_filename_component(filename ${filepath} NAME)
    
    configure_file (
      "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/${filepath}"
      "${CMAKE_CURRENT_BINARY_DIR}/sphinx/${filepath}" COPYONLY)
  endforeach()
  
  # css
  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/_static/biblio.css"
    "${CMAKE_CURRENT_BINARY_DIR}/sphinx/_static/biblio.css" COPYONLY)

  # !! conf.py.in must be copied later, after doxygen setup !!
  # Done during finalize_doc call.
  
  # --- Build a target for sphinx ---
  # --> use make html or make latex ... to build documentation.
  set(SPHINX_EXECUTABLE sphinx-build)
  set(SPHINX_LOG_FILE ${CMAKE_BINARY_DIR}/docs/sphinx_warnings.log)    
  set(SPHINX_PARAMETERS
    -w ${SPHINX_LOG_FILE}
    -b html # html output
    -d build/doctrees # path to doctree files
    ${CMAKE_CURRENT_BINARY_DIR}/sphinx  # path to rst source
    )

  # - Target to copy figures required by sphinx from source to binary dir -
  #  make copy-sphinx-figures
  add_custom_target(copy-sphinx-figures
    COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/figures ${CMAKE_CURRENT_BINARY_DIR}/sphinx/figures
    )

  # -- Target to generate rst files from sources (C/C++/python) --
  #  1. creates rst files in a breathe-ready format from xml outputs of doxygen run on C/C++ headers
  #  2. creates rst files from 'autodoc' outputs, docstrings, or others python doc. things.
  add_custom_target(rst_api
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
    "from gendoctools.generate_api import build_rst_api as f; f('${CMAKE_CURRENT_BINARY_DIR}/sphinx')"
    VERBATIM
    COMMENT "Create main rst files for C/C++ and python API.")
  
  # -- To build html doc : make html --
  add_custom_target(html
    COMMAND ${SPHINX_EXECUTABLE} ${SPHINX_PARAMETERS} 2> errors.txt ${CMAKE_CURRENT_BINARY_DIR}/build/html
    VERBATIM
    DEPENDS copy-sphinx-figures rst_api
    )

  # post command to make html : filter sphinx warnings to a log file.
  add_custom_command(TARGET html POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
    "from gendoctools.filters import filter_sphinx_warnings as f; f('${SPHINX_LOG_FILE}')"
    VERBATIM
    COMMENT "Filtering sphinx logfile into ${SPHINX_LOG_FILE}.")
  
  add_custom_command(TARGET html POST_BUILD
    COMMENT "Done. Sphinx doc has been generated in ${DOC_ROOT_DIR}/html")
  
  # -- To build latex doc : make latex -- 
  set(SPHINX_LATEX_PARAMETERS
    -b latex # latex output
    -d build/doctrees # path to doctree files
    ${CMAKE_CURRENT_BINARY_DIR}/sphinx  # path to rst source
    )
  
  add_custom_target(latex
    COMMAND ${SPHINX_EXECUTABLE} ${SPHINX_LATEX_PARAMETERS}
    ${CMAKE_CURRENT_BINARY_DIR}/build/latex
    )

  add_custom_target(doc DEPENDS html)

  message("------ End global setup for documentation ------\n")


endif()
