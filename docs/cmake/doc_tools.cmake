# =======================================
# Macros and functions related to doxygen
# 
# =======================================

# ---------------------------------------------
# For a given component (numerics, kernel ...)
# get the list of directories containing
# sources to be used to generate doxygen doc.
# Result is saved in CACHE variable DOXY_INPUTS
# This variables is unique and stands for all
# components!
#
# NOT USEFUL ANYMORE ?
# ---------------------------------------------
# macro(update_doxygen_inputs COMP)
#   # Scan all dirs of current component and append
#   # them to DOXY_INPUTS/
#   # Do not include dirs matching 'test' and dirs listed
#   # in <COMP>_EXCLUDE_DOXY.
#   foreach(_dir IN LISTS ${COMP}_DIRS)
#     list(FIND ${COMP}_EXCLUDE_DOXY ${_dir} check_dir)
#     if(NOT ${_dir} MATCHES test AND ${check_dir} EQUAL -1)
#       list(APPEND DOXY_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
#     endif()
#   endforeach()
#   if(DOXY_INPUTS)
#     list(REMOVE_DUPLICATES DOXY_INPUTS)
#   endif()
#   # convert cmake list to a single string
#   #foreach(_dir ${DOXY_INPUTS})
#   #  set(_INPUTS "${_INPUTS} ${_dir}")
#   #endforeach()

#   # Save doxy_inputs to cache.
#   set(${COMP}_DOXYGEN_INPUTS ${DOXY_INPUTS} CACHE INTERNAL "List of inputs (directories) used by doxygen to generate doc for <COMP>.")
# endmacro()


# --------------------------------------------------
# Prepare config and targets for
# documentaition (xml-->rst-->sphinx)
#
# For a given component (numerics, kernel ...)
# generate rst files (for sphinx/breathe)
# from xml outputs (doxygen).
# 
# ---------------------------------------------
function(doxy2rst_sphinx COMPONENT)

  set(multiValueArgs HEADERS)
  cmake_parse_arguments(component "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  
  # --- Target : create rst files from xml outputs (doxygen) ---
  # run : make doxy2rst
  # depends : xml4rst
  #
  # Call python script to create rst files that can be parsed by sphinx/breathe
  # to create documentation.
  
  # output path, required to generate conf.py (breathe part)
  # This is the place where xml files created by doxygen
  # will be generated, for a later use by sphinx.
  set(DOXYGEN_4_RST ${DOXYGEN_OUTPUT}/xml4rst CACHE INTERNAL "doxy (xml) output path")
  
  # Doxygen conf for xml outputs for breathe. It might be different
  # from the one used for xml outputs for swig.
  file(MAKE_DIRECTORY ${DOXYGEN_4_RST}/${COMPONENT})

  include(doxycommon)
  set(DOXYGEN_GENERATE_HTML NO)
  set(DOXYGEN_GENERATE_XML YES)
  set(DOXYGEN_XML_OUTPUT xml4rst/${COMPONENT})
  doxygen_add_docs(
    ${COMPONENT}-doxy2xml ${CMAKE_SOURCE_DIR}/${COMPONENT}/src
    COMMENT "Generate xml/doxygen files for ${COMPONENT} (conf: ${DOXY_CONFIG_XML})."
    )
  
  if(WITH_${COMPONENT}_DOCUMENTATION)
    # Path where rst files will be generated.
    set(SPHINX_DIR "${CMAKE_BINARY_DIR}/docs/sphinx")
    
    # Create a new target used to create sphinx inputs (rst) from doxygen outputs (xml).
    # It calls a python function defined in gendoctools (create_breathe_files)
    add_custom_target(${COMPONENT}-xml2rst
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
      "from gendoctools.cpp2rst import create_breathe_files as f; f('${component_HEADERS}', '${CMAKE_SOURCE_DIR}', '${COMPONENT}', '${SPHINX_DIR}','${DOXYGEN_OUTPUT}/${DOXYGEN_XML_OUTPUT}')"
      VERBATIM
      DEPENDS ${COMPONENT}-doxy2xml
      )
    add_dependencies(rst_api ${COMPONENT}-xml2rst)
  
    add_custom_command(TARGET  ${COMPONENT}-xml2rst POST_BUILD
      COMMENT "${COMPONENT} : rst (c++ API) files have been generated in ${SPHINX_DIR}/reference/cpp.")
  endif()
endfunction()

# --------------------------------------
# Call this macro when configuration
# process is done for all components.
# 
# It will :
#  - create doxygen configs
# (for breathe, docstrings and so on)
#  - create targets related to
#   documentation
#    --> make doxygen-html to generate html doxygen doc from sources
#    --> make doxypng2sphinx to convert doxygen diagrams (class ...) to rst files
# --------------------------------------
macro(finalize_doc)
  if(WITH_DOCUMENTATION)

    # --- Target : html documentation from doxygen ---
    # run : make doxygen-html
    #  verbose mode : always off, since warnings
    # may be obtained with 'doxygen_warnings' target.
    # == set configuration file for doxygen doc from sources ==
    #  - Results in binary_dir/docs/
    #  - input config from config/doxy.config
    #  - only html output
    # 'inputs' are updated by each component, during call to update_doxygen_inputs macro.
    # Doc will be built when 'make doxygen' is called.
    # config file name

    # Build list of all dirs taken into accound by doxygen to build
    # doc, from each component own list.
    set(DOXYGEN_INPUTS)
    foreach(COMP IN LISTS COMPONENTS)
      list(APPEND DOXYGEN_INPUTS ${CMAKE_SOURCE_DIR}/${COMP}/src)
    endforeach()
    
    include(doxycommon)
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_HTML_OUTPUT ${DOC_ROOT_DIR}/html/doxygen)
    set(DOXYGEN_GENERATE_XML NO)
    doxygen_add_docs(
      doxygen-html ${DOXYGEN_INPUTS}
      COMMENT "Generate doxygen html doc ...")
    
    add_custom_command(TARGET doxygen-html POST_BUILD
      COMMENT "Doxygen documentation has been built in : \n - ${DOXYGEN_OUTPUT} (xml) \n - ${DOC_ROOT_DIR}/html/doxygen (html).")
    
    # --- Target : create class diagrams from doxygen outputs, for sphinx. ---
    # run : make doxypng2sphinx
    # depends : doxygen-html
    # Call a python function defined in gendoctools (find_doxygen_diagrams)
    add_custom_target(doxypng2sphinx
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
      "from gendoctools.generate_api import find_doxygen_diagrams as f; f('${CMAKE_BINARY_DIR}/docs/build/html/doxygen', '${CMAKE_BINARY_DIR}/docs/sphinx/reference')"
      VERBATIM
      DEPENDS doxygen-html
      COMMENT "Browse doxygen outputs (graphs ...)")

    add_custom_command(TARGET doxypng2sphinx POST_BUILD
      COMMENT "Done. Generate class_diagrams.rst file in ${CMAKE_BINARY_DIR}/docs/sphinx/reference")

    add_dependencies(html doxypng2sphinx)

    foreach(COMP ${COMPONENTS})
      add_dependencies(html ${COMP}-xml2rst)
    endforeach()

    # --- Generates conf.py, to describe sphinx setup ---
    # !! Should be call after doxygen setup
    # to have a correct DOXYGEN_INPUT value.
    configure_file (
      "${CMAKE_SOURCE_DIR}/docs/sphinx/conf.py.in"
      "${CMAKE_BINARY_DIR}/docs/sphinx/conf.py" @ONLY)

  endif()
endmacro()


# create a target to generate sphinx (rst) documentation
# from docstrings in python files (sphinx autodoc stuff)
# Note that docstrings are automatically generated using -doxygen option from swig.
function(docstrings2rst)
  set(oneValueArgs PATH NAME)

  cmake_parse_arguments(module "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(pymodule_fullname ${SICONOS_SWIG_ROOT_DIR}/${module_PATH}/${module_NAME}.py)
  
  message("Start setup for generation of rst files from python docstrings for ${pymodule_fullname}")
  
  # Path where rst files (docstrings --> rst) will be written.
  set(SPHINX_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs/sphinx/)
  
  # A new target to convert python/swig docstrings (swig outputs) into rst files (sphinx inputs)
  # Calls a python function defined in gendoctools (module_docstrings2rst)
  # --> make <comp>_autodoc
  add_custom_target(${module_NAME}_autodoc
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share:${CMAKE_BINARY_DIR}/wrap ${PYTHON_EXECUTABLE} -c
    "from gendoctools.python2rst import docstrings2rst as f; f('${module_PATH}', '${module_NAME}', '${SPHINX_OUTPUT_DIR}')"
    VERBATIM
    DEPENDS ${SWIG_MODULE_${module_NAME}_REAL_NAME}
    COMMENT "Create rst files from python docstrings for module siconos..${module_NAME}")

  # rst_api needs autodoc.
  add_dependencies(rst_api ${module_NAME}_autodoc)
  
endfunction()
  
