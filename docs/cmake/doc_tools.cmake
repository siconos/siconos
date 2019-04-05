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
# ---------------------------------------------
macro(update_doxygen_inputs COMP)
  # Scan all dirs of current component and append
  # them to DOXY_INPUTS/
  # Do not include dirs matching 'test' and dirs listed
  # in <COMP>_EXCLUDE_DOXY.
  foreach(_dir IN LISTS ${COMP}_DIRS)
    list(FIND ${COMP}_EXCLUDE_DOXY ${_dir} check_dir)
    if(NOT ${_dir} MATCHES test AND ${check_dir} EQUAL -1)
      list(APPEND DOXY_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
    endif()
  endforeach()
  if(DOXY_INPUTS)
    list(REMOVE_DUPLICATES DOXY_INPUTS)
  endif()
  # convert cmake list to a single string
  foreach(_dir ${DOXY_INPUTS})
    set(_INPUTS "${_INPUTS} ${_dir}")
  endforeach()

  # Save doxy_inputs to cache.
  set(${COMP}_DOXYGEN_INPUTS ${_INPUTS} CACHE INTERNAL "List of inputs (directories) used by doxygen to generate doc for <COMP>.")
endmacro()


# --------------------------------------------------
# Prepare config and targets for
# documentaition (xml-->rst-->sphinx)
#
# For a given component (numerics, kernel ...)
# generate rst files (for sphinx/breathe)
# from xml outputs (doxygen).
# 
# ---------------------------------------------
macro(doxy2rst_sphinx COMP)

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
    set(DOXY_QUIET "YES")
    set(DOXY_WARNINGS "NO")
    set(GENERATE_HTML NO)
    set(GENERATE_XML YES)
    set(EXTRACT_ALL NO)
    set(EXTRACT_PRIVATE NO)
    set(XML_OUTPUT xml4rst/${COMP})
    file(MAKE_DIRECTORY ${DOXYGEN_4_RST}/${COMP})
    # Set config file name
    set(DOXY_CONFIG_XML "${CMAKE_BINARY_DIR}/docs/config/${COMP}doxy-xml.config")

    # Get list of inputs
    set(DOXYGEN_INPUTS ${${COMP}_DOXYGEN_INPUTS})
    
    configure_file(${CMAKE_SOURCE_DIR}/docs/config/doxyxml2sphinx.config.in ${DOXY_CONFIG_XML} @ONLY)

    # Create a new target used to create doxygen outputs (xml).
    # Required for documentation and/or for serialization.
    add_custom_target(${COMP}-doxy2xml
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG_XML}
      OUTPUT_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy4rst.log ERROR_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy4rst.log
      COMMENT "Generate xml/doxygen files for ${COMP} (conf: ${DOXY_CONFIG_XML}).") 

    if(WITH_${COMPONENT}_DOCUMENTATION)
      # Path where rst files will be generated.
      set(SPHINX_DIR "${CMAKE_BINARY_DIR}/docs/sphinx")
      
      # Create a new target used to create sphinx inputs (rst) from doxygen outputs (xml).
      # It calls a python function defined in gendoctools (create_breathe_files)
      add_custom_target(${COMP}-xml2rst
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
        "from gendoctools.cpp2rst import create_breathe_files as f; f('${${COMP}_HDRS}', '${CMAKE_SOURCE_DIR}', '${COMP}', '${SPHINX_DIR}','${DOXY_CONFIG_XML}')"
        VERBATIM
        DEPENDS ${COMP}-doxy2xml
        )
      
      add_custom_command(TARGET  ${COMP}-xml2rst POST_BUILD
        COMMENT "${COMP} : rst (c++ API) files have been generated in ${SPHINX_DIR}/reference/cpp.")
    endif()
endmacro()

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
    set(DOXY_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/config/doxy.config" CACHE INTERNAL "Doxygen configuration file : used to produce html (doxygen only) and xml files for sphinx (breathe).")
    # The generation of the config and the creation of the target will be done later,
    # after the update of inputs by each component, with macro 'finalize_doc'
    set(DOXY_QUIET "YES")
    set(DOXY_WARNINGS "NO")
    set(GENERATE_HTML YES)
    set(GENERATE_XML NO)
    set(EXTRACT_ALL NO)
    if(USE_DEVEL_DOXYGEN) # OFF  by default. Activate to extract all.
      set(EXTRACT_ALL YES)
    endif()
    # Build list of all dirs taken into accound by doxygen to build
    # doc, from each component own list.
    set(DOXYGEN_INPUTS)
    foreach(COMP IN LISTS COMPONENTS)
      set(DOXYGEN_INPUTS "${DOXYGEN_INPUTS} ${${COMP}_DOXYGEN_INPUTS}")
    endforeach()
    configure_file(${CMAKE_SOURCE_DIR}/docs/config/doxy.config.in ${DOXY_CONFIG} @ONLY)

    # A new target to create  html files (doxygen) from code sources.
    # One target for all components.
    add_custom_target(doxygen-html
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}
      COMMENT "Run doxygen ${DOXY_CONFIG} ...")
    
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



function(docstrings2rst module_path module_name)
  set(pymodule_name ${SICONOS_SWIG_ROOT_DIR}/${module_path}/${module_name}.py)
  
  message("Start setup for generation of rst files from python docstrings for ${pymodule_name}")
  # A target to postprocess latex forms in docstrings into
  # something readable by sphinx.
  # Calls a python function defined in gendoctools (replace_latex)
  add_custom_target(${module_name}_replace_latex
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
    "from gendoctools.common import replace_latex as f; f('${pymodule_name}', '${SICONOS_SWIG_ROOT_DIR}/tmp_${COMPONENT}/')"
    VERBATIM
    DEPENDS ${SWIG_MODULE_${module_name}_REAL_NAME}
    COMMENT "Insert latex into docstrings.")
  
  # Path where rst files (docstrings --> rst) will be written.
  set(SPHINX_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs/sphinx/)
  # python modules for previous components are required to apidoc (e.g. kernel.py for control).
  # So we get this last comp and add a dependency.
  list(APPEND PROCESSED_PYTHON_MODULES ${SWIG_MODULE_${module_name}_REAL_NAME})
  list(REMOVE_DUPLICATES PROCESSED_PYTHON_MODULES)
  set(PROCESSED_PYTHON_MODULES ${PROCESSED_PYTHON_MODULES} CACHE INTERNAL "python modules for siconos")
  
  # A new target to convert python/swig docstrings (swig outputs) into rst files (sphinx inputs)
  # Calls a python function defined in gendoctools (module_docstrings2rst)
  # --> make <comp>_autodoc
  add_custom_target(${module_name}_autodoc
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share:${CMAKE_BINARY_DIR}/wrap ${PYTHON_EXECUTABLE} -c
    "from gendoctools.python2rst import docstrings2rst as f; f('${COMPONENT}', '${module_path}', '${module_name}', '${SPHINX_OUTPUT_DIR}', '${SICONOS_SWIG_ROOT_DIR}')"
    VERBATIM
    DEPENDS ${module_name}_replace_latex
    COMMENT "Create rst files from python docstrings for module siconos.${module_name}")
  
  # Create dependency between autodoc target and siconos python modules.
  foreach(dep IN LISTS PROCESSED_PYTHON_MODULES)
    add_dependencies(${module_name}_autodoc ${dep})
  endforeach()
  
  # rst_api needs autodoc.
  add_dependencies(rst_api ${module_name}_autodoc)
  
endfunction()
  
