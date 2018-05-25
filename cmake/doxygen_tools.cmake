# =======================================
# Macros and functions related to doxygen
# 
# =======================================

# ---------------------------------------------
# For a given component (numerics, kernel ...)
# get the list of directories containing
# sources to be used to generate doxygen doc.
# Result is saved in CACHE variable DOXY_INPUTS
# ---------------------------------------------
macro(update_doxy_config_file COMP)
  # Scan all dirs of current component and append
  # them to DOXY_INPUTS
  foreach(_dir ${${COMP}_DIRS})
    list(FIND ${COMP}_EXCLUDE_DOXY ${_dir} check_dir)
    if(NOT ${_dir} MATCHES test AND ${check_dir} EQUAL -1)
      list(APPEND DOXY_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
    endif()
  endforeach()
  list(REMOVE_DUPLICATES DOXY_INPUTS)
  # Save doxy_inputs to cache.
  set(DOXY_INPUTS ${DOXY_INPUTS} CACHE INTERNAL "doxy inputs")
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
    set(DOXYGEN_4_RST ${DOXYGEN_OUTPUT}/xml4rst CACHE INTERNAL "doxy (xml) output path")

    # xml outputs for breathe (or exhale). Conf might be different from
    # xml outputs for swig.
    set(DOXY_QUIET "YES")
    set(DOXY_WARNINGS "NO")
    set(GENERATE_HTML NO)
    set(GENERATE_XML YES)
    set(EXTRACT_ALL NO)
    set(EXTRACT_PRIVATE NO)
    set(DOXY_CONFIG_XML "${CMAKE_BINARY_DIR}/docs/config/${COMP}doxy-xml2rst.config")
    set(XML_OUTPUT xml4rst/${COMP})
    file(MAKE_DIRECTORY ${DOXYGEN_4_RST}/${COMP})
    # Set config file name
    set(DOXY_CONFIG_XML "${CMAKE_BINARY_DIR}/docs/config/${COMP}doxy-xml.config")

    # Add all subdirectories related to the current component into DOXYGEN_INPUTS
    set(XML_INPUTS)
    foreach(_dir ${${COMP}_DIRS})
      list(FIND ${COMP}_EXCLUDE_DOXY ${_dir} check_dir)
      if(NOT ${_dir} MATCHES test AND ${check_dir} EQUAL -1)
	list(APPEND XML_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
      endif()
    endforeach()
    if(XML_INPUTS)
      list(REMOVE_DUPLICATES XML_INPUTS)
    endif()
    set(DOXYGEN_INPUTS)
    foreach(_dir ${XML_INPUTS})
      set(DOXYGEN_INPUTS "${DOXYGEN_INPUTS} ${_dir}")
    endforeach()
    
    configure_file(${CMAKE_SOURCE_DIR}/docs/config/doxyxml2sphinx.config.in
      ${DOXY_CONFIG_XML} @ONLY)

    add_custom_target(${COMP}-doxy2xml
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG_XML}
      OUTPUT_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy4rst.log ERROR_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy4rst.log
      COMMENT "Generate xml/doxygen files for ${COMP} (conf: ${DOXY_CONFIG_XML}).") 

    
    set(SPHINX_API_DIR "${CMAKE_BINARY_DIR}/docs/sphinx/reference/")
    
    add_custom_target(${COMP}-doxy2rst
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
      "import buildtools; buildtools.create_breathe_files('${${COMP}_HDRS}', '${CMAKE_SOURCE_DIR}', '${COMP}', '${SPHINX_API_DIR}','${DOXY_CONFIG_XML}')"
      VERBATIM
      DEPENDS ${COMP}-doxy2xml
      )
    
    add_custom_command(TARGET  ${COMP}-doxy2rst POST_BUILD
      COMMENT "${COMP} : rst (c++ API) files have been generated in ${SPHINX_API_DIR}.")

endmacro()

# --------------------------------------
# Call this macro when configuration
# process is done for all components.
# 
# It will :
#  - create doxygen configs
# (for breathe, docstrings and so on)
#  - create all the targets related to
#   documentation
# --------------------------------------
macro(finalize_doc)
  if(WITH_DOCUMENTATION)

    # --- Target : html documentation from doxygen ---
    # run : make doxygen-html
    #  verbose mode : always off, since warnings
    # may be obtained with 'doxygen_warnings' target.
    set(DOXY_QUIET "YES")
    set(DOXY_WARNINGS "NO")
    set(GENERATE_HTML YES)
    set(GENERATE_XML NO)
    set(EXTRACT_ALL NO)
    if(USE_DEVEL_DOXYGEN) # OFF  by default. Activate to extract all.
      set(EXTRACT_ALL YES)
    endif()
    foreach(_dir ${DOXY_INPUTS})
      set(DOXYGEN_INPUTS "${DOXYGEN_INPUTS} ${_dir}")
    endforeach()
    configure_file(${CMAKE_SOURCE_DIR}/docs/config/doxy.config.in ${DOXY_CONFIG} @ONLY)
    add_custom_target(doxygen-html
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}
      COMMENT "Run doxygen ${DOXY_CONFIG} ...")
    
    add_custom_command(TARGET doxygen-html POST_BUILD
      COMMENT "Doxygen documentation has been built in : \n - ${DOXYGEN_OUTPUT} (xml) \n - ${DOC_ROOT_DIR}/html/doxygen (html).")
    
    # --- Target : create class diagrams from doxygen outputs, for sphinx. ---
    # run : make doxypng2sphinx
    # depends : doxygen-html
    add_custom_target(doxypng2sphinx
      COMMAND  ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/docs/find_doxygen_diagrams.py
      DEPENDS doxygen-html
      COMMENT "Browse doxygen outputs (graphs ...)")
    add_custom_command(TARGET doxypng2sphinx POST_BUILD
      COMMENT "Done. Generate class_diagrams.rst file in ${CMAKE_BINARY_DIR}/docs/sphinx/reference")

    add_dependencies(html doxypng2sphinx)

    # if(HAS_DOXYREST) # obsolete
    #   add_dependencies(doxyrest doxygen-html)
    #   add_dependencies(html doxyrest)
    # endif()

    ## TEST PURPOSE : maybe this tool is better than doxy2swig/exhale??
    # if(DOXY2RST)
    #   set(DOXYGEN_4_RST ${DOXYGEN_OUTPUT}/xml4rst)
    #   add_custom_target(doxy2rst
    # 	COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
    # 	"import buildtools; buildtools.parse_doxygen_wrapper('${DOXYGEN_4_RST}', '${CMAKE_BINARY_DIR}/docs/sphinx/reference')"
    # 	VERBATIM
    # 	)
    # endif()

    foreach(COMP ${COMPONENTS})
      message("ADD DEP TO ${COMP}")
      add_dependencies(html ${COMP}-doxy2rst)
    endforeach()

    # --- Generates conf.py, to describe sphinx setup ---
    # !! Should be call after doxygen setup
    # to have a correct DOXYGEN_INPUT value.
    configure_file (
      "${CMAKE_SOURCE_DIR}/docs/sphinx/conf.py.in"
      "${CMAKE_BINARY_DIR}/docs/sphinx/conf.py" @ONLY)

  endif()
  
  if(WITH_DOXYGEN_WARNINGS)
    add_custom_target(filter_warnings
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
      "import buildtools; buildtools.filter_doxygen_warnings_files('${CMAKE_BINARY_DIR}/doxygen_warnings', 'SUMMARY.warnings')"
      VERBATIM
      COMMENT "Filter doxygen warnings (result : SUMMARY.warnings)."
      )

  endif()


endmacro()
