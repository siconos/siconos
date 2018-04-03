# =======================================
# Macros and functions related to doxygen
# 
# =======================================

# --------------------------------
# Update doxy.config to take
# into account sources/headers files
# of COMP.
# --------------------------------
macro(update_doxy_config_file COMP)
  foreach(_dir ${${COMP}_DIRS})
    list(FIND ${COMP}_EXCLUDE_DOXY ${_dir} check_dir)
    if(NOT ${_dir} MATCHES test AND ${check_dir} EQUAL -1)
      list(APPEND DOXY_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
    endif()
  endforeach()
  list(REMOVE_DUPLICATES DOXY_INPUTS)
  set(DOXY_INPUTS ${DOXY_INPUTS} CACHE INTERNAL "doxy inputs")
endmacro()


# --------------------------------
# Update doxy.config to take
# into account sources/headers files
# of COMP.
# !! Generate only xml output !!
# Required for python-docstrings
# --------------------------------
macro(update_xml_doxy_config_file COMP)
  set(XML_INPUTS)
  set(DOXY_CONFIG_XML "${CMAKE_BINARY_DIR}/docs/config/${COMP}doxy.config.xml")
  foreach(_dir ${${COMP}_DIRS})
    list(FIND ${COMP}_EXCLUDE_DOXY ${_dir} check_dir)
    if(NOT ${_dir} MATCHES test AND ${check_dir} EQUAL -1)
      list(APPEND XML_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
    endif()
  endforeach()
  list(REMOVE_DUPLICATES XML_INPUTS)
  set(DOXYGEN_INPUTS)
  foreach(_dir ${XML_INPUTS})
    set(DOXYGEN_INPUTS "${DOXYGEN_INPUTS} ${_dir}")
  endforeach()
  set(GENERATE_HTML NO)
  set(GENERATE_XML YES)
  configure_file(${CMAKE_SOURCE_DIR}/docs/config/doxy.config.in ${DOXY_CONFIG_XML} @ONLY)
endmacro()


# --------------------------------
# Run doxygen to build documentation
# See CMAKE_BINARY_DIR/docs/<component_name>doxy.log
# for errors and warnings.
# --------------------------------
macro(build_doc_xml COMP)
  set(confname ${CMAKE_BINARY_DIR}/docs/config/${COMP}doxy.config.xml)
  execute_process(COMMAND ${DOXYGEN_EXECUTABLE} ${confname}
    OUTPUT_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy.log ERROR_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy.log)
  message(" -- Build xml doc for component ${COMP} ...")
endmacro()


# --------------------------------
# Call this macro after each
# component setup, to prepare
# doxygen config with a uptodate
# list of inputs.
# --------------------------------
macro(finalize_doxygen)
  if(WITH_DOCUMENTATION)
    # verbose mode.
    #  Always off, since warnings may be obtained with 'doxygen_warnings' target.
    set(DOXY_QUIET "YES")
    set(DOXY_WARNINGS "NO")
    set(GENERATE_HTML YES)
    set(GENERATE_XML YES)
    set(EXTRACT_ALL NO)
    if(USE_DEVEL_DOXYGEN)
      set(EXTRACT_ALL YES)
    endif()
    foreach(_dir ${DOXY_INPUTS})
      set(DOXYGEN_INPUTS "${DOXYGEN_INPUTS} ${_dir}")
    endforeach()
    configure_file(${CMAKE_SOURCE_DIR}/docs/config/doxy.config.in ${DOXY_CONFIG} @ONLY)
    add_custom_target(doxygen
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}
      COMMENT "Run doxygen ${DOXY_CONFIG} ...")
    add_custom_command(TARGET doxygen POST_BUILD
      COMMENT "Doxygen documentation has been built in : \n - ${DOXYGEN_OUTPUT} (xml) \n - ${DOC_ROOT_DIR}/html/doxygen (html).")
    
    add_custom_target(doxypng2sphinx
      COMMAND  ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/docs/find_doxygen_diagrams.py
      DEPENDS doxygen
      COMMENT "Browse doxygen outputs (graphs ...)")
    add_custom_command(TARGET doxypng2sphinx POST_BUILD
      COMMENT "Done. Generate class_diagrams.rst file in ${CMAKE_BINARY_DIR}/docs/sphinx/reference")
    add_dependencies(html doxypng2sphinx)
    add_dependencies(html doxygen)
    if(HAS_DOXYREST)
      add_dependencies(doxyrest doxygen)
      add_dependencies(html doxyrest)
    endif()
  endif()
endmacro()
