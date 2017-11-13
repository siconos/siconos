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
  set(DOXY_CONFIG_XML "${CMAKE_BINARY_DIR}/Docs/config/${COMP}doxy.config.xml")
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
  configure_file(${CMAKE_SOURCE_DIR}/Docs/config/doxy.config.in ${DOXY_CONFIG_XML} @ONLY)
endmacro()


# --------------------------------
# Run doxygen to build documentation
# See CMAKE_BINARY_DIR/Docs/<component_name>doxy.log
# for errors and warnings.
# --------------------------------
macro(build_doc_xml COMP)
  set(confname ${CMAKE_BINARY_DIR}/Docs/config/${COMP}doxy.config.xml)
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
    foreach(_dir ${DOXY_INPUTS})
      set(DOXYGEN_INPUTS "${DOXYGEN_INPUTS} ${_dir}")
    endforeach()
    configure_file(${CMAKE_SOURCE_DIR}/Docs/config/doxy.config.in ${DOXY_CONFIG} @ONLY)
    configure_file(${CMAKE_SOURCE_DIR}/Docs/doxygen_layout/header.html.in
      Docs/doxygen_layout/header.html)
    add_custom_target(doxygen ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG})
    add_dependencies(html doxygen)

  endif()
endmacro()
