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

# --------------------------------------
# Call this macro when all documentation
# process are over (doxygen, sphinx ...)
# to prepare doxygen config with
# an uptodate list of inputs.
# --------------------------------
macro(finalize_doc)
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
    # --- Generates conf.py, to describe sphinx setup ---
    # !! Should be call after doxygen setup
    # to have a correct DOXYGEN_INPUT value.
    configure_file (
      "${CMAKE_SOURCE_DIR}/docs/sphinx/conf.py.in"
      "${CMAKE_BINARY_DIR}/docs/sphinx/conf.py" @ONLY)
    configure_file(
      "${DOXY_CONFIG}"
      "${CMAKE_BINARY_DIR}/docs/sphinx/Doxyfile" @ONLY)
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
