# ===========================================================================
# Enable doxygen warnings (using a scan of headers of the current component)
#
# If WITH_COMPONENT_DOXYGEN_WARNINGS is ON,
# during build (i.e. make component), 
# generates a list of warnings produced by doxygen
# for each header of the current component
#
# A config file is generated for each header, from ${DOXY_WARNINGS_CONFIG} file
# set docs/CMakeLists.txt
#
# Use -DWITH_DOXYGEN_WARNINGS_INFILE=ON to save outputs in files.
# Default = OFF.
# ===========================================================================

if(WITH_${COMPONENT}_DOXYGEN_WARNINGS)
  set(WITH_DOXYGEN_WARNINGS_INFILE "False" CACHE INTERNAL "Generate Doxygen warnings into a file.")
  foreach(_F IN LISTS ${COMPONENT}_SRCS)
    get_filename_component(_FP ${_F} PATH)
    get_filename_component(_FWE1 ${_F} NAME_WE)
    set(_FWE ${_FP}/${_FWE1})
    if(EXISTS ${_FWE}.hpp)
      set(CURRENT_SICONOS_HEADER ${_FWE}.hpp)
    else()
      if(EXISTS ${_FWE}.h)
        set(CURRENT_SICONOS_HEADER ${_FWE}.h)
      endif()
    endif()
    if(WITH_DOXYGEN_WARNINGS_INFILE)
      set(DOXYGEN_WARN_FILE ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings)
    else()
      set(DOXYGEN_WARN_FILE)
    endif()
    #   MUST BE set TO set(OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings)
    set(DOXY_QUIET "YES")
    set(DOXY_WARNINGS "YES")
    set(DOXYGEN_INPUTS ${CURRENT_SICONOS_HEADER})
    set(GENERATE_HTML NO)
    set(GENERATE_XML YES)
    set(XML_OUTPUT xml)
    set(EXTRACT_ALL NO)
    if(USE_DEVEL_DOXYGEN) # OFF  by default. Activate to extract all.
      set(EXTRACT_ALL YES)
    endif()
    set(EXTRACT_PRIVATE NO)

    configure_file(${DOXY_WARNINGS_CONFIG} 
      ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.config)

    add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings
      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.config
      #COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings
      DEPENDS ${_F}
      DEPENDS ${CURRENT_SICONOS_HEADER}
      )
    set_source_files_properties(${_F} PROPERTIES OBJECT_DEPENDS ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings)
    set_source_files_properties(${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.config PROPERTIES GENERATED TRUE)
    set_source_files_properties(${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings PROPERTIES GENERATED TRUE)
  endforeach()
endif()
