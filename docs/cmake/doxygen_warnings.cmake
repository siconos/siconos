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
# Default = ON.
# ===========================================================================

if(WITH_${COMPONENT}_DOXYGEN_WARNINGS)
  set(WITH_DOXYGEN_WARNINGS_INFILE "False" CACHE INTERNAL "Generate Doxygen warnings into a file.")
  foreach(_F IN LISTS ${COMPONENT}_SRCS)
    get_filename_component(_FP ${_F} PATH)
    get_filename_component(_FWE1 ${_F} NAME_WE)
    SET(_FWE ${_FP}/${_FWE1})
    IF(EXISTS ${_FWE}.hpp)
      SET(CURRENT_SICONOS_HEADER ${_FWE}.hpp)
    ELSE()
      IF(EXISTS ${_FWE}.h)
        SET(CURRENT_SICONOS_HEADER ${_FWE}.h)
      ENDIF()
    ENDIF()
    IF(WITH_DOXYGEN_WARNINGS_INFILE)
      SET(DOXYGEN_WARN_FILE ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings)
    ELSE()
      SET(DOXYGEN_WARN_FILE)
    ENDIF()
    CONFIGURE_FILE(${DOXY_WARNINGS_CONFIG} 
      ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.config)
    
    ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings
      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.config
      #COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings
      DEPENDS ${_F}
      DEPENDS ${CURRENT_SICONOS_HEADER}
      )
    SET_SOURCE_FILES_PROPERTIES(${_F} PROPERTIES OBJECT_DEPENDS ${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings)
    SET_SOURCE_FILES_PROPERTIES(${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.config PROPERTIES GENERATED TRUE)
    SET_SOURCE_FILES_PROPERTIES(${CMAKE_BINARY_DIR}/doxygen_warnings/${_FWE1}.warnings PROPERTIES GENERATED TRUE)
  endforeach()
endif()
