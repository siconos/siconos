# ===========================================================================
# Enable doxygen warnings (using a scan of headers of the current component)
#
# If WITH_COMPONENT_DOXYGEN_WARNINGS is ON,
# during build (i.e. make component), 
# generates a list of warnings produced by doxygen
# for each header of the current component
# 
# Usage:
#
# make doxygen_warnings to generate warnings for all headers
# 
# make doxygen_warnings_<component>_<path>_<header> for a single header
# e.g. make doxygen_warnings_kernel_modeling_LagrangianR
# 
# ===========================================================================


function(doxygen_warnings COMPONENT)

  set(multiValueArgs HEADERS)
  cmake_parse_arguments(component "" "" "${multiValueArgs}" ${ARGN})

  if(NOT WITH_${COMPONENT}_DOXYGEN_WARNINGS)
    return()
  endif()

  option(WITH_DOXYGEN_WARNINGS_INFILE "Generate Doxygen warnings into a file." ON)

  set(doxy_warnings_dir ${CMAKE_BINARY_DIR}/doxygen_warnings)

  include(doxycommon)
  set(DOXYGEN_PROJECT_NAME "Siconos - Doxygen Warnings")
  set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings)
  set(DOXYGEN_WARN_LAYOUT_FILE YES)
  set(DOXYGEN_WARN_AS_ERROR NO)
  set(DOXYGEN_WARN_FORMAT "$file:$line: $text")
  set(DOXYGEN_WARN_LINE_FORMAT "at line $line of file $file")
  set(DOXYGEN_GENERATE_HTML NO)
  set(DOXYGEN_GENERATE_XML YES)
  set(DOXYGEN_XML_OUTPUT xml)
  set(DOXYGEN_XML_PROGRAMLISTING YES)
  set(DOXYGEN_QUIET YES)
  set(DOXYGEN_WARNINGS YES)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
  set(DOXYGEN_WARN_IF_DOC_ERROR YES)
  set(DOXYGEN_WARN_IF_INCOMPLETE_DOC YES)
  set(DOXYGEN_WARN_NO_PARAMDOC YES)
  #set(DOXYGEN_EXTRACT_ALL NO)
  #if(USE_DEVEL_DOXYGEN)
    set(DOXYGEN_EXTRACT_ALL YES)
  #endif()
  set(DOXYGEN_EXTRACT_PRIVATE NO)
  if(NOT TARGET doxygen_warnings)
    add_custom_target(doxygen_warnings)
  endif()

  foreach(_F IN LISTS component_HEADERS)
    get_filename_component(_FWE1 ${_F} NAME_WE)
    get_filename_component(_FP ${_F} PATH)
    get_filename_component(_FP ${_FP} NAME_WE)
    if(WITH_DOXYGEN_WARNINGS_INFILE)
      set(DOXYGEN_WARN_LOGFILE ${doxy_warnings_dir}/${_FP}_${_FWE1}.warnings)
    else()
      unset(DOXYGEN_WARN_LOGFILE)
    endif()

    set(_target doxygen_warnings_${COMPONENT}_${_FP}_${_FWE1})
    doxygen_add_docs(${_target}
      ${_F}
      WORKING_DIRECTORY ${doxy_warnings_dir}
      USE_STAMP_FILE
      COMMENT "Checking Doxygen warnings for ${_F}"
      )

    if(TARGET ${COMPONENT})
      add_dependencies(doxygen_warnings ${_target})
    endif()
  endforeach()

endfunction()