#====================================================================
#
# Define and build a target for the current component
#
# Usage:
#
# library_project_setup(component_name)
#
#
# <COMPONENT>_SRCS : [optional] Project SRCS on per files basis
# <COMPONENT>_LIBS_NAME : [optional] libraries name. if it is empty the libs name are the same as COMPONENT.
# <COMPONENT>_DIRS : sources directories
# <COMPONENT>_Unstable_SRCS : built only if -DWITH_UNSTABLE=ON
# <COMPONENT>_VERSION : version of the library
# <COMPONENT>_HDRS : installation headers  (if none all headers)
# <COMPONENT>_HDRS_EXCLUDE_DIR : exclude headers from this dir from installation
# <COMPONENT>_LINKER_LANGUAGE : the language used to link the whole librarie

macro(LIBRARY_PROJECT_SETUP)

  # --- Collect source files from given directories ---
  # --> set ${COMPONENT}_SRCS
  get_sources("${${COMPONENT}_DIRS}")
  set(${COMPONENT}_SRCS ${${COMPONENT}_SRCS} ${SOURCES_FILES})
  # Unstable sources
  if(NOT WITH_${COMPONENT}_UNSTABLE)
    if(${COMPONENT}_Unstable_SRCS)
      foreach(_FILE ${${COMPONENT}_Unstable_SRCS})
        file(GLOB _GFILE ${_FILE})
        if(_GFILE)
          message("--  Source file excluded : ${_GFILE}")
          list(REMOVE_ITEM ${COMPONENT}_SRCS ${_GFILE})
        else()
          message("WARNING : Unstable file NOT FOUND : ${_FILE}")
	endif()
      endforeach()
    endif()
  endif()

  if(${COMPONENT}_EXCLUDE_SRCS)
    foreach(_FILE ${${COMPONENT}_EXCLUDE_SRCS})
      file(GLOB _GFILE ${_FILE})
      if(_GFILE)
        list(REMOVE_ITEM ${COMPONENT}_SRCS ${_GFILE})
      else()
        message("WARNING : file to be excluded NOT FOUND : ${_FILE}")
      endif()
    endforeach()
  endif()
  # --- Collect headers ---
  # --> set ${COMPONENT}_HDRS (for installation)
  set(HDRS_DIRS ${${COMPONENT}_DIRS})
  foreach(_DIR ${${COMPONENT}_HDRS_EXCLUDE_DIR})
    list(REMOVE_ITEM HDRS_DIRS ${_DIR})
  endforeach()

  get_headers("${HDRS_DIRS}")
  set(${COMPONENT}_HDRS "${HDRS_FILES}")

  # --- remove excluded headers if any ---
  # if not found, error to indicate that it should be removed.
  foreach(_HDR ${${COMPONENT}_HDRS_EXCLUDE})
    file(GLOB _EXCL_HDR ${_HDR})
    if (_EXCL_HDR)
      list(REMOVE_ITEM ${COMPONENT}_HDRS ${_EXCL_HDR})
    else()
      message(FATAL_ERROR "Excluded header ${_HDR} not found.")
    endif()
  endforeach()

  # --- include dirs for library ---
  # --> local includes (for build) : ${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
  # --> includes for external packages : SICONOS_INCLUDE_DIRECTORIES

  # --- Append local includes to ${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES ---
  remember_local_include_directories("${${COMPONENT}_DIRS}")

  set(${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
    ${${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES}
    ${CMAKE_CURRENT_BINARY_DIR}
    CACHE INTERNAL "Include directories for external dependencies.")

  include(doxygen_tools)
  # --- doxygen warnings ---
  include(doxygen_warnings)

  # --- doxygen documentation ---
  # xml files for python docstrings ...
  # xml files are required to build docstrings target
  # and so they must be build during cmake run.
  include(swig_python_tools)
  doxy2swig_docstrings(${COMPONENT})

  # update the main doxy file, without building the doc
  if(WITH_${COMPONENT}_DOCUMENTATION)
    update_doxy_config_file(${COMPONENT})
  endif()

  # -- include --
  # for local headers ...
  include_directories(${${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES})
  # and for headers of external libraries
  include_directories(${SICONOS_INCLUDE_DIRECTORIES})

  if(BUILD_SHARED_LIBS AND NOT BUILD_${COMPONENT}_STATIC)
    add_library(${COMPONENT} SHARED ${${COMPONENT}_SRCS})
  else()
    add_library(${COMPONENT} STATIC ${${COMPONENT}_SRCS})
    set_property(TARGET ${COMPONENT} PROPERTY POSITION_INDEPENDENT_CODE ON)
  endif()

  list(APPEND installed_targets ${COMPONENT})
  list(REMOVE_DUPLICATES installed_targets)
  set(installed_targets ${installed_targets}
    CACHE INTERNAL "Include directories for external dependencies.")
  set_target_properties(${COMPONENT} PROPERTIES 
    OUTPUT_NAME "${COMPONENT_LIBRARY_NAME}"
    VERSION "${SICONOS_SOVERSION}"
    SOVERSION "${SICONOS_SOVERSION_MAJOR}"
    CLEAN_DIRECT_OUTPUT 1 # no clobbering
    LINKER_LANGUAGE ${${COMPONENT}_LINKER_LANGUAGE})

  # windows stuff ...
  include(WindowsLibrarySetup)
  windows_library_extra_setup(${COMPONENT_LIBRARY_NAME} ${COMPONENT})
  # Link target with external libs ...
  target_link_libraries(${COMPONENT} ${PRIVATE} ${${COMPONENT}_LINK_LIBRARIES})
  
  if(BUILD_SHARED_LIBS)
    if(LINK_STATICALLY) # static linking is a nightmare
      set(REVERSE_LIST ${${COMPONENT}_LINK_LIBRARIES})
      LIST(REVERSE REVERSE_LIST)
      target_link_libraries(${COMPONENT} ${PRIVATE} ${REVERSE_LIST})
    endif()
  endif()    

  # ---- Installation ---
  # Headers
  install(FILES ${${COMPONENT}_HDRS} DESTINATION include/${PROJECT_NAME})
  
  # libraries
  install(TARGETS ${COMPONENT}
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION ${_install_lib}
    LIBRARY DESTINATION ${_install_lib})
  install(EXPORT ${PROJECT_NAME}Targets
    DESTINATION share/${PROJECT_NAME}/cmake)

  # --- python bindings ---
  if(WITH_${COMPONENT}_PYTHON_WRAPPER)
    add_subdirectory(swig)
  endif()

  
endmacro()

