include(SiconosTools)

#====================================================================
#
# Define and setup build process for a target for the current component
# 
# Usage:
#
# create_siconos_component(COMPONENT)
#
# The following variables must be properly set before any call to this function
# - <COMPONENT>_DIRS : list of directories (path relative to CMAKE_SOURCE_DIR)
#    that contain source files.
#
# This function:
#   creates a target <component> from all sources files in <component>_DIRS
#   excluding files listed after the keyword EXCLUDE.
#
function(create_siconos_component COMPONENT)

  set(multiValueArgs SRCDIRS  EXCLUDE)
  cmake_parse_arguments(component "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  # --- Collect source files from given directories ---

  # --> Scan source directories and return a list of files
  # to be compiled.
  get_sources(${COMPONENT} DIRS ${${COMPONENT}_DIRS} EXCLUDE ${component_EXCLUDE})
  
  # Create the library
  if(BUILD_SHARED_LIBS)
    add_library(${COMPONENT} SHARED ${${COMPONENT}_SRCS})
  else()
    add_library(${COMPONENT} STATIC ${${COMPONENT}_SRCS})
    set_property(TARGET ${COMPONENT} PROPERTY POSITION_INDEPENDENT_CODE ON)
  endif()

  # Set compiler options
  # reminder : WARNINGS_LEVEL=0 -> no warnings, =1, developers mininmal set of warnings,
  # =2 : strict mode, warnings to errors.
  apply_compiler_options(${COMPONENT} DIAGNOSTICS_LEVEL ${WARNINGS_LEVEL})
  
  # Append component source dirs to include directories
  # (Private : only to build current component).
  foreach(dir IN LISTS ${COMPONENT}_DIRS)
    target_include_directories(${COMPONENT} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}>)
  endforeach()

  # Add current component include dirs that should be propagated through
  # interface, in build tree.
  # WARNING : includes for install interface are handled later
  # in component install, and may be different.
  foreach(dir IN LISTS ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${COMPONENT} INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}>)
  endforeach()
    
  include(SiconosVersion)

  set_target_properties(${COMPONENT} PROPERTIES
    OUTPUT_NAME "siconos_${COMPONENT}"
    VERSION "${SICONOS_SOVERSION}"
    SOVERSION "${SICONOS_SOVERSION_MAJOR}")

  # windows stuff : this should be reviewed.
  include(WindowsLibrarySetup)
  windows_library_extra_setup(siconos_${COMPONENT} ${COMPONENT})
  
  if(WITH_GENERATION)
    # includes to be sent to python script for serialization/generation ...
    # Unstable and to be reviewed.
    foreach(_dir IN LISTS ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES)
      list(APPEND ${COMPONENT}_GENERATED_INCLUDES -I${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
    endforeach()
    set(GENERATED_INCLUDES "${GENERATED_INCLUDES};${${COMPONENT}_GENERATED_INCLUDES}"
      CACHE INTERNAL "")
  endif()
  
  
endfunction()


#====================================================================
#
# Define and setup documentation generation process
# for a the current component
# 
# Usage:
#
# configure_component_documentation(COMPONENT)
#
# The following variables must be properly set before any call to this function
# - <COMPONENT>_DIRS : list of directories (path relative to CMAKE_SOURCE_DIR)
#    that contain source files.
# - <COMPONENT>_EXCLUDE_DOXY : list of directories to exclude from doc
#   for this component.
#
#
function(configure_component_documentation COMPONENT)

  set(multiValueArgs HEADERS)
  cmake_parse_arguments(component "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  include(doc_tools)
  # --- doxygen warnings ---
  include(doxygen_warnings)
  
  # update the main doxy file, without building the doc
  if(WITH_${COMPONENT}_DOCUMENTATION  OR WITH_SERIALIZATION)
    # Prepare target to generate rst files from xml
    doxy2rst_sphinx(${COMPONENT} HEADERS ${component_HEADERS})
  endif()
endfunction()

#====================================================================
#
# Function to define and setup install process for a target for the current component
# 
# Usage:
#
# siconos_component_install_setup(<COMPONENT>)
#
#
#
function(siconos_component_install_setup COMPONENT)
  
  set(options NODOC)
  cmake_parse_arguments(component_install "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # libraries
  install(TARGETS ${COMPONENT}
    EXPORT siconosTargets
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION include/siconos
    )

  # Required for SiconosConfig.h
  target_include_directories(${COMPONENT} INTERFACE
    $<INSTALL_INTERFACE:include/siconos>)

  # Setup the list of all headers to be installed.
  foreach(dir IN LISTS ${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES)

    file(GLOB _headers CONFIGURE_DEPENDS
      LIST_DIRECTORIES false ${_FILE} ${dir}/*.h ${dir}/*.hpp)
    list(APPEND _all_headers ${_headers})
    
    # And each include path in install interface must obviously be installed ...
    # Note FP :  maybe we should have an explicit list of headers to be installed,
    # for each component, instead of a list of dirs?
  endforeach()

  if(_all_headers)
    # Do not install files listed in ${COMPONENT}_HDRS_EXCLUDE
    foreach(_file IN LISTS ${COMPONENT}_HDRS_EXCLUDE)
      list(REMOVE_ITEM _all_headers "${CMAKE_CURRENT_SOURCE_DIR}/${_file}")
    endforeach()
    # install files collected in _all_headers
    install(
      FILES ${_all_headers}
      DESTINATION include/siconos/${COMPONENT}
      )
    
    # Add include dirs in target interface 
    target_include_directories(${COMPONENT} INTERFACE
      $<INSTALL_INTERFACE:include/siconos/${COMPONENT}>)

  endif()
  # prepare documentation
  if(NOT component_install_NODOC)
    configure_component_documentation(${COMPONENT} HEADERS ${_all_headers})
  endif()
  # Set installed_targets list (for siconos-config.cmake file)
  list(APPEND installed_targets ${COMPONENT})
  list(REMOVE_DUPLICATES installed_targets)
  set(installed_targets ${installed_targets}
    CACHE INTERNAL "List of all exported components (targets).")
  
endfunction()


