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
# - <COMPONENT>_EXCLUDE_SRCS : list of files to exclude from build process
#   for this component.
#
# This function:
#   creates a target <component> from all sources files in <component>_DIRS
#   excluding files from <component>_EXCLUDE_SRCS
#
function(create_siconos_component COMPONENT)

  set(multiValueArgs SRCDIRS  EXCLUDE)
  cmake_parse_arguments(component "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  # --- Collect source files from given directories ---

  # --> Scan source directories and return a list of files
  # to be compiled.
  get_sources(${COMPONENT} DIRS ${${COMPONENT}_DIRS} EXCLUDE ${${COMPONENT}_EXCLUDE_SRCS})
  
  # Create the library
  if(BUILD_SHARED_LIBS)
    add_library(${COMPONENT} SHARED ${${COMPONENT}_SRCS})
  else()
    add_library(${COMPONENT} STATIC ${${COMPONENT}_SRCS})
    set_property(TARGET ${COMPONENT} PROPERTY POSITION_INDEPENDENT_CODE ON)
  endif()

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
  
  configure_component_documentation(${COMPONENT})
  
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
# This function:
#   creates a target <component> from all sources files in <component>_DIRS
#   excluding files from <component>_EXCLUDE_SRCS
#
function(configure_component_documentation COMPONENT)
  
  include(doc_tools)
  # --- doxygen warnings ---
  include(doxygen_warnings)

  # --- documentation ---
  if(WITH_DOCUMENTATION OR WITH_DOXY2SWIG)
    # Update list of source directories to be taken
    # into account by doxygen for the current component
    # --> set CACHE var ${COMPONENT}_DOXYGEN_INPUTS
    # Required by doxy2swig_docstrings and doxy2rst_sphinx.
    update_doxygen_inputs(${COMPONENT})
  endif()
  
  # xml files for python docstrings ...
  # xml files are required to build docstrings target
  # and so they must be built during cmake run.
  if(WITH_PYTHON_WRAPPER)
    include(doxy2swig_docstrings)
    doxy2swig_docstrings(${COMPONENT})
  endif()
  
  # update the main doxy file, without building the doc
  if(WITH_${COMPONENT}_DOCUMENTATION)
    # Prepare target to generate rst files from xml
    doxy2rst_sphinx(${COMPONENT})
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
# This function 
#   creates a target <component> from all sources files in <component>_DIRS
#   excluding files from <component>_EXCLUDE_SRCS
#
function(siconos_component_install_setup COMPONENT)
  
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

    if(${CMAKE_VERSION} VERSION_GREATER "3.12.0")
      file(GLOB _headers CONFIGURE_DEPENDS
        LIST_DIRECTORIES false ${_FILE} ${dir}/*.h ${dir}/*.hpp)
    else()
      file(GLOB _headers
        LIST_DIRECTORIES false ${_FILE} ${_FILE} ${dir}/*.h ${dir}/*.hpp)
    endif()
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

  # Set installed_targets list (for siconos-config.cmake file)
  list(APPEND installed_targets ${COMPONENT})
  list(REMOVE_DUPLICATES installed_targets)
  set(installed_targets ${installed_targets}
    CACHE INTERNAL "List of all exported components (targets).")
  
endfunction()


