#================================================================
# cmake utilities to build mechanisms component
#================================================================
#
# The objective is to call component_setup to create the target <COMPONENT>.
# Before, it's necessary to set:
# 
# - COMPONENT component name
# - <COMPONENT>_DIRS: the list of paths (relative to CMAKE_CURRENT_SOURCE_DIR) that
#   contain source files
# - <COMPONENT>_EXCLUDE_SRCS: the list of files, in <COMPONENT>_DIRS, that must be excluded
#   from build.
# - <COMPONENT>_INTERFACE_INCLUDE_DIRECTORIES: a list of directories
#   to populate the interface of the target for include directories at build time

set(COMPONENT mechanisms)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# ------ source directories for current component ------
# What is needed by component to compile ?
# List here all directories that contain sources files
# for current component.
# Path must be relative to component path (i.e. to CMAKE_CURRENT_SOURCE_DIR)

set(${COMPONENT}_DIRS
  src/.
  src/CADMBTB
  src/MBTB
  src/TIMERMBTB)


# -- Documentation --
# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
# set(${COMPONENT}_EXCLUDE_DOXY)

# ------ include interface ------
# What is needed at build time
# by other targets to compile with current component.
# 
# It means that a call to
#  target_link_libraries(truc PRIVATE mechanics)
# will imply -I<dirs> with dirs listed in
# ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES.
set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES ${${COMPONENT}_DIRS})
  
# -- create/setup component target --
include(ComponentSetup)
create_siconos_component(${COMPONENT})

# --- Extra setup for the component ---

# Links with other Siconos components
target_link_libraries(mechanisms PUBLIC
  externals
  numerics 
  kernel
  mechanics)

# Links with non-Siconos libraries
# --- Search component dependencies ---
if(NOT WITH_OCE)
  message(FATAL_ERROR "Please run cmake again with option WITH_OCE=ON (required for mechanisms) or do not setup Siconos with mechanisms component.")
endif()
include(oce_setup)
# Add oce headers to search path
# 'private' add limited to current target would be better.
# But those headers are missing during swig process. Todo : find a proper
# way to deal with this.
target_include_directories(${COMPONENT} PUBLIC ${OCE_INCLUDE_DIRS})
# Add oce libs to components links
target_link_libraries(${COMPONENT} PUBLIC ${OCE_LIBRARIES})

# --- python bindings ---
if(WITH_${COMPONENT}_PYTHON_WRAPPER)
  add_subdirectory(swig)
endif()

if(WITH_PYTHON_WRAPPER)
  add_dependencies(${COMPONENT} ${COMPONENT}_docstrings)
endif()  

# ---- Installation ----
# Call siconos_component_install_setup(<COMPONENT>)
# to prepare installation of the current target.
#
# Before, it's necessary to set:
# 
# - <COMPONENT>_INSTALL_INTERFACE_INCLUDE_DIRECTORIES with all directories
#    that contain headers files that must be installed.
# 
set(${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES
  ${${COMPONENT}_DIRS} # All .hpp are installed
  )

siconos_component_install_setup(${COMPONENT})

# --- tests ---
include(${COMPONENT}_tests.cmake)
