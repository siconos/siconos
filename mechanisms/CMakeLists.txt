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

include(occ_setup)

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

MESSAGE("################################# ${OpenCASCADE_INCLUDE_DIRS}" ${OpenCASCADE_INCLUDE_DIR})

target_link_libraries(${COMPONENT} PUBLIC $<BUILD_INTERFACE:OpenCASCADE::OpenCASCADE>)
target_compile_options(${COMPONENT} PRIVATE "-Wno-zero-as-null-pointer-constant")

# --- python bindings ---
if(WITH_${COMPONENT}_PYTHON_WRAPPER)
  add_subdirectory(swig)
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
