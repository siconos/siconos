
# cmake utilities to build mechanics component
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

# This component has two optional sub-parts:
# - collision/bullet which requires bullet physics software  (and thus -DWITH_BULLET=ON cmake option)
# - occ which requires open-cascade community edition (and thus -DWITH_OpenCASCADE=ON cmake option)
# Those parts are automatically built and installed if bullet and/or opencascade are available and properly detected by cmake.


set(COMPONENT mechanics)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# ------ source directories for current component ------
# What is needed by component to compile ?
# List here all directories that contain sources files
# for current component.
# Path must be relative to component path (i.e. to CMAKE_CURRENT_SOURCE_DIR)

# First, 'native/minimal' mechanics part ---
set(${COMPONENT}_DIRS
  src/.
  src/collision
  src/collision/native
  src/collision/native/bodies
  src/joints
  )

# Then, optional parts ---

# -- Bullet --
# - add collision/bullet to the build process
if(WITH_BULLET)
  list(APPEND ${COMPONENT}_DIRS src/collision/bullet)
endif()

# -- OPenCascade --
# - add extra mechanics sources, depending on opencascade, to the build process
if(WITH_OpenCASCADE)
  list(APPEND ${COMPONENT}_DIRS src/occ)
endif()

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
target_link_libraries(mechanics PUBLIC
  externals # for boost/numerics/bindings
  numerics  # for debug.h ...
  kernel)

# Links with non-Siconos libraries
# dependency upon boost because of line 140 in SiconosPointers.hpp
# --> need boost::shared_array for HEM5 and lsodar.
# This should be reviewed (?)

target_include_directories(mechanics PRIVATE ${Boost_INCLUDE_DIRS})


# Bullet setup
# Different options are available. Check BulletSetup.cmake file
include(bullet_setup)

# -- Open Cascade --
if(WITH_OpenCASCADE)
  include(occ_setup)
  # Add opencascade headers and libs to components build
  target_link_libraries(${COMPONENT} PUBLIC $<BUILD_INTERFACE:OpenCASCADE::OpenCASCADE>)
  target_link_libraries(${COMPONENT} PUBLIC $<INSTALL_INTERFACE:OpenCASCADE::OpenCASCADE>)
endif()


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
