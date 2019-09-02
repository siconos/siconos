#================================================================
# cmake utilities to build numerics component
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

# Component name (i.e. target name)
set(COMPONENT numerics)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# --- Numerics optional dependencies ---
include(gams_setup)

# ------ source directories for current component ------
# What is needed by component to compile ?
# List here all directories that contain sources files
# for current component.
# Path must be relative to component path (i.e. to CMAKE_CURRENT_SOURCE_DIR)
set(${COMPONENT}_DIRS
  src/.
  src/tools
  src/tools/InterfaceToPathFerris
  src/tools/internal
  src/AVI
  src/LCP
  src/MLCP
  src/SOCP
  src/QP
  src/Relay
  src/FrictionContact
  src/FrictionContact/Generated
  src/NCP
  src/MCP
  src/GenericMechanical
  src/VI
)

# - Extras, optional -
# Static verifications (ACSL)
if(WITH_FRAMA_C)
  include(FramaC)
  file(GLOB frama_c_files src/FrictionContact/Generated/*.c)
  foreach(f ${frama_c_files})
    add_frama_c_test(${f} ENTRY main)
  endforeach()
  add_frama_c_test(${CMAKE_CURRENT_SOURCE_DIR}/src/tools/test/test_static_op3x3.c ENTRY main INCLUDES ${CMAKE_SOURCE_DIR}/externals/tools)
endif()

# -- Documentation --
# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
set(${COMPONENT}_EXCLUDE_DOXY src/Unstable_or_deprecated)

# ---- Final setup for the library ----

# Windows stuff
include(WindowsNumericsSetup)

### list(APPEND ${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
# List of specific headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE src/FrictionContact/Generated/funcodegen.h)
 

include(LibraryProjectSetup)
library_project_setup()


# --- Extra setup for the component ---
if(MLCPSIMPLEX_FOUND OR BUILD_AS_CPP) 
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE CXX)
else()
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE C)
endif()

# Links with other Siconos components
target_link_libraries(numerics PRIVATE externals)

# Links with non-Siconos libraries

# Blas/Lapack from externals
target_include_directories(numerics PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>/externals/blas_lapack)

if(CMAKE_DL_LIBS)
  # CMAKE_DL_LIBS needed for dlopen or equivalent
  target_link_libraries(numerics PRIVATE ${CMAKE_DL_LIBS})
endif()

# -- Optional deps --
# Fclib
if(WITH_FCLIB)
  find_package(FCLIB 3.0.0 CONFIG REQUIRED)
  if(FCLIB_NOT_FOUND)
    message(FATAL_ERROR "Can not find fclib package. Try to set FCLIB_ROOT.")
  endif()
  target_link_libraries(numerics PUBLIC FCLIB::fclib)
endif()

# SuiteSparse
# numerics possibly needs suite sparse installed on the system.
if(WITH_SYSTEM_SUITESPARSE)
  #   Note : the CSparse data structures are referred to in
  #   kernel, but the functions are only called from numerics, so it is
  #   not a link-time dependency for kernel.
  find_package(SuiteSparse COMPONENTS CXSparse)
  if (NOT SuiteSparse_FOUND OR NOT SuiteSparse_CXSparse_FOUND)
    message(FATAL_ERROR "\n!!! Cannot find suite sparse and/or cxsparse. Please install it or use internal (provided by Siconos) suite sparse by setting WITH_SYSTEM_SUITESPARSE to OFF !!!\n")
  endif()
  # We should turn public to private later (?)
  target_link_libraries(numerics PUBLIC SuiteSparse::CXSparse)
endif()

target_link_libraries(numerics INTERFACE $<BUILD_INTERFACE:LAPACK::LAPACK>)


if(BUILD_AS_CPP)
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX)
endif()

# mumps and mpi
if(WITH_MPI)
  target_include_directories(numerics PRIVATE ${MPI_C_INCLUDE_DIRS})
  # Note FP : PUBLIC for the moment since it's required by tests.
  # This will be improved later.
  target_link_libraries(numerics PRIVATE ${MPI_C_LIBRARIES})
  #   target_include_directories(numerics PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
  #   target_link_libraries(numerics PRIVATE ${MPI_Fortran_LIBRARIES})
endif()

if(WITH_MUMPS)
  find_package(MUMPS REQUIRED)
  target_link_libraries(numerics PRIVATE MUMPS::MUMPS)
endif()

# --- tests ---
include(${COMPONENT}_tests.cmake)

if(WITH_MPI AND WITH_${COMPONENT}_TESTING)
  # Note FP : temporary fix, to deal with PRIVATE deps of some components.
  # This will be reviewed later.
  target_include_directories(numerics-test PUBLIC ${MPI_C_INCLUDE_DIRS})
  target_link_libraries(numerics-test PUBLIC ${MPI_C_LIBRARIES})
endif()

