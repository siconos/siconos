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

# ------ include interface ------
# What is needed at build time
# by other targets to compile with current component.
# 
# It means that a call to
#  target_link_libraries(truc PRIVATE numerics)
# will imply -I<dirs> with dirs listed in
# ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES.
set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES
  src
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

# ---- Final setup for the library ----

# Windows stuff
include(WindowsNumericsSetup)

# -- create/setup component target --
include(ComponentSetup)
create_siconos_component(${COMPONENT})

# --- Extra setup for the component ---
if(MLCPSIMPLEX_FOUND OR BUILD_AS_CPP) 
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE CXX)
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX)
else()
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE C)
endif()

# Links with other Siconos components
target_link_libraries(numerics PRIVATE externals)

# Links with non-Siconos libraries

# Blas/Lapack from externals
target_include_directories(numerics PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>/externals/blas_lapack)
target_link_libraries(numerics PUBLIC $<BUILD_INTERFACE:LAPACK::LAPACK>)
#target_link_libraries(numerics INTERFACE LAPACK::LAPACK)

if(CMAKE_DL_LIBS)
  # CMAKE_DL_LIBS needed for dlopen or equivalent
  target_link_libraries(numerics PRIVATE ${CMAKE_DL_LIBS})
endif()

# -- Optional deps --

# Always searched and activated if found :
#    lpsolve, mlcpsimplex, pthread, path (Ferris and VI)
# Searched if asked in config (WITH_...=TRUE), error if not found :
#    fclib, suite-sparse, gams

# - LpSolve -
find_package(LPSOLVE QUIET)
if(LPSOLVE_FOUND)
  # to add more tests in numerics
  set(HAS_ONE_LP_SOLVER TRUE CACHE INTERNAL "True if a LP Solver has been found and is activated.")
  # to add more tests in kernel and control
  set(HAS_EXTREME_POINT_ALGO TRUE CACHE INTERNAL "True if lpsolve has been found and is activated.")
  target_link_libraries(numerics PRIVATE LPSOLVE::LPSOLVE)
endif()

# - Thread -
# If Pthread is to be chosen in any case, set(CMAKE_THREAD_PREFER_PTHREAD 1)
find_package(Threads)
# https://cmake.org/cmake/help/v3.7/module/FindThreads.html?highlight=pthread
if(Threads_FOUND)
  target_link_libraries(numerics PRIVATE Threads::Threads)
endif()

# - Fclib - 
if(WITH_FCLIB)
  find_package(FCLIB 3.0.0 CONFIG REQUIRED)
  if(FCLIB_NOT_FOUND)
    message(FATAL_ERROR "Can not find fclib package. Try to set FCLIB_ROOT.")
  endif()
  target_link_libraries(numerics PUBLIC FCLIB::fclib)
endif()

# - SuiteSparse -
# numerics possibly needs suite sparse installed on the system.
# -- Check size_t value --> choice of CSparse/CXSparse integer size
if(NOT DEFINED SICONOS_INT64)
  if(NOT SIZE_OF_CSI)
    include(CheckTypeSize)
    check_type_size("size_t" SIZE_OF_CSI)
    if(NOT SIZE_OF_CSI) # This won't happen with c11
      message(FATAL_ERROR
        "Could not get size of size_t, please specify SICONOS_INT64.")
    endif()
  endif()

  if("${SIZE_OF_CSI}" EQUAL 8)
    set(SICONOS_INT64 TRUE CACHE INTERNAL "")
  else()
    set(SICONOS_INT64 FALSE CACHE INTERNAL "")
  endif()
endif()

if(WITH_SYSTEM_SUITESPARSE)
  #   Note : the CSparse data structures are referred to in
  #   kernel, but the functions are only called from numerics, so it is
  #   not a link-time dependency for kernel.
  find_package(SuiteSparse COMPONENTS CXSparse REQUIRED)
  # We should turn public to private later (?)
  # SuiteSparse headers are required by kernel
  target_link_libraries(numerics PUBLIC $<BUILD_INTERFACE:SuiteSparse::CXSparse>)
  target_include_directories(numerics PUBLIC $<INSTALL_INTERFACE:${SuiteSparse_CXSparse_INCLUDE_DIR}>)
endif()

# - MLCP Simplex -
find_package(MlcpSimplex QUIET)
if(MlcpSimplex_FOUND)
  set(HAVE_MLCPSIMPLEX TRUE CACHE INTERNAL "True if mlcp simplex solver API has been found and is activated.")
  target_link_libraries(numerics PRIVATE $<BUILD_INTERFACE:MlcpSimplex::MlcpSimplex>)
endif()

# - GAMS, Path -
# Checks GAMS, Path (Ferris, VI). See details in gams_path_setup.cmake file.
include(gams_path_setup)

# - UMFPack -
if(WITH_UMFPACK)
  find_package(UMFPACK REQUIRED)
  target_link_libraries(numerics PUBLIC UMFPACK::UMFPACK)
endif()

# - SuperLU -

# Serial version
if(WITH_SUPERLU)
  find_package(SuperLU REQUIRED)
endif()
# Multithreaded version
if(WITH_SUPERLU_MT)
  find_package(SuperLU_MT REQUIRED COMPONENTS MT)
endif()
# MPI/OpenMP version
if(WITH_SUPERLU_DIST)
  find_package(SuperLU_MT REQUIRED COMPONENTS DIST)
endif()
if(SuperLU_FOUND)
  
  target_link_libraries(numerics PRIVATE SuperLU::SuperLU)
endif()

# mumps and mpi
if(WITH_MPI)
  target_include_directories(numerics PUBLIC ${MPI_C_INCLUDE_DIRS})
  # Note FP : PUBLIC for the moment since it's required by tests.
  # This will be improved later.
  target_link_libraries(numerics PUBLIC ${MPI_C_LIBRARIES})
  #   target_include_directories(numerics PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
  #   target_link_libraries(numerics PRIVATE ${MPI_Fortran_LIBRARIES})
endif()

if(WITH_MUMPS)
  if(WITH_MPI)
    find_package(MUMPS REQUIRED COMPONENTS PAR)
  else()
    find_package(MUMPS REQUIRED COMPONENTS SEQ)
  endif()
  target_link_libraries(numerics PRIVATE MUMPS::MUMPS)
endif()

# -- HDF5 --
# For logging
if(WITH_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C HL)
  target_include_directories(numerics PRIVATE ${HDF5_C_INCLUDE_DIRS})
  target_link_libraries(numerics PRIVATE ${HDF5_C_LIBRARIES})
  target_link_libraries(numerics PRIVATE ${HDF5_HL_LIBRARIES})
endif()

# -- OpenMP --
if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  # Notice : cmake  >= 3.9.6 is advised for a proper handling of openmp
  if(${CMAKE_VERSION} VERSION_GREATER "3.9.6")
    target_link_libraries(numerics  PRIVATE OpenMP::OpenMP_C)
    target_link_libraries(numerics  PRIVATE OpenMP::OpenMP_CXX)
  else()
    message("Openmp flags : ${OpenMP_C_FLAGS}")
    target_compile_options(numerics PRIVATE  ${OpenMP_C_FLAGS})
  endif()
endif()

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
  src
  src/tools
  #   src/tools/internal
  src/AVI
  src/LCP
  src/MCP
  src/MLCP
  src/Relay
  src/FrictionContact
  # src/FrictionContact/Generated
  # src/GenericMechanical
  )

# List header files (in dir above) that must not be installed 
set(${COMPONENT}_HDRS_EXCLUDE src/FrictionContact/Generated/funcodegen.h)

siconos_component_install_setup(${COMPONENT})

# --- tests ---
include(${COMPONENT}_tests.cmake)

if(WITH_MPI AND WITH_${COMPONENT}_TESTING)
  # Note FP : temporary fix, to deal with PRIVATE deps of some components.
  # This will be reviewed later.
  target_include_directories(numerics-test PUBLIC ${MPI_C_INCLUDE_DIRS})
  target_link_libraries(numerics-test PUBLIC ${MPI_C_LIBRARIES})
endif()

