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
  src/core
  src/linalg
  src/optimization
  src/io
  src/utils
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
  src/FrictionContact/nsn_formulation
  src/FrictionContact/nsn_formulation/Generated
  src/RollingFrictionContact
  src/NCP
  src/MCP
  src/GenericMechanical
  src/VI
  src/Plasticity
  src/NonSmoothSolvers
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
  src/core
  src/linalg
  src/optimization
  src/io
  src/utils
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
  src/RollingFrictionContact
  src/NonSmoothSolvers
  src/NCP
  src/MCP
  src/GenericMechanical
  src/VI
  src/Plasticity	
  )

# ---- Final setup for the library ----

# -- Add parallel solvers
if(WITH_OPENMP AND WITH_PETSC)
  # Add parallel tools (graph coloring with PETSc)
  list(APPEND ${COMPONENT}_DIRS src/tools/parallel)
  list(APPEND ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES src/tools/parallel)

  # Add parallel LCP solvers
  list(APPEND ${COMPONENT}_DIRS src/LCP/parallel)

  # Add parallel FrictionContact solvers
  list(APPEND ${COMPONENT}_DIRS src/FrictionContact/parallel)

endif()

# -- create/setup component target --
include(ComponentSetup)


create_siconos_component(${COMPONENT})

# Windows stuff
include(WindowsNumericsSetup)

# --- Extra setup for the component ---
if(MLCPSIMPLEX_FOUND) 
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE CXX)
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX)
else()
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE C)
endif()

# Links with other Siconos components
target_link_libraries(numerics PRIVATE externals)

# Timers
# See also WITH_TIMERS option in advanced_options.cmake
include(CheckIncludeFile)
check_include_file("time.h"      HAVE_TIME_H)
check_include_file("sys/times.h" HAVE_SYS_TIMES_H)

# ---- Specific compilation options ---
set_source_files_properties(src//FrictionContact/FischerBurmeisterGenerated.c
  PROPERTIES COMPILE_FLAGS "-Wno-unused-value")

if(HAVE_QL0001) # From externals setup
  set_source_files_properties(src/LCP/lcp_nsqp.c src/LCP/lcp_qp.c PROPERTIES COMPILE_DEFINITIONS "HAVE_QL0001")  
endif()


if(WITH_MA57)
  target_sources(numerics PRIVATE src/linalg/optional/NM_MA57.c)
  list(APPEND ${COMPONENT}_HDRS_EXCLUDE src/linalg/optional/NM_MA57.h)
  target_include_directories(numerics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/linalg/optional>)
endif()

# - sort -
# --- Sort ---
if(HAVE_SORT)
    set_source_files_properties(src/linalg/NumericsSparseMatrix.c PROPERTIES COMPILE_DEFINITIONS "HAVE_SORT")  
    set_source_files_properties(src/tools/NumericsSparseMatrix.c
      PROPERTIES COMPILE_OPTIONS "-Wno-unused-function") # To remove all warnings due to include "sort.h"
endif()

# -- fc problems specific dump option
if(DUMP_PROBLEM)
  set_source_files_properties(
    src/FrictionContact/gfc2d_driver.c 
    src/FrictionContact/gfc3d_nonsmooth_Newton_AlartCurnier.c PROPERTIES COMPILE_DEFINITIONS "DUMP_PROBLEM")  
endif()

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
#    suite-sparse, gams

# - LpSolve -
find_package(LPSOLVE QUIET)
if(LPSOLVE_FOUND)
  # to add more tests in numerics
  set(HAS_ONE_LP_SOLVER TRUE CACHE INTERNAL "True if a LP Solver has been found and is activated.")
  # to add more tests in kernel and control
  target_link_libraries(numerics PRIVATE LPSOLVE::LPSOLVE)
  set_source_files_properties(src/utils/vertex_extraction.c PROPERTIES COMPILE_DEFINITIONS "HAS_ONE_LP_SOLVER") 
endif()

# - Thread -
# If Pthread is to be chosen in any case, set(CMAKE_THREAD_PREFER_PTHREAD 1)
find_package(Threads)
# https://cmake.org/cmake/help/v3.7/module/FindThreads.html?highlight=pthread
if(Threads_FOUND)
  target_link_libraries(numerics PRIVATE Threads::Threads)
endif()

include(FetchContent)

if(WITH_FCLIB)
  target_link_libraries(numerics PRIVATE fclib::fclib)
endif()

# - SuiteSparse -
if(WITH_PETSC)
  include(PetscSetup)
  target_link_libraries(numerics PRIVATE PkgConfig::PETSC)
#  target_link_libraries(numerics PRIVATE MPI::MPI_C)
  target_compile_definitions(numerics PUBLIC SICONOS_PETSC)
endif()


#   Note : the CSparse data structures are referred to in
#   kernel, but the functions are only called from numerics, so it is
#   not a link-time dependency for kernel.
find_package(SuiteSparse COMPONENTS CXSparse REQUIRED)
# We should turn public to private later (?)
# SuiteSparse headers are required by kernel
target_link_libraries(numerics PUBLIC $<BUILD_INTERFACE:SuiteSparse::CXSparse>)
target_include_directories(numerics PUBLIC $<INSTALL_INTERFACE:${SuiteSparse_CXSparse_INCLUDE_DIR}>)

find_package(SuiteSparse COMPONENTS LDL REQUIRED)
target_link_libraries(numerics PUBLIC $<BUILD_INTERFACE:SuiteSparse::LDL>)
target_include_directories(numerics PUBLIC $<INSTALL_INTERFACE:${SuiteSparse_LDL_INCLUDE_DIR}>)
# numerics possibly needs suite sparse installed on the system.
# -- Check size_t value --> choice of CSparse/CXSparse integer size
include(CheckTypeSize)
check_type_size(size_t SIZE_OF_SIZE_T LANGUAGE C)
# check_type_size("size_t" SIZE_OF_CSIX LANGUAGE CXX)
if("${SIZE_OF_SIZE_T}" EQUAL 8)
  # Use CS_LONG
  target_compile_definitions(numerics PUBLIC CS_LONG) # used to select CXSparse version (depends on int size), for CS_INT and so on
  target_compile_definitions(numerics PUBLIC LDL_LONG) # we use long interface of ldl (see ldl.h from suitesparse)
else()

endif()

target_compile_definitions(numerics PUBLIC NCOMPLEX) # We don't need the complex interface of CXSparse
#  set(NCOMPLEX 1 CACHE INTERNAL "") # for SiconosConfig.h Done with target_compile ...

# - MLCP Simplex -
find_package(MlcpSimplex QUIET)
if(MlcpSimplex_FOUND)
  set(HAVE_MLCPSIMPLEX TRUE CACHE INTERNAL "True if mlcp simplex solver API has been found and is activated.")
  target_link_libraries(numerics PRIVATE $<BUILD_INTERFACE:MlcpSimplex::MlcpSimplex>)
  set_source_files_properties(src/MLCP/mlcp_simplex.c PROPERTIES COMPILE_DEFINITIONS "HAVE_MLCPSIMPLEX")  
endif()

# - GAMS, Path -
# Checks GAMS, Path (Ferris, VI). See details in gams_path_setup.cmake file.
include(gams_path_setup)

# - UMFPack -
if(WITH_UMFPACK)
  find_package(UMFPACK REQUIRED)
  target_sources(numerics PRIVATE src/linalg/optional/NM_UMFPACK.c)
  target_link_libraries(numerics PUBLIC UMFPACK::UMFPACK)
  target_include_directories(numerics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/linalg/optional>)
endif()

# - SuperLU -

# Serial version
if(WITH_SUPERLU)
  find_package(SuperLU REQUIRED)
  target_sources(numerics PRIVATE src/linalg/optional/NM_SuperLU.c)
  target_include_directories(numerics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/linalg/optional>)
endif()

# Multithreaded version
if(WITH_SUPERLU_MT)
  find_package(SuperLU_MT REQUIRED COMPONENTS MT)
  target_sources(numerics PRIVATE src/linalg/optional/NM_SuperLU_MT.c)
  target_include_directories(numerics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/linalg/optional>)
endif()
# MPI/OpenMP version
if(WITH_SUPERLU_DIST)
  find_package(SuperLU_MT REQUIRED COMPONENTS DIST)
endif()
if(SuperLU_FOUND)
  
  target_link_libraries(numerics PRIVATE SuperLU::SuperLU)
endif()

# mumps and mpi

if(SICONOS_HAS_MPI)
  target_sources(numerics PRIVATE src/linalg/optional/NM_MPI.c)
  target_include_directories(numerics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/linalg/optional>)
  target_include_directories(numerics PUBLIC ${MPI_C_INCLUDE_DIRS})
  # Note FP : PUBLIC for the moment since it's required by tests.
  # This will be improved later.
  target_link_libraries(numerics PUBLIC ${MPI_C_LIBRARIES})
endif()

if(WITH_MUMPS)
  target_sources(numerics PRIVATE src/linalg/optional/NM_MUMPS.c)
  if(SICONOS_HAS_MPI)
    find_package(MUMPS REQUIRED COMPONENTS PAR)
  else()
    find_package(MUMPS REQUIRED COMPONENTS SEQ)
  endif()
  target_link_libraries(numerics PUBLIC MUMPS::MUMPS)
endif()

if(WITH_OPENSSL)
  find_package(OpenSSL REQUIRED COMPONENTS Crypto)
  target_link_libraries(numerics PRIVATE OpenSSL::Crypto)
  target_compile_definitions(numerics PUBLIC WITH_OPENSSL)  
endif()

if(WITH_FACTORISATION_DEBUG)
  set_source_files_properties(src/linalg/NumericsMatrix.c PROPERTIES COMPILE_DEFINITIONS "FACTORIZATION_DEBUG")  
endif()

# -- HDF5 --
if(WITH_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C HL)
  if(HDF5_IS_PARALLEL)
      message(STATUS "HDF5 was built with MPI support")
      find_package(MPI REQUIRED)
      target_link_libraries(numerics PRIVATE MPI::MPI_C HDF5::HDF5)
  else()
      target_link_libraries(numerics PRIVATE HDF5::HDF5)
  endif()
  target_compile_definitions(numerics PRIVATE WITH_HDF5)

endif()

# -- OpenMP --
if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(numerics  PRIVATE OpenMP::OpenMP_C)
  target_link_libraries(numerics  PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(numerics PUBLIC SICONOS_OMP)
endif()

# -- CUDA --
if(WITH_CUDA)

  # CUDA architectures (adjust based on target GPU)
  # 60=Pascal, 70=Volta, 75=Turing, 80=Ampere, 86=Ampere(RTX30xx), 89=Ada
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      message(FATAL_ERROR 
      " CMAKE_CUDA_ARCHITECTURES is not defined!\n"
      " Please specify the target GPU architecture(s) using:\n"
      "   -DCMAKE_CUDA_ARCHITECTURES=native (for local machine, if CMake version >= 3.24)\n"
      "   -DCMAKE_CUDA_ARCHITECTURES=80     (for Ampere/A100)\n"
      "   -DCMAKE_CUDA_ARCHITECTURES=\"80;90\" (for multiple)\n"
      " See: https://cmake.org/cmake/help/latest/prop_tgt/CUDA_ARCHITECTURES.html"
      )
  endif()
  
  enable_language(CUDA)


  # Required for cuSPARSE
  find_package(CUDAToolkit REQUIRED)

  # Link CUDA libraries
  target_link_libraries(${COMPONENT} PRIVATE
    CUDA::cudart
    CUDA::cusparse
    CUDA::cublas)

  target_compile_definitions(numerics PRIVATE SICONOS_CUDA)

  # Set CUDA standard
  set_target_properties(${COMPONENT} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CUDA_STANDARD ${CMAKE_CXX_STANDARD}
    CUDA_STANDARD_REQUIRED ON
    CUDA_SEPARABLE_COMPILATION ON)

  # Add CUDA files
  add_subdirectory(src/FrictionContact/parallel/gpu)
  target_sources(${COMPONENT} PRIVATE ${CUDA_FILES}) 

  message(STATUS "CUDA support enabled for numerics")
  message(STATUS "  CUDA Toolkit version: ${CUDAToolkit_VERSION}")
  message(STATUS "  CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
endif()


# --- tests ---
include(${COMPONENT}_tests.cmake)

# --- python bindings ---
if(WITH_PYB11_WRAPPER)
  add_subdirectory(python)
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
# set(${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES
#   src
#   src/tools
#   #   src/tools/internal
#   src/AVI
#   src/LCP
#   src/MCP
#   src/MLCP
#   src/Relay
#   src/FrictionContact
#   # src/FrictionContact/Generated
#   src/GenericMechanical
#   src/FrictionContact
#   )

set(${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES
${${COMPONENT}_DIRS} # All .hpp are installed
)

# List header files (in dir above) that must not be installed 
list(APPEND ${COMPONENT}_HDRS_EXCLUDE src/FrictionContact/Generated/funcodegen.h)

siconos_component_install_setup(${COMPONENT})
