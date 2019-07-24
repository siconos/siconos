#================================================================
#
# Look for siconos dependencies, common to several components.
#
#================================================================

# --- Blas Lapack ---
# include(BlasLapackSetup)
compile_with(BLAS REQUIRED SICONOS_COMPONENTS kernel numerics externals)
compile_with(LAPACK REQUIRED SICONOS_COMPONENTS kernel numerics externals)
if(NOT BLAS_INCLUDE_DIRS)
  message(FATAL_ERROR "cannot find blas include directories")
endif()
if(NOT LAPACK_INCLUDE_DIRS)
  message(FATAL_ERROR "cannot find lapack include directories")
endif()

# --- Numerics optional dependencies ---
compile_with(MlcpSimplex SICONOS_COMPONENTS numerics)
compile_with(Pthread SICONOS_COMPONENTS numerics)
IF(GAMS_DIR)
  SET(GAMS_C_API_FIND_REQUIRED TRUE)
ENDIF(GAMS_DIR)
COMPILE_WITH(GamsCApi SICONOS_COMPONENTS numerics)
IF(GAMSCAPI_FOUND)
  # needed for siconosconfig.h
  IF(NOT GAMS_DIR)
    SET(GAMS_DIR " ")
  ENDIF(NOT GAMS_DIR)
  SET(GAMS_MODELS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/numerics/share/gams")
  SET(GAMS_MODELS_SHARE_DIR "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/gams")
  #  IF(UNIX AND NOT APPLE)
  #    SET(SICONOS_DEFAULT_LINKER_OPTIONS "${SICONOS_DEFAULT_LINKER_OPTIONS} -Wl,-rpath-link,${GAMS_DIR}")
  #  ENDIF(UNIX AND NOT APPLE)
ENDIF(GAMSCAPI_FOUND)

# --- SuiteSparse ---
# Look for system-installed SuiteSparse/CSparse
if (WITH_SYSTEM_SUITESPARSE)
  compile_with(SuiteSparse COMPONENTS CXSparse
    SICONOS_COMPONENTS externals numerics)
  # Note on the above: The CSparse data structures are referred to in
  # kernel, but the functions are only called from numerics, so it is
  # not a link-time dependency for kernel.
  if (NOT SuiteSparse_FOUND OR NOT SuiteSparse_CXSparse_FOUND)
    set(_sys_CXSparse FALSE)
    message(STATUS "System SuiteSparse was requested (WITH_SYSTEM_SUITESPARSE=${WITH_SYSTEM_SUITESPARSE})\ 
    but not found! Using the internal copy of suitesparse")
  else()
    set(_sys_CXSparse TRUE)
  endif()
  set(USE_SYSTEM_SUITESPARSE ${_sys_CXSparse} CACHE INTERNAL "flag to check to systemwide SuiteSparse install")
endif()

# --- Other solvers ---
compile_with(PathFerris SICONOS_COMPONENTS numerics)
compile_with(PathVI SICONOS_COMPONENTS numerics)
compile_with(LpSolve SICONOS_COMPONENTS numerics)
if(LpSolve_FOUND)
  set(HAS_ONE_LP_SOLVER TRUE)
  set(HAS_EXTREME_POINT_ALGO TRUE)
  set(WITH_LPSOLVE TRUE)
endif(LpSolve_FOUND)

# --- Sort ---
if(EXISTS "${CMAKE_SOURCE_DIR}/externals/sort/sort.h")
  if(EXISTS "${CMAKE_SOURCE_DIR}/externals/sort/sort_common.h")
    set(HAVE_SORT TRUE)
  endif()
endif()

# --- ql0001 ---
if(EXISTS "${CMAKE_SOURCE_DIR}/externals/optim_misc/ql0001/ql0001.f")
  set(HAVE_QL0001 TRUE)
endif()

# --- Mumps ---
# if(WITH_MUMPS)
#   if(MPI_FOUND)
#     # Fedora allow parallel install of MPI and vanilla version of MUMPS.
#     # This shouldn't hurt in any case ...
#     if(NOT MUMPS_LIBRARY_DIRECTORY)
#       get_filename_component(MUMPS_LIBRARY_DIRECTORY "${MPI_LIBRARY}" PATH)
#     endif()
#   endif()
#   compile_with(MUMPS REQUIRED SICONOS_COMPONENTS numerics)
# endif()

# --- UMFPACK ---
if(WITH_UMFPACK)
  compile_with(Umfpack REQUIRED SICONOS_COMPONENTS numerics)
endif()

# --- SUPERLU ---
IF (WITH_SUPERLU AND WITH_SUPERLU_MT)
  message(FATAL_ERROR "Both SuperLU and SuperLU_MT are enabled. Due to symbol collision, both cannot be enabled at the same time")
ENDIF()

if(WITH_SUPERLU)
  compile_with(SuperLU REQUIRED SICONOS_COMPONENTS numerics)
endif()

if(WITH_SUPERLU_MT)
  compile_with(SuperLU_MT REQUIRED SICONOS_COMPONENTS numerics)
endif()

# not ready yet
#if(WITH_SUPERLU_dist)
#  compile_with(SuperLU_dist REQUIRED SICONOS_COMPONENTS numerics)
#endif()

# --- Fclib ---
IF(WITH_FCLIB)
  COMPILE_WITH(FCLIB REQUIRED SICONOS_COMPONENTS numerics)
  IF(FCLib_FCLIB_HEADER_ONLY)
    COMPILE_WITH(HDF5 REQUIRED COMPONENTS C HL SICONOS_COMPONENTS numerics )
  ELSE()
    APPEND_C_FLAGS("-DFCLIB_NOT_HEADER_ONLY")
ENDIF()
  IF(FCLIB_NOTFOUND)
    # try the package stuff
    # need FCLib_DIR !!
    COMPILE_WITH(FCLib 1.0 REQUIRED SICONOS_COMPONENTS numerics)
  ENDIF()
ENDIF()

# GMP
compile_with(GMP REQUIRED SICONOS_COMPONENTS kernel)

IF(WITH_CXX)
  # --- Boost ---
  compile_with(Boost 1.47 REQUIRED)
ENDIF(WITH_CXX)

#SET(WITH_BOOST_LOG TRUE)
IF(WITH_BOOST_LOG)
  compile_with(Boost 1.47 COMPONENTS log  REQUIRED)
  APPEND_CXX_FLAGS("-DBOOST_LOG_DYN_LINK")
ENDIF()

# -- VTK --
IF(WITH_VTK)
  COMPILE_WITH(VTK SICONOS_COMPONENTS mechanics)
  IF(VTK_FOUND)
    MESSAGE(STATUS "Found vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
    SET(SICONOS_HAVE_VTK TRUE)
    IF(VTK_USE_FILE)
      INCLUDE(${VTK_USE_FILE})
    ENDIF()
  ENDIF()
ENDIF()

# -- FreeCAD --
# For python bindings
if(WITH_FREECAD)
  compile_with(FreeCAD COMPONENTS Part REQUIRED SICONOS_COMPONENTS mechanics)
endif()


# -- HDF5 --
# For logging in Numerics
IF(WITH_HDF5)
  COMPILE_WITH(HDF5 REQUIRED COMPONENTS C HL SICONOS_COMPONENTS numerics)
ENDIF(WITH_HDF5)

#
# -- Serialization --
#
include(serialization_vector_test)
if(WITH_SERIALIZATION)
  COMPILE_WITH(Boost 1.47
    COMPONENTS serialization filesystem REQUIRED)
  if (Boost_VERSION GREATER 106100)
    # If boost is recent enough, prefer system boost serialization to
    # the one included in "externals/boost_serialization".
    set(WITH_SYSTEM_BOOST_SERIALIZATION ON)
  endif()
  TEST_SERIALIZATION_VECTOR_BUG()
endif()

# -- Python bindings --
if(WITH_PYTHON_WRAPPER)
  compile_with(Numpy REQUIRED)
  # trick (required with nix stuff) to force swig3 rather than swig2, if available.
  find_program(SWIG_EXECUTABLE NAMES swig swig3.0 swig2.0 PATHS ENV PATH)
  compile_with(SWIG 2.0.11 REQUIRED)
  include(${SWIG_USE_FILE})
endif()

# See if help2man is available
find_program(HELP2MAN help2man)

#
# Fedora13 https://fedoraproject.org/wiki/UnderstandingDSOLinkChange
if(UNIX)
  # add -lm to linker
  CHECK_C_COMPILER_FLAG("-lm" C_HAVE_LINKER_M)
  if(C_HAVE_LINKER_M)
    set(SICONOS_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES} "m" CACHE INTERNAL "List of external libraries")
  endif()
endif()

# man pages
IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/man)
  CONFIGURE_FILE(man/siconos.1.in man/siconos.1)
  INSTALL(FILES ${CMAKE_BINARY_DIR}/man/siconos.1 DESTINATION share/man/man1)
ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/man)

set(${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
  ${${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES}
  ${CMAKE_BINARY_DIR}
  CACHE INTERNAL "Include directories for external dependencies.")
