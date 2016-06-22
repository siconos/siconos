#================================================================
#
# Look for siconos dependencies, common to several components.
#
#================================================================

# --- Blas Lapack ---
# include(BlasLapackSetup)
compile_with(BLAS REQUIRED)
compile_with(LAPACK REQUIRED)
if(NOT BLAS_INCLUDE_DIRS)
  message(FATAL_ERROR "cannot find blas include directories")
endif()
if(NOT LAPACK_INCLUDE_DIRS)
  message(FATAL_ERROR "cannot find lapack include directories")
endif()

# --- Numerics optional dependencies ---
compile_with(MlcpSimplex)
compile_with(Pthread)
IF(GAMS_DIR)
  SET(GAMS_C_API_FIND_REQUIRED TRUE)
  COMPILE_WITH(GamsCApi)
  # needed for siconosconfig.h
  SET(GAMS_MODELS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/numerics/share/gams")
  SET(GAMS_MODELS_SHARE_DIR "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/gams")
ENDIF(GAMS_DIR)
compile_with(PathFerris)
compile_with(LpSolve)
if(LpSolve_FOUND)
  set(HAS_ONE_LP_SOLVER TRUE)
  set(HAS_EXTREME_POINT_ALGO TRUE)
  set(WITH_LPSOLVE TRUE)
  if(WITH_CXX)
    string(REPLACE "-Werror=conversion" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  endif()
  string(REPLACE "-Werror=conversion" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
endif(LpSolve_FOUND)

# --- Mumps ---
if(WITH_MUMPS)
  if(NOT IDONTWANTMPI)
    compile_with(MPI REQUIRED)
  endif(NOT IDONTWANTMPI)
  if(MPI_FOUND)
    set(HAVE_MPI TRUE)
    # Fedora allow parallel install of MPI and vanilla version of MUMPS.
    # This shouldn't hurt in any case ...
    if(NOT MUMPS_LIBRARY_DIRECTORY)
      get_filename_component(MUMPS_LIBRARY_DIRECTORY "${MPI_LIBRARY}" PATH)
    endif()
  endif()
  compile_with(MUMPS REQUIRED)
endif()

# --- UMFPACK ---
if(WITH_UMFPACK)
  compile_with(Umfpack REQUIRED)
endif()

# --- Fclib ---
IF(WITH_FCLIB)
  COMPILE_WITH(FCLIB REQUIRED)   
  IF(FCLIB_NOTFOUND)
    # try the package stuff
    # need FCLib_DIR !!
    COMPILE_WITH(FCLib 1.0 REQUIRED)
  ENDIF()
ENDIF()

# --- Boost ---
compile_with(Boost 1.47 REQUIRED)

#SET(WITH_BOOST_LOG TRUE)
IF(WITH_BOOST_LOG)
  compile_with(Boost 1.47 COMPONENTS log  REQUIRED)
  APPEND_CXX_FLAGS("-DBOOST_LOG_DYN_LINK")
ENDIF()
# --- Bullet ---
IF(WITH_BULLET)
  COMPILE_WITH(Bullet REQUIRED)
  IF(BULLET_FOUND)
    SET(HAVE_BULLET TRUE)
    MESSAGE( STATUS " Bullet include dirs : ${BULLET_INCLUDE_DIRS}" )
    IF(BULLET_USE_DOUBLE_PRECISION)
      APPEND_CXX_FLAGS("-DBT_USE_DOUBLE_PRECISION")
    ENDIF(BULLET_USE_DOUBLE_PRECISION)
  ENDIF(BULLET_FOUND)
ENDIF(WITH_BULLET)

# --- OCC ---
IF(WITH_OCC)
  if(NOT WITH_MECHANISMS)
    COMPILE_WITH(OCE 0.15 REQUIRED ONLY mechanics)
    SET(HAVE_OCC TRUE)
    LIST(REMOVE_ITEM SICONOS_LINK_LIBRARIES DRAWEXE)
  endif()
ENDIF()

# --- Mechanisms (Saladyn?) ---
IF(WITH_MECHANISMS)
  SET(OCE_TOOLKITS "TKernel"  "TKMath" "TKService" "TKV3d"  "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")

  message(STATUS "Searching for OCE ....")
  compile_with(OCE COMPONENTS ${OCE_TOOLKITS} ONLY mechanics)
  if(OCE_VERSION VERSION_LESS 0.16)
    # DRAWEXE link fails on some systems and must be removed.
    list(REMOVE_ITEM mechanics_LINK_LIBRARIES DRAWEXE)
    list(REMOVE_ITEM SICONOS_LINK_LIBRARIES DRAWEXE)
    set(SICONOS_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES} "m" CACHE INTERNAL "List of external libraries")
    set(mechanics_LINK_LIBRARIES ${mechanics_LINK_LIBRARIES} "m" CACHE INTERNAL "List of external libraries")
  endif()
  if(OCE_FOUND)
    message(STATUS "Found OCE version ${OCE_VERSION}")
    if(NOT OCE_ALL_FOUND)
      set(OCE_FOUND false)
      message(WARNING "Ignoring OCE installation due to missing toolkit(s): ${OCE_MISSING_TOOLKITS}")
    endif(NOT OCE_ALL_FOUND)
  endif(OCE_FOUND)

  if(OCE_FOUND)
    # Include files reside in ${OCE_INCLUDE_DIRS};
    #    include_directories(${OCE_INCLUDE_DIRS})
    # We do not need library path, they will be automatically imported.
  else(OCE_FOUND)
    # OCE not found; either it is not found and user
    # has to set OCE_DIR to the directory containing
    # OCEConfig.cmake, or OCE is not installed and we
    # try to find OpenCascade files.
    message(STATUS "OCE not found.  Try to find OpenCascade files.")

    FIND_PACKAGE(OpenCASCADE REQUIRED COMPONENTS ${OCE_TOOLKITS})
    COMPILE_WITH(OpenCASCADE)

    IF(OpenCASCADE_FOUND)
      message(STATUS "OpenCASCADE_INCLUDE_DIR = " ${OpenCASCADE_INCLUDE_DIR})
      message(STATUS "OpenCASCADE_LIBRARIES = " ${OpenCASCADE_LIBRARIES})
      message(STATUS "OpenCASCADE_LINK_DIRECTORY = " ${OpenCASCADE_LINK_DIRECTORY})
      include_directories(${OpenCASCADE_INCLUDE_DIR})
    ELSE(OpenCASCADE_FOUND)
      MESSAGE(STATUS "OpenCascade Libraries not found in standard paths.")
    ENDIF(OpenCASCADE_FOUND)

  endif(OCE_FOUND)
  SET(HAVE_MECHANISMS TRUE)
endif()

# -- VTK --
IF(WITH_VTK)
  COMPILE_WITH(VTK)
  IF(VTK_FOUND)
    MESSAGE(STATUS "Found vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
    SET(HAVE_VTK TRUE)
    IF(VTK_USE_FILE)
      INCLUDE(${VTK_USE_FILE})
    ENDIF()
  ENDIF()
ENDIF()

# -- FreeCAD --
# For python bindings
if(WITH_FREECAD)
  compile_with(FreeCAD COMPONENTS Part REQUIRED)
endif()


# -- HDF5 --
# For logging in Numerics
IF(WITH_HDF5)
  COMPILE_WITH(HDF5 REQUIRED)
ENDIF(WITH_HDF5)

#
# -- Serialization --
#
include(serialization_vector_test)
if(WITH_SERIALIZATION)
  COMPILE_WITH(Boost 1.47 COMPONENTS serialization filesystem REQUIRED)
  TEST_SERIALIZATION_VECTOR_BUG()
endif()

# -- Python bindings --
if(WITH_PYTHON_WRAPPER)
  compile_With(Numpy REQUIRED)
  find_package(SWIG 2.0.7 REQUIRED)
  include(${SWIG_USE_FILE})
  include(FindPythonModule)
  if(NOT NO_RUNTIME_BUILD_DEP)
    find_python_module(scipy REQUIRED) # for sparse
    IF(scipy_VERSION VERSION_LESS "0.14.0")
      SET(SICONOS_FORCE_NPY_INT32 TRUE)
      MESSAGE(STATUS "Old version of scipy detected. Support for 64 bits int was added in 0.14.0, forcing 32 bits int for python bindings")
    ENDIF(scipy_VERSION VERSION_LESS "0.14.0")
  ENDIF(NOT NO_RUNTIME_BUILD_DEP)
  find_python_module(pyhull)
  if(NOT pyhull_FOUND)
    message(STATUS "Warning, python pyhull package not found. Some examples may not run properly. Try pip install pyhull.")
  endif()
endif()


#
# Fedora13 https://fedoraproject.org/wiki/UnderstandingDSOLinkChange
if(UNIX)
  # add -lm to linker
  CHECK_C_COMPILER_FLAG("-lm" C_HAVE_LINKER_M)
  if(C_HAVE_LINKER_M)
    set(SICONOS_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES} "m" CACHE INTERNAL "List of external libraries")
  endif()
endif()


# SiconosConfig.h generation and include
if(EXISTS ${CMAKE_SOURCE_DIR}/config.h.cmake)
  configure_file(${CMAKE_SOURCE_DIR}/config.h.cmake
    ${CMAKE_BINARY_DIR}/SiconosConfig.h)
endif()

set(${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
  ${${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES}
  ${CMAKE_BINARY_DIR}
  CACHE INTERNAL "Include directories for external dependencies.")
