#================================================================
#
# Look for siconos dependencies, common to several components.
#
#================================================================

# --- Blas Lapack ---
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
include(gams_setup)
# the following library may be found in a GAMS install -> we search for them
# after
compile_with(PathFerris)
compile_with(Cplex)
compile_with(LpSolve)
if(LpSolve_FOUND)
  set(HAS_ONE_LP_SOLVER TRUE)
  set(HAS_EXTREME_POINT_ALGO TRUE)
  set(WITH_LPSOLVE TRUE)
endif(LpSolve_FOUND)


# --- Mumps ---
if(WITH_MUMPS)
  compile_with(MPI REQUIRED)
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
compile_with(Boost REQUIRED)
compile_with(GMP REQUIRED)

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
    COMPILE_WITH(OCE REQUIRED)
    SET(HAVE_OCC TRUE)
    LIST(REMOVE_ITEM SICONOS_LINK_LIBRARIES DRAWEXE)
  endif()
ENDIF()

# --- Mechanisms (Saladyn?) ---
IF(WITH_MECHANISMS)
  SET(OCE_TOOLKITS "TKernel"  "TKMath" "TKService" "TKV3d"  "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")

  message(STATUS "Searching for OCE ....")
  compile_with(OCE COMPONENTS ${OCE_TOOLKITS})
  LIST(REMOVE_ITEM SICONOS_LINK_LIBRARIES DRAWEXE)

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

# -- Python bindings --
if(WITH_PYTHON_BINDINGS)
  compile_With(Python_Numpy REQUIRED)
  find_package(SWIG 2.0.7 REQUIRED)
  include(${SWIG_USE_FILE})
endif()