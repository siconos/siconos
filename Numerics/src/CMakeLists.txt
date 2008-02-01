#
# Here comes Numerics
#
PROJECT(SiconosNumerics)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosNumerics_VERSION "${VERSION}")

#
# Published headers
#
SET(SiconosNumerics_HDRS
  Solver_Options.h
  SiconosNumerics.h
)

#
# Project directories 
#
SET(SiconosNumerics_DIRS
  netlib/odepack
  netlib/dftemplates
  NonSmoothSolvers
  NonSmoothSolvers/LCP
  NonSmoothSolvers/MLCP
#  NonSmoothSolvers/InterfaceToPathFerris
  NonSmoothSolvers/pfc3D
  NonSmoothSolvers/QP
  NonSmoothSolvers/Relay
  NonSmoothSolvers/FrictionContact2D
  NonSmoothSolvers/FrictionContact3D
  NonSmoothSolvers/NCP
  NonSmoothSolvers/NSSTools
  utils)

#
# Unstable sources
#
SET(SiconosNumerics_Unstable_SRCS
  NonSmoothSolvers/NSSTools/freeSpBlMat.c
  NonSmoothSolvers/test/main_pfc_2D.c
  NonSmoothSolvers/test/main_lcp.c
  NonSmoothSolvers/test/main_dfc_2D.c
  NonSmoothSolvers/test/main_pr.c
  NonSmoothSolvers/test/main_lcp_new.c
  NonSmoothSolvers/test/main_dr.c)

#
# Installation layout
#
SET(SiconosNumerics_include_dir include/Siconos)

#
# compilation setup 
#
SET(SiconosNumerics_LINK ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

LIBRARY_PROJECT_SETUP()

#ADD_DEFINITIONS(-std=c99) : warning with f95
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

# F95 Cmake Bug
CONFIGURE_FILE(SiconosNumerics.h  DPCG.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DPCGS.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DSOLBT.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DUSOL.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  if.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  in.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  the.mod.proxy COPYONLY) 

# test directories
SUBDIRS(NonSmoothSolvers/test)
