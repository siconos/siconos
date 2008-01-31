#
# Here comes Numerics
#
PROJECT(Numerics)

#
# Major, Minor and patch versions defined in top-level CMakeLists.txt
#
SET(Numerics_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")

#
# Published headers
#
SET(Numerics_HDRS
  Solver_Options.h
  SiconosNumerics.h
)

#
# Project directories 
#
SET(Numerics_DIRS
  netlib/odepack
  netlib/dftemplates
  NonSmoothSolvers
  NonSmoothSolvers/LCP
  NonSmoothSolvers/MLCP
  NonSmoothSolvers/InterfaceToPathFerris
  NonSmoothSolvers/pfc3D
  NonSmoothSolvers/QP
  NonSmoothSolvers/Relay
  NonSmoothSolvers/FrictionContact2D
  NonSmoothSolvers/FrictionContact3D
  NonSmoothSolvers/NCP
  NonSmoothSolvers/NSSTools
  NonSmoothSolvers/test
  utils)

#
# Unstable sources
#
SET(Numerics_Unstable_SRCS
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
SET(Numerics_include_dir include/Siconos)

#
# compilation setup 
#
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
