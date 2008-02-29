#
# Here comes Numerics
#
PROJECT(SiconosNumerics C CXX Fortran)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosNumerics_VERSION "${VERSION}")

#
# Published headers
#
# if <PROJECT>_HDRS unset -> all headers under <PROJECT>_DIRS

#
# Project directories 
#
SET(SiconosNumerics_DIRS
  .
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
  utils)

#
# Unstable sources
#
SET(SiconosNumerics_Unstable_SRCS
  NonSmoothSolvers/LCP/lcp_solver_block_pred_vec.c
  NonSmoothSolvers/LCP/lcp_solver_pred.c
  NonSmoothSolvers/InterfaceToPathFerris/SimpleLCP.h
  NonSmoothSolvers/InterfaceToPathFerris/SimpleLCP_Path.c
)

#
# Installation layout
#
SET(SiconosNumerics_INSTALL_INCLUDE_DIR include/Siconos/Numerics)

#
# compilation setup 
#
SET(SiconosNumerics_LINK_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FORTRAN_LIBRARIES})

LIBRARY_PROJECT_SETUP()

IF(C_HAVE_C99)
  APPEND_FLAG(CMAKE_C_FLAGS "-std=c99")
ENDIF(C_HAVE_C99)

# F95 Cmake Bug
CONFIGURE_FILE(SiconosNumerics.h  DPCG.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DPCGS.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DSOLBT.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DUSOL.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  if.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  in.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  the.mod.proxy COPYONLY) 

# test directories
SUBDIRS(netlib/odepack/test)
SUBDIRS(NonSmoothSolvers/test)
