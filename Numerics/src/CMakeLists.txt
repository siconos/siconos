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
  ..
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
# Atlas is not required
COMPILE_WITH(ATLAS)

COMPILE_WITH(LAPACK REQUIRED)
COMPILE_WITH(BLAS REQUIRED)

IF(C_HAVE_C99)
  APPEND_FLAGS(CMAKE_C_FLAGS "-std=c99")
ENDIF(C_HAVE_C99)

LIBRARY_PROJECT_SETUP()



# F95 Cmake Bug
CONFIGURE_FILE(SiconosNumerics.h  DPCG.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DPCGS.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DSOLBT.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  DUSOL.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  if.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  in.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  the.mod.proxy COPYONLY) 
CONFIGURE_FILE(SiconosNumerics.h  unit.mod.proxy COPYONLY) 

# test directories
#SUBDIRS(netlib/odepack/test)
#SUBDIRS(NonSmoothSolvers/test)

TEST_DIRECTORY(netlib/odepack/test)

SET_TEST(test1 DLSODES-test.f)
SET_TEST(test2 DLSODAR-test.f)
SET_TEST(test3 DLSODI-test.f)
SET_TEST(test4 DLSODPK-test.f)
SET_TEST(test5 DLSODA-test.f)
SET_TEST(test6 DLSODE-test.f)
SET_TEST(test7 DLSODIS-test.f)
SET_TEST(test8 DLSODKR-test.f)
SET_TEST(test9 DLSOIBT-test.f)
SET_TEST(test10 test-funcC-inC.cpp funC.cpp)

TARGET_LINK_LIBRARIES(netlib/odepack/test/test10 ${FORTRAN_LIBRARIES})

TEST_DIRECTORY(NonSmoothSolvers/test)
#SET_TEST(NumericsMatrixTest main_NumericsMatrix.c)
#SET_TEST(LCPtest main_lcp.c)
#SET_TEST(MLCPtest main_mlcp.c)
#SET_TESTS_PROPERTIES(NonSmoothSolvers/test/MLCPtest PROPERTIES WILL_FAIL TRUE)

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)
