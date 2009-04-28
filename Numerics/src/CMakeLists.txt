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
  NonSmoothSolvers/MCP
  NonSmoothSolvers/NSSTools
  utils)

#
# Unstable sources
#
SET(SiconosNumerics_Unstable_SRCS
)

#
# Installation layout
#
SET(SiconosNumerics_INSTALL_INCLUDE_DIR include/Siconos/Numerics)

# NumericsConfig.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/NumericsConfig.h DESTINATION include/Siconos/Numerics)

#
# compilation setup 
#
# Atlas is not required
COMPILE_WITH(ATLAS)

COMPILE_WITH(LAPACK REQUIRED COMPLETE)
COMPILE_WITH(BLAS REQUIRED)

COMPILE_WITH(PathFerris)
COMPILE_WITH(MlcpSimplex)
COMPILE_WITH(Pthread)
COMPILE_WITH(Cplex)
IF(C_HAVE_C99)
  APPEND_FLAGS(CMAKE_C_FLAGS "-std=c99")
ENDIF(C_HAVE_C99)

APPEND_FLAGS(CMAKE_CXX_FLAGS "-D_NUMERICS_INTERNAL_CXX_")

# CXX may be needed with MlcpSimplex
SET(SiconosNumerics_LINKER_LANGUAGE C)
IF(MlcpSimplex_FOUND) 
SET(SiconosNumerics_LINKER_LANGUAGE CXX)
ENDIF(MlcpSimplex_FOUND)


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

# The library tests

# wrapper are not needed
SET(TEST_WRAP)

BEGIN_TEST(netlib/odepack/test)
 
NEW_TEST(test1 DLSODES-test.f)
NEW_TEST(test2 DLSODAR-test.f)
NEW_TEST(test3 DLSODI-test.f)
NEW_TEST(test4 DLSODPK-test.f)
NEW_TEST(test5 DLSODA-test.f)
NEW_TEST(test6 DLSODE-test.f)
NEW_TEST(test7 DLSODIS-test.f)
NEW_TEST(test8 DLSODKR-test.f)
NEW_TEST(test9 DLSOIBT-test.f)
NEW_TEST(test10 test-funcC-inC.cpp funC.cpp)
 
END_TEST()
 
BEGIN_TEST(NonSmoothSolvers/test)
 
NEW_TEST(NumericsMatrixTest main_NumericsMatrix.c)
NEW_TEST(LCPtest main_lcp.c)

##SET(MLCPtest_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(MLCPtest main_mlcp.cpp)

NEW_TEST(Relaytest main_relay.c)

END_TEST()

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)



IF(PRINT_VAR)
  MESSAGE(STATUS)
  PRINT_VAR(${PRINT_VAR})
  MESSAGE(STATUS)
ENDIF(PRINT_VAR)
