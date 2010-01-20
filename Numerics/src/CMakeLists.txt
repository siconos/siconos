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
  NonSmoothSolvers/LCP/test-utils
  NonSmoothSolvers/MLCP
  NonSmoothSolvers/LinearSystem
  NonSmoothSolvers/InterfaceToPathFerris
  NonSmoothSolvers/QP
  NonSmoothSolvers/Relay
  NonSmoothSolvers/Relay/test-utils
  NonSmoothSolvers/FrictionContact
  NonSmoothSolvers/FrictionContact/test-utils
  NonSmoothSolvers/NCP			
  NonSmoothSolvers/MCP
  NonSmoothTools
  NonSmoothTools/test-utils
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

IF(C_HAVE_WALL)
  APPEND_FLAGS(CMAKE_C_FLAGS "-Wall")
ENDIF(C_HAVE_WALL)

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
NEW_TEST(odepacktest1 DLSODES-test.f)
NEW_TEST(odepacktest2 DLSODAR-test.f)
NEW_TEST(odepacktest3 DLSODI-test.f)
NEW_TEST(odepacktest4 DLSODPK-test.f)
NEW_TEST(odepacktest5 DLSODA-test.f)
NEW_TEST(odepacktest6 DLSODE-test.f)
NEW_TEST(odepacktest7 DLSODIS-test.f)
NEW_TEST(odepacktest8 DLSODKR-test.f)
NEW_TEST(odepacktest9 DLSOIBT-test.f)
NEW_TEST(odepacktest10 test-funcC-inC.cpp funC.cpp)
END_TEST()
 

BEGIN_TEST(NonSmoothTools/test) 
#NEW_TEST(NumericsMatrixTest main_NumericsMatrix.c)
NEW_TEST(NumericsMatrixTest0 NumericsMatrix_test0.c)
NEW_TEST(NumericsMatrixTest1 NumericsMatrix_test1.c)
NEW_TEST(NumericsMatrixTest2 NumericsMatrix_test2.c)
NEW_TEST(NumericsMatrixTest3 NumericsMatrix_test3.c)
NEW_TEST(NumericsMatrixTest4 NumericsMatrix_test4.c)
NEW_TEST(NumericsMatrixTest5 NumericsMatrix_test5.c)
NEW_TEST(NumericsMatrixTest6 NumericsMatrix_test6.c)
END_TEST()

BEGIN_TEST(NonSmoothSolvers/LCP/test)

NEW_TEST(LCPtest1 lcp_test1.c)
NEW_TEST(LCPtest2 lcp_test2.c)
NEW_TEST(LCPtest3 lcp_test3.c)
NEW_TEST(LCPtest4 lcp_test4.c)
NEW_TEST(LCPtest5 lcp_test5.c)
NEW_TEST(LCPtest6 lcp_test6.c)
NEW_TEST(LCPtest7 lcp_test7.c)
NEW_TEST(LCPtest8 lcp_test8.c)
NEW_TEST(LCPtest9 lcp_test9.c)
#SET(LCPtest10_PROPERTIES WILL_FAIL TRUE)
#NEW_TEST(LCPtest10 lcp_test10.c)
IF(PathFerris_FOUND)
NEW_TEST(LCPtest11 lcp_test11.c)
ENDIF(PathFerris_FOUND)

NEW_TEST(LCPtest20 lcp_test20.c)
NEW_TEST(LCPtest21 lcp_test21.c)
SET(LCPtest22_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(LCPtest22 lcp_test22.c)

NEW_TEST(LCPtest23 lcp_test23.c)
NEW_TEST(LCPtest24 lcp_test24.c)
NEW_TEST(LCPtest25 lcp_test25.c)
NEW_TEST(LCPtest26 lcp_test26.c)
NEW_TEST(LCPtest27 lcp_test27.c)
NEW_TEST(LCPtest28 lcp_test28.c)
NEW_TEST(LCPtest29 lcp_test29.c)
IF(PathFerris_FOUND)
NEW_TEST(LCPtest30 lcp_test30.c)
ENDIF(PathFerris_FOUND)

NEW_TEST(LCPtest40 lcp_test40.c)
NEW_TEST(LCPtest41 lcp_test41.c)
NEW_TEST(LCPtest42 lcp_test42.c)
NEW_TEST(LCPtest43 lcp_test43.c)
NEW_TEST(LCPtest44 lcp_test44.c)
NEW_TEST(LCPtest45 lcp_test45.c)
NEW_TEST(LCPtest46 lcp_test46.c)
NEW_TEST(LCPtest47 lcp_test47.c)
NEW_TEST(LCPtest48 lcp_test48.c)



NEW_TEST(LCPtest50 lcp_test50.c)
NEW_TEST(LCPtest51 lcp_test51.c)
NEW_TEST(LCPtest52 lcp_test52.c)
NEW_TEST(LCPtest53 lcp_test53.c)
NEW_TEST(LCPtest54 lcp_test54.c)
NEW_TEST(LCPtest55 lcp_test55.c)
NEW_TEST(LCPtest56 lcp_test56.c)
NEW_TEST(LCPtest57 lcp_test57.c)
NEW_TEST(LCPtest58 lcp_test58.c)

NEW_TEST(LCPtest60 lcp_test60.c)
NEW_TEST(LCPtest61 lcp_test61.c)
NEW_TEST(LCPtest62 lcp_test62.c)
NEW_TEST(LCPtest63 lcp_test63.c)
NEW_TEST(LCPtest64 lcp_test64.c)
NEW_TEST(LCPtest65 lcp_test65.c)
NEW_TEST(LCPtest66 lcp_test66.c)
NEW_TEST(LCPtest67 lcp_test67.c)
NEW_TEST(LCPtest68 lcp_test68.c)

NEW_TEST(LCPtest70 lcp_test70.c)
NEW_TEST(LCPtest71 lcp_test71.c)
SET(LCPtest72_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(LCPtest72 lcp_test72.c)
NEW_TEST(LCPtest73 lcp_test73.c)
NEW_TEST(LCPtest74 lcp_test74.c)
NEW_TEST(LCPtest75 lcp_test75.c)
NEW_TEST(LCPtest76 lcp_test76.c)
NEW_TEST(LCPtest77 lcp_test77.c)
NEW_TEST(LCPtest78 lcp_test78.c)

NEW_TEST(LCPtest80 lcp_test80.c)
NEW_TEST(LCPtest81 lcp_test81.c)
NEW_TEST(LCPtest82 lcp_test82.c)
NEW_TEST(LCPtest83 lcp_test83.c)
NEW_TEST(LCPtest84 lcp_test84.c)
NEW_TEST(LCPtest85 lcp_test85.c)
NEW_TEST(LCPtest86 lcp_test86.c)
NEW_TEST(LCPtest87 lcp_test87.c)
NEW_TEST(LCPtest88 lcp_test88.c)

END_TEST()

BEGIN_TEST(NonSmoothSolvers/Relay/test)

NEW_TEST(Relaytest1 relay_test1.c)
NEW_TEST(Relaytest2 relay_test2.c)

IF(PathFerris_FOUND)
NEW_TEST(Relaytest3 relay_test3.c)
ENDIF(PathFerris_FOUND)

NEW_TEST(Relaytest10 relay_test10.c)
NEW_TEST(Relaytest11 relay_test11.c)
NEW_TEST(Relaytest12 relay_test12.c)
END_TEST()

BEGIN_TEST(NonSmoothSolvers/MLCP/test)
NEW_TEST(MLCPtest main_mlcp.cpp)
END_TEST()
BEGIN_TEST(NonSmoothSolvers/FrictionContact/test)
 
#NEW_TEST(FrictionContact_Problemtest main_FC3D.c)
 NEW_TEST(FC3D_DefaultSolverOptionstest FrictionContact3D_DefaultSolverOptions_test.c)


# (see NonSmoothSolvers/FrictionContact/test/README for short details)
NEW_TEST(FC3Dtest1 FrictionContact3D_test1.c)
NEW_TEST(FC3Dtest2 FrictionContact3D_test2.c)
NEW_TEST(FC3Dtest3 FrictionContact3D_test3.c)
NEW_TEST(FC3Dtest4 FrictionContact3D_test4.c)
SET(FC3Dtest5_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest5 FrictionContact3D_test5.c)

IF(PathFerris_FOUND)
NEW_TEST(FC3Dtest6 FrictionContact3D_test6.c)
ENDIF(PathFerris_FOUND)

SET(FC3Dtest7_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest7 FrictionContact3D_test7.c)





NEW_TEST(FC3Dtest10 FrictionContact3D_test10.c)
NEW_TEST(FC3Dtest11 FrictionContact3D_test11.c)
NEW_TEST(FC3Dtest12 FrictionContact3D_test12.c)
SET(FC3Dtest13_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest13 FrictionContact3D_test13.c)

NEW_TEST(FC3Dtest20 FrictionContact3D_test20.c)
NEW_TEST(FC3Dtest21 FrictionContact3D_test21.c)
NEW_TEST(FC3Dtest22 FrictionContact3D_test22.c)
SET(FC3Dtest23_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest23 FrictionContact3D_test23.c)

NEW_TEST(FC3Dtest30 FrictionContact3D_test30.c)
NEW_TEST(FC3Dtest31 FrictionContact3D_test31.c)
NEW_TEST(FC3Dtest32 FrictionContact3D_test32.c)

SET(FC3Dtest50_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest50 FrictionContact3D_test50.c)


NEW_TEST(FC3Dtest60 FrictionContact3D_test60.c)

SET(FC3Dtest50_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest70 FrictionContact3D_test70.c)


NEW_TEST(FC3Dtest100 FrictionContact3D_test100.c)
SET(FC3Dtest110_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest110 FrictionContact3D_test110.c)
SET(FC3Dtest120_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest120 FrictionContact3D_test120.c)
SET(FC3Dtest130_PROPERTIES WILL_FAIL TRUE)
NEW_TEST(FC3Dtest130 FrictionContact3D_test130.c)


## test 2D dense on two differents files
#NEW_TEST(FrictionContact2D_Problemtest1 FrictionContact2D_test1.c)
#NEW_TEST(FrictionContact2D_Problemtest2 FrictionContact2D_test2.c)
#NEW_TEST(FrictionContact2D_Problemtest3 FrictionContact2D_test3.c)


#NEW_TEST(FrictionContact2D_Problemtest10 FrictionContact2D_test10.c)
#NEW_TEST(FrictionContact2D_Problemtest11 FrictionContact2D_test11.c)
#NEW_TEST(FrictionContact2D_Problemtest12 FrictionContact2D_test12.c)


## test 2D sparse on 4 differents files
#NEW_TEST(FrictionContact2D_Problemtest20 FrictionContact2D_test20.c)
#NEW_TEST(FrictionContact2D_Problemtest21 FrictionContact2D_test21.c)
#NEW_TEST(FrictionContact2D_Problemtest22 FrictionContact2D_test22.c)
#NEW_TEST(FrictionContact2D_Problemtest23 FrictionContact2D_test23.c)

NEW_TEST(PFC3D_test1 PrimalFrictionContact3D_test1.c)
NEW_TEST(PFC3D_test2 PrimalFrictionContact3D_test2.c)
NEW_TEST(PFC3D_test3 PrimalFrictionContact3D_test3.c)
NEW_TEST(PFC3D_test4 PrimalFrictionContact3D_test4.c)
NEW_TEST(PFC3D_test5 PrimalFrictionContact3D_test5.c)

END_TEST()


#BEGIN_TEST(utils/test)
#NEW_TEST(pinvtest testpinv.c)
#END_TEST()




IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)



IF(PRINT_VAR)
  MESSAGE(STATUS)
  PRINT_VAR(${PRINT_VAR})
  MESSAGE(STATUS)
ENDIF(PRINT_VAR)
