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
  hairer/
  netlib/odepack
  netlib/dftemplates
  NonSmoothSolvers
  NonSmoothSolvers/LCP
  NonSmoothSolvers/LCP/test-utils
  NonSmoothSolvers/MLCP
  NonSmoothSolvers/LinearSystem
  NonSmoothSolvers/LinearSystem/test-utils
  NonSmoothSolvers/InterfaceToPathFerris
  NonSmoothSolvers/QP
  NonSmoothSolvers/Relay
  NonSmoothSolvers/Relay/test-utils
  NonSmoothSolvers/FrictionContact
  NonSmoothSolvers/FrictionContact/test-utils
  NonSmoothSolvers/NCP			
  NonSmoothSolvers/MCP
  NonSmoothSolvers/GenericMechanical
  NonSmoothSolvers/GenericMechanical/test-utils
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

if(USE_MKL)
  compile_with(MKL REQUIRED)
else()
 if(NOT WITHOUT_ATLAS)
  COMPILE_WITH(ATLAS)
 endif(NOT WITHOUT_ATLAS)
  COMPILE_WITH(LAPACK REQUIRED COMPLETE)
  COMPILE_WITH(BLAS REQUIRED)
endif()

COMPILE_WITH(PathFerris)
COMPILE_WITH(MlcpSimplex)
COMPILE_WITH(Pthread)
COMPILE_WITH(Cplex)

IF(WITH_MUMPS)
  COMPILE_WITH(MPI REQUIRED)
  COMPILE_WITH(MUMPS REQUIRED)
ENDIF(WITH_MUMPS)

IF(WITH_FCLIB)
  COMPILE_WITH(FCLIB REQUIRED)   
  IF(FCLIB_NOTFOUND)
    # try the package stuff
    # need FCLib_DIR !!
    COMPILE_WITH(FCLib 1.0 REQUIRED)
  ENDIF()
  #  APPEND_FLAGS(CMAKE_C_FLAGS "-DWITH_FCLIB")	
ENDIF()


IF(C_HAVE_C99)
  APPEND_FLAGS(CMAKE_C_FLAGS "-std=c99")
ENDIF(C_HAVE_C99)

IF(C_HAVE_WALL)
  APPEND_FLAGS(CMAKE_C_FLAGS "-Wall")
ENDIF(C_HAVE_WALL)

IF(C_HAVE_WEXTRA)
 APPEND_FLAGS(CMAKE_CXX_FLAGS "-Wextra -Wno-unused-parameter")
ENDIF(C_HAVE_WEXTRA)

#
# Fedora13 https://fedoraproject.org/wiki/UnderstandingDSOLinkChange
#
IF(C_HAVE_LINKER_M)
  FIND_LIBRARY(M_FOUND "m")
  SET(M_LIBRARIES ${M_FOUND})
  COMPILE_WITH(M STANDARD)
ENDIF(C_HAVE_LINKER_M)  



APPEND_FLAGS(CMAKE_CXX_FLAGS "-D_NUMERICS_INTERNAL_CXX_")

# CXX may be needed with MlcpSimplex
SET(SiconosNumerics_LINKER_LANGUAGE C)
IF(MlcpSimplex_FOUND) 
SET(SiconosNumerics_LINKER_LANGUAGE CXX)
ENDIF(MlcpSimplex_FOUND)

# do not skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

LIBRARY_PROJECT_SETUP()

# F95 Cmake Bug. No more needed? 
#CONFIGURE_FILE(SiconosNumerics.h  DPCG.mod.proxy COPYONLY) 
#CONFIGURE_FILE(SiconosNumerics.h  DPCGS.mod.proxy COPYONLY) 
#CONFIGURE_FILE(SiconosNumerics.h  DSOLBT.mod.proxy COPYONLY) 
#CONFIGURE_FILE(SiconosNumerics.h  DUSOL.mod.proxy COPYONLY) 
#CONFIGURE_FILE(SiconosNumerics.h  if.mod.proxy COPYONLY) 
#CONFIGURE_FILE(SiconosNumerics.h  in.mod.proxy COPYONLY) 
#CONFIGURE_FILE(SiconosNumerics.h  the.mod.proxy COPYONLY) 
#CONFIGURE_FILE(SiconosNumerics.h  unit.mod.proxy COPYONLY) 

# The library tests

# wrapper are not needed
SET(TEST_WRAP)

IF(WITH_TESTING)
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

  BEGIN_TEST(hairer/test)
  NEW_TEST(dr_iso1 dr_iso.f)
  NEW_TEST(dr_iso1sp dr_isosp.f)
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
  NEW_TEST(SBMTest1 SBM_test1.c)
  NEW_TEST(SBMTest2 SBM_test2.c)
  NEW_TEST(SBMTest3 SBM_test3.c)
  NEW_TEST(SBMTest4  SBM_test4.c)
 END_TEST()
  
  BEGIN_TEST(NonSmoothSolvers/LCP/test)
  NEW_LCP_TEST(SICONOS_LCP_PGS lcp_mmc.dat)
  NEW_LCP_TEST(SICONOS_LCP_RPGS lcp_mmc.dat)
  NEW_LCP_TEST(SICONOS_LCP_CPG lcp_mmc.dat)
  NEW_LCP_TEST(SICONOS_LCP_LEMKE lcp_mmc.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN lcp_mmc.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN_W lcp_mmc.dat)
  NEW_LCP_TEST(SICONOS_LCP_QP lcp_mmc.dat)
  NEW_LCP_TEST(SICONOS_LCP_NSQP lcp_mmc.dat)
  #SET(test-LCP_NEWTONMIN-lcp_mmc_PROPERTIES WILL_FAIL TRUE)
  #NEW_LCP_TEST(SICONOS_LCP_NEWTONMIN lcp_mmc.dat)
  
  IF(PathFerris_FOUND)
    NEW_LCP_TEST(SICONOS_LCP_PATH lcp_mmc.dat)
  ENDIF(PathFerris_FOUND)
 
  NEW_LCP_TEST(SICONOS_LCP_PGS lcp_deudeu.dat)
  NEW_LCP_TEST(SICONOS_LCP_RPGS lcp_deudeu.dat)
  SET(test-LCP_CPG-lcp_deudeu_PROPERTIES WILL_FAIL TRUE)
  NEW_LCP_TEST(SICONOS_LCP_CPG lcp_deudeu.dat)
  NEW_LCP_TEST(SICONOS_LCP_LEMKE lcp_deudeu.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN lcp_deudeu.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN_W lcp_deudeu.dat)
  NEW_LCP_TEST(SICONOS_LCP_QP lcp_deudeu.dat)
  NEW_LCP_TEST(SICONOS_LCP_NSQP lcp_deudeu.dat)
  NEW_LCP_TEST(SICONOS_LCP_NEWTONMIN lcp_deudeu.dat)
  IF(PathFerris_FOUND)
    NEW_LCP_TEST(SICONOS_LCP_PATH lcp_deudeu.dat)
  ENDIF(PathFerris_FOUND)


  NEW_LCP_TEST(SICONOS_LCP_PGS lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_RPGS lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_LEMKE lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_LATIN lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_LATIN_W lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_QP lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_NSQP lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_NEWTONMIN lcp_deudeu_block.dat 1)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_deudeu_block.dat 1)

  NEW_LCP_TEST(SICONOS_LCP_PGS lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_RPGS lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_CPG lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_LEMKE lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN_W lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_QP lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_NSQP lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_NEWTONMIN lcp_trivial.dat)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_trivial.dat)

  NEW_LCP_TEST(SICONOS_LCP_PGS lcp_ortiz.dat)
  SET(test-LCP_RPGS-lcp_ortiz WILL_FAIL TRUE)
  NEW_LCP_TEST(SICONOS_LCP_RPGS lcp_ortiz.dat)
  SET(test-LCP_CPG-lcp_ortiz_PROPERTIES WILL_FAIL TRUE)
  NEW_LCP_TEST(SICONOS_LCP_CPG lcp_ortiz.dat)
  NEW_LCP_TEST(SICONOS_LCP_LEMKE lcp_ortiz.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN lcp_ortiz.dat)
  NEW_LCP_TEST(SICONOS_LCP_LATIN_W lcp_ortiz.dat)
  NEW_LCP_TEST(SICONOS_LCP_QP lcp_ortiz.dat)
  NEW_LCP_TEST(SICONOS_LCP_NSQP lcp_ortiz.dat)
  NEW_LCP_TEST(SICONOS_LCP_NEWTONMIN lcp_ortiz.dat)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_ortiz.dat)

  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol.dat)
  SET(test-LCP_ENUM-lcp_Pang_isolated_sol_perturbed_PROPERTIES WILL_FAIL TRUE)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol_perturbed.dat)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_inf_sol_perturbed.dat)
  NEW_TEST(LCP_DefaultSolverOptionstest LinearComplementarity_DefaultSolverOptions_test.c)
  
  
  END_TEST()

  BEGIN_TEST(NonSmoothSolvers/LinearSystem/test)

  NEW_LS_TEST(SICONOS_LS_0 ls_trivial.dat)
  SET(test-LS_0-ls_inf_sol_perturbed_PROPERTIES WILL_FAIL TRUE)
  NEW_LS_TEST(SICONOS_LS_0 ls_inf_sol_perturbed.dat)
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
  NEW_TEST(Relaytest20 relay_test20.c)
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
    SET(FC3Dtest6_PROPERTIES WILL_FAIL TRUE)
  ENDIF(PathFerris_FOUND)
  
  SET(FC3Dtest7_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(FC3Dtest7 FrictionContact3D_test7.c)
  
  



  NEW_TEST(FC3Dtest10 FrictionContact3D_test10.c)
  NEW_TEST(FC3Dtest11 FrictionContact3D_test11.c)
  NEW_TEST(FC3Dtest12 FrictionContact3D_test12.c)
  NEW_TEST(FC3Dtest13 FrictionContact3D_test13.c)
  NEW_TEST(FC3Dtest14 FrictionContact3D_test14.c)

  SET(NSGS_TOL 1e-6)
  SET(NSGS_NB_IT 10000)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Example1_Fc3D.dat 1e-5 ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_AlartCurnierNewton 1e-18 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover1039.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover1040.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover1041.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover3865.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4144.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4396.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4493.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4516.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4609.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4613.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4622.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover9770.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover11035.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover11211.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_QUARTIC 1e-6 10)
  SET(FC3Dtest20_PROPERTIES WILL_FAIL TRUE)

  
  NEW_TEST(FC3Dtest20 FrictionContact3D_test20.c)
  NEW_TEST(FC3Dtest21 FrictionContact3D_test21.c)
  SET(FC3Dtest22_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(FC3Dtest22 FrictionContact3D_test22.c)
  NEW_TEST(FC3Dtest23 FrictionContact3D_test23.c)
  
  NEW_TEST(FC3Dtest30 FrictionContact3D_test30.c)
  NEW_TEST(FC3Dtest31 FrictionContact3D_test31.c)
  NEW_TEST(FC3Dtest32 FrictionContact3D_test32.c)
  
  SET(FC3Dtest50_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(FC3Dtest50 FrictionContact3D_test50.c)
  

  NEW_TEST(FC3Dtest60 FrictionContact3D_test60.c)

  NEW_TEST(FC3Dtest70 FrictionContact3D_test70.c)
  NEW_TEST(FC3Dtest71 FrictionContact3D_test71.c)	
  NEW_TEST(FC3Dtest72 FrictionContact3D_test72.c)	

  NEW_TEST(FC3Dtest100 FrictionContact3D_test100.c)
  NEW_TEST(FC3Dtest110 FrictionContact3D_test110.c)
  SET(FC3Dtest120_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(FC3Dtest120 FrictionContact3D_test120.c)
  NEW_TEST(FC3Dtest130 FrictionContact3D_test130.c)

  NEW_TEST(FC3Dtest200 FrictionContact3D_test200.c)
  NEW_TEST(FC3Dtest210 FrictionContact3D_test210.c)

  # Global Alart Curnier + Rover ok.
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Example1_Fc3D.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover1039.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover1040.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover1041.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover3865.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover4144.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover4396.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover4493.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover4516.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover4609.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover4613.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover4622.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover9770.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover11035.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC Rover11211.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_QUARTIC FrictionContact3D_1c.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_QUARTIC FrictionContact3D_RR_1c.dat)

  # Newton Euler 3D Spheres : some failures with global Alart Curnier
  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC NESpheres_10_1.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_GLOBALAC NESpheres_30_1.dat)

  ## test from the Rover Example
  NEW_TEST(FC3Dtest501 FrictionContact3D_Rover2.c) ## test ok without LineSearch
  NEW_TEST(FC3Dtest503 FrictionContact3D_Rover3.c) ## LineSearch is needed 


  ## test 2D dense on two differents files
  NEW_TEST(FrictionContact2D_Problemtest1 FrictionContact2D_test1.c)
  NEW_TEST(FrictionContact2D_Problemtest2 FrictionContact2D_test2.c)
  NEW_TEST(FrictionContact2D_Problemtest3 FrictionContact2D_test3.c)
  
  
  NEW_TEST(FrictionContact2D_Problemtest10 FrictionContact2D_test10.c)
  NEW_TEST(FrictionContact2D_Problemtest11 FrictionContact2D_test11.c)
  NEW_TEST(FrictionContact2D_Problemtest12 FrictionContact2D_test12.c)
  

  ## test 2D sparse on 4 differents files
  NEW_TEST(FrictionContact2D_Problemtest20 FrictionContact2D_test20.c)
  NEW_TEST(FrictionContact2D_Problemtest21 FrictionContact2D_test21.c)
  NEW_TEST(FrictionContact2D_Problemtest22 FrictionContact2D_test22.c)
  NEW_TEST(FrictionContact2D_Problemtest23 FrictionContact2D_test23.c)

  ## test 2D dense with Lemke NSGS failed on it
  NEW_TEST(FrictionContact2D_Problemtest30 FrictionContact2D_test30.c)
  NEW_TEST(FrictionContact2D_Problemtest31 FrictionContact2D_test31.c)
  #SET(FrictionContact2D_Problemtest32 WILL_FAIL TRUE)
  #NEW_TEST(FrictionContact2D_Problemtest32 FrictionContact2D_test32.c)

  ## test 2D dense with Enum Lemke failed on 41 !!
  NEW_TEST(FrictionContact2D_Problemtest40 FrictionContact2D_test40.c)
  NEW_TEST(FrictionContact2D_Problemtest41 FrictionContact2D_test41.c)


  NEW_TEST(PFC3D_test1 PrimalFrictionContact3D_test1.c)
  NEW_TEST(PFC3D_test2 PrimalFrictionContact3D_test2.c)
  NEW_TEST(PFC3D_test3 PrimalFrictionContact3D_test3.c)
  NEW_TEST(PFC3D_test4 PrimalFrictionContact3D_test4.c)
  NEW_TEST(PFC3D_test5 PrimalFrictionContact3D_test5.c)


  SET(PFC3D_test6_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(PFC3D_test6 PrimalFrictionContact3D_test6.c)

  NEW_TEST(PFC3D_test7 PrimalFrictionContact3D_test7.c)

  NEW_TEST(PFC3D_test10 PrimalFrictionContact3D_test10.c)
  NEW_TEST(PFC3D_test11 PrimalFrictionContact3D_test11.c)
  NEW_TEST(PFC3D_test12 PrimalFrictionContact3D_test12.c)
  NEW_TEST(PFC3D_test13 PrimalFrictionContact3D_test13.c)
  NEW_TEST(PFC3D_test14 PrimalFrictionContact3D_test14.c)
  NEW_TEST(PFC3D_test14bis PrimalFrictionContact3D_test14bis.c)
  NEW_TEST(PFC3D_test15 PrimalFrictionContact3D_test15.c)

  
  ## Alart Curnier functions
  NEW_TEST(AlartCurnierFunctions_test FrictionContact3D_AlartCurnierFunctions_test.c)
  IF(WITH_FCLIB)
  NEW_TEST(FCLIB_test1 FrictionContact3D_writefclib_local_test.c)
  ENDIF(WITH_FCLIB)


  END_TEST()
  
  BEGIN_TEST(NonSmoothSolvers/GenericMechanical/test)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP0.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP1.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP2.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP3.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP4.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP5.dat)
  SET(test-GMP-REDUCED0_3D_QUARTIC-GMP6_PROPERTIES WILL_FAIL TRUE)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP6.dat)
  
  IF(COMPLETE_LAPACK_LIBRARIES)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP0.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP1.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP2.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP3.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP4.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP5.dat 0 0 0 0 0 1)
#    SET(test-GMP-REDUCED1_3D_QUARTIC-GMP6_PROPERTIES WILL_FAIL TRUE)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP6.dat 0 0 0 0 0 1)
  ENDIF(COMPLETE_LAPACK_LIBRARIES)
  
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP0.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP1.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP2.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP3.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP4.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP5.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP6.dat 0 0 0 0 0 2)
  
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP0.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP1.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP2.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP3.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP6.dat 0 0 0 0 0 3)

  END_TEST()
  #BEGIN_TEST(NonSmoothSolvers/GenericMechanical/test)
  #NEW_TEST(GMP_FAILED GenericMechanical_test1.c)
  #END_TEST()
  
  BEGIN_TEST(utils/test)
  #NEW_TEST(pinvtest testpinv.c)
  NEW_TEST(test_op3x3 test_op3x3.c)
  NEW_TEST(test_timers_interf test_timers_interf.c)
  END_TEST()

ENDIF(WITH_TESTING)

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)

IF(PRINT_VAR)
  MESSAGE(STATUS)
  PRINT_VAR(${PRINT_VAR})
  MESSAGE(STATUS)
ENDIF(PRINT_VAR)
