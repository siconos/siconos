include(tools4tests)

# wrapper are not needed
set(TEST_WRAP)

if(WITH_${COMPONENT}_TESTING)
  
  BEGIN_TEST(src/tools/test)
  if(HAS_LAPACK_DGESVD)
    NEW_TEST(pinvtest testpinv.c)
  endif()
  NEW_TEST(test_op3x3 test_op3x3.c)
  NEW_TEST(test_timers_interf test_timers_interf.c)
  NEW_TEST(test_cblas test_cblas.c)
  NEW_TEST(test_dgesv test_dgesv.c)
  if(HAS_LAPACK_DGESVD)
    NEW_TEST(test_gesvd test_gesvd.c)
  endif()
  if(HAS_LAPACK_DGELS)
    NEW_TEST(test_dgels test_dgels.c)
  endif()
  NEW_TEST(test_dpotrf test_dpotrf.c)
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
  NEW_TEST(SBMTest4 SBM_test4.c)
  NEW_TEST(SBMTest5 SBM_test5.c)
  NEW_TEST(SparseMatrix0 SparseMatrix_test0.c)
  IF(HAS_ONE_LP_SOLVER)
   NEW_TEST(Vertex_extraction vertex_problem.c)
  ENDIF(HAS_ONE_LP_SOLVER)
 END_TEST()

 BEGIN_TEST2(src/LCP/test)

  MACRO(SET_LCP_TEST_AS_FAILED DATASET_LCP_DIAG FAILING_ALGO)
   FOREACH(_DS ${DATASET_LCP_DIAG})
    FOREACH(_SOLVER ${FAILING_ALGO})
     SET(test-LCP_${_SOLVER}-lcp_${_DS}_PROPERTIES WILL_FAIL TRUE)
    ENDFOREACH()
   ENDFOREACH()
  ENDMACRO()

  SET(DATASET_LCP "lcp_mmc.dat;lcp_deudeu.dat;lcp_trivial.dat;lcp_ortiz.dat;lcp_enum_fails.dat")
  LIST(APPEND DATASET_LCP
   "lcp_exp_murty.dat;lcp_exp_murty2.dat;lcp_CPS_1.dat;lcp_CPS_2.dat;lcp_CPS_3.dat;lcp_CPS_4.dat;lcp_CPS_4bis.dat;lcp_CPS_5.dat")
  SET(DATASET_BLOCK_LCP "lcp_deudeu_block.dat")
  # PSOR is not working :(
  SET(SICONOS_LCP_SOLVERS
   "ENUM;LEMKE;CPG;PGS;RPGS;LATIN;LATIN_W;QP;NSQP;AVI_CAOFERRIS;NEWTONMIN;NEWTON_FBLSA;NEWTON_MINFBLSA;BARD;MURTY;PIVOT;PIVOT_LUMOD;PATHSEARCH")
  IF(HAVE_PATHFERRIS)
   LIST(APPEND SICONOS_LCP_SOLVERS "PATH")
  ENDIF()
  IF(HAVE_GAMS_C_API)
   LIST(APPEND SICONOS_LCP_SOLVERS "GAMS")
  ENDIF(HAVE_GAMS_C_API)
  FOREACH(_DS ${DATASET_LCP})
    FOREACH(_SOLVER ${SICONOS_LCP_SOLVERS})
     NEW_LCP_TEST(SICONOS_LCP_${_SOLVER} ${_DS})
    ENDFOREACH()
  ENDFOREACH()
  FOREACH(_DS ${DATASET_BLOCK_LCP})
   FOREACH(_SOLVER ${SICONOS_LCP_SOLVERS})
    NEW_LCP_TEST(SICONOS_LCP_${_SOLVER} ${_DS} 1)
    ENDFOREACH()
  ENDFOREACH()

  # CPG does not work everywhere
  SET(test-LCP_CPG-lcp_deudeu_PROPERTIES WILL_FAIL TRUE)
  SET(test-LCP_CPG-lcp_exp_murty_PROPERTIES WILL_FAIL TRUE)
  SET(test-LCP_CPG-lcp_CPS_2_PROPERTIES WILL_FAIL TRUE)
  SET(test-LCP_CPG-lcp_CPS_4_PROPERTIES WILL_FAIL TRUE)
  SET(test-LCP_CPG-lcp_CPS_4bis_PROPERTIES WILL_FAIL TRUE)

  # problem with Cholesky here
  SET_LCP_TEST_AS_FAILED("exp_murty;exp_murty2" "LATIN;LATIN_W")
  RM_TEST2(SICONOS_LCP_LATIN "lcp_ortiz.dat")
  RM_TEST2(SICONOS_LCP_LATIN_W "lcp_ortiz.dat")

  # QP reformulation does not always work when the matrix is not symmetric
  # Use NSQP
  SET_LCP_TEST_AS_FAILED("exp_murty;exp_murty2;ortiz;enum_fails;CPS_2;CPS_3;CPS_4;CPS_4bis" "QP")

  # NEWTONMIN has no backup descent dir -> problem in DGESV -> GAME OVER !
  SET(test-LCP_NEWTONMIN-lcp_CPS_1_PROPERTIES WILL_FAIL TRUE)
  SET(test-LCP_NEWTONMIN-lcp_CPS_2_PROPERTIES WILL_FAIL TRUE)
  SET(test-LCP_NEWTONMIN-lcp_CPS_5_PROPERTIES WILL_FAIL TRUE)

  # those test cannot be solved with an algorithm that requires non-zero
  # diagonal elements, that is PGS, BARD, MURTY, LATIN and LATIN_W
  SET_LCP_TEST_AS_FAILED("enum_fails;CPS_2;CPS_3;CPS_4;CPS_4bis" "PGS;BARD;MURTY;LATIN;LATIN_W")
  # suprinsingly this works ...
  SET(test-LCP_MURTY-lcp_enum_fails_PROPERTIES WILL_FAIL FALSE)

  # those test cannot be solved with Lemke-based solvers (CPS_3 is for Lemke-Howson)
  SET_LCP_TEST_AS_FAILED("CPS_3" "LEMKE;AVI_CAOFERRIS;PIVOT;PIVOT_LUMOD;PATHSEARCH")

  # PSD matrices and those algo does not seem to be a good idea
  SET_LCP_TEST_AS_FAILED("CPS_2;CPS_3" "NSQP;RPGS")

  # lcp_mmc is of size 26, way too much for enum
  RM_TEST2(SICONOS_LCP_ENUM "lcp_mmc.dat")
  # this LCP was put here to show that enum does not work on every LCP, likely
  # due to numerical problems, but works on some system ...
  RM_TEST2(SICONOS_LCP_ENUM "lcp_enum_fails.dat")

  # TODO backup path when GDESV fails
  SET(test-LCP_NEWTON_FBLSA-lcp_CPS_1_PROPERTIES WILL_FAIL TRUE)

  # special tests
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol.dat)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol_perturbed.dat)
  SET(test-LCP_ENUM-lcp_Pang_isolated_sol_perturbed_PROPERTIES WILL_FAIL TRUE)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_inf_sol_perturbed.dat)

  # TODO refinment of solution
  # NEW_LCP_TEST(SICONOS_LCP_LEMKE lcp_tobenna.dat)
  # NEW_LCP_TEST(SICONOS_LCP_PIVOT lcp_tobenna.dat)
  #  NEW_LCP_TEST(SICONOS_LCP_PIVOT_LUMOD lcp_tobenna.dat)
  # LUMOD is not ready for prime time now
  SET(test-LCP_PIVOT_LUMOD-lcp_mmc_PROPERTIES WILL_FAIL TRUE)
  IF(DEV_MODE)
   SET(test-LCP_PIVOT-lcp_tobenna_PROPERTIES WILL_FAIL FALSE)
   #   SET(test-LCP_PIVOT_LUMOD-lcp_tobenna_PROPERTIES WILL_FAIL FALSE)
  ENDIF(DEV_MODE)

  IF(HAVE_PATHFERRIS)
   NEW_LCP_TEST(SICONOS_LCP_PATH lcp_tobenna.dat)
  ENDIF(HAVE_PATHFERRIS)
  IF(HAVE_GAMS_C_API)
   NEW_LCP_TEST(SICONOS_LCP_GAMS lcp_tobenna.dat)
  ENDIF(HAVE_GAMS_C_API)

  NEW_TEST(LCP_DefaultSolverOptionstest LinearComplementarity_DefaultSolverOptions_test.c)

  END_TEST(LCP/test)

  BEGIN_TEST(src/LinearSystem/test)

  NEW_LS_TEST(SICONOS_LS_0 ls_trivial.dat)
  SET(test-LS_0-ls_inf_sol_perturbed_PROPERTIES WILL_FAIL TRUE)
  NEW_LS_TEST(SICONOS_LS_0 ls_inf_sol_perturbed.dat)
  END_TEST()

  BEGIN_TEST2(src/Relay/test)

  SET(DATA_SET "relay1.dat;relay_2x2.dat;relay_4x4.dat;relay_simple2.dat;step_1x1.dat;step_2x2.dat;step_4x4.dat")
  SET(SICONOS_RELAY_SOLVERS "ENUM;LEMKE;PGS;AVI_CAOFERRIS")
  IF(HAS_ONE_LP_SOLVER)
   LIST(APPEND SICONOS_RELAY_SOLVERS "AVI_CAOFERRIS_TEST")
  ENDIF()

  IF(HAVE_PATHFERRIS)
   LIST(APPEND SICONOS_RELAY_SOLVERS "PATH")
  ENDIF()
  FOREACH(_DS ${DATA_SET})
    FOREACH(_SOLVER ${SICONOS_RELAY_SOLVERS})
      NEW_RELAY_TEST(SICONOS_RELAY_${_SOLVER} ${_DS})
    ENDFOREACH()
  ENDFOREACH()

  # ENUM on an LCP of size 30 is a bad idea ...
  RM_TEST2(RELAY_ENUM "relay1.dat")

  NEW_TEST(Relaytest1 relay_test1.c)
  NEW_TEST(Relaytest2 relay_test2.c)

  IF(HAVE_PATHFERRIS)
    NEW_TEST(Relaytest3 relay_test3.c)
  ENDIF(HAVE_PATHFERRIS)

  NEW_TEST(Relaytest10 relay_test10.c)
  NEW_TEST(Relaytest11 relay_test11.c)
  NEW_TEST(Relaytest12 relay_test12.c)
  NEW_TEST(Relaytest13 relay_test13.c)
  NEW_TEST(Relaytest20 relay_test20.c)
  NEW_TEST(Steptest1 step_test1.c)
  NEW_TEST(Steptest2 step_test2.c)
  NEW_TEST(Steptest3 step_test3.c)
  NEW_TEST(Steptest4 step_test4.c)

  END_TEST()


  BEGIN_TEST(src/MLCP/test)
  IF(HAVE_SYSTIMES_H)
    if(WITH_CXX)
      NEW_TEST(MLCPtest main_mlcp.cpp)
    endif()
  ENDIF(HAVE_SYSTIMES_H)
  NEW_TEST(ReadWrite_MLCPtest MixedLinearComplementarity_ReadWrite_test.c)
  END_TEST()

  BEGIN_TEST(src/MCP/test)
  NEW_TEST(MCPtest MCP_test.c)
  NEW_TEST(MCPtest1 MCP_test1.c)
  END_TEST()

  BEGIN_TEST(src/NCP/test)
  SET(SICONOS_NCP_SOLVERS "NEWTON_FBLSA;NEWTON_MINFBLSA;PATHSEARCH")
  IF(HAVE_PATHFERRIS)
    LIST(APPEND SICONOS_NCP_SOLVERS "PATH")
  ENDIF(HAVE_PATHFERRIS)

  SET(SICONOS_NCP_TEST_PROBLEMS "NCP_ZI1")
  IF(DEV_MODE)
   LIST(APPEND SICONOS_NCP_TEST_PROBLEMS "NCP_ZIT1")
  ENDIF(DEV_MODE)

  FOREACH(_PB ${SICONOS_NCP_TEST_PROBLEMS})
   FOREACH(_SOLVER ${SICONOS_NCP_SOLVERS})
    NEW_NCP_TEST(${_PB} SICONOS_NCP_${_SOLVER})
   ENDFOREACH()
  ENDFOREACH()

  IF(NOT DEV_MODE)
   SET(NCP_NEWTON_FBLSA-NCP_ZI1_PROPERTIES WILL_FAIL TRUE)
  ENDIF(NOT DEV_MODE)

  END_TEST() # NCP

  BEGIN_TEST(src/FrictionContact/test)
  
  #NEW_TEST(FrictionContact_Problemtest main_FC3D.c)
  NEW_TEST(FC3D_DefaultSolverOptionstest fc3d_DefaultSolverOptions_test.c)
  
  NEW_TEST(FC3D_sparse_test fc3d_sparse_test.c)

  # (see FrictionContact/test/README for short details)
  NEW_TEST(test_fc3d_1 fc3d_test1.c)
  NEW_TEST(test_fc3d_2 fc3d_test2.c)
  NEW_TEST(test_fc3d_3 fc3d_test3.c)
  NEW_TEST(test_fc3d_4 fc3d_test4.c)
  SET(test_fc3d_5_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(test_fc3d_5 fc3d_test5.c)
  
  IF(HAVE_PATHFERRIS)
    NEW_TEST(test_fc3d_6 fc3d_test6.c)
    SET(test_fc3d_6_PROPERTIES WILL_FAIL TRUE)
  ENDIF(HAVE_PATHFERRIS)
  
  SET(test_fc3d_7_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(test_fc3d_7 fc3d_test7.c)
  
  
  NEW_TEST(test_fc3d_9 fc3d_test9.c)



  NEW_TEST(test_fc3d_10 fc3d_test10.c)
  NEW_TEST(test_fc3d_11 fc3d_test11.c)
  NEW_TEST(test_fc3d_12 fc3d_test12.c)
  NEW_TEST(test_fc3d_13 fc3d_test13.c)
  NEW_TEST(test_fc3d_14 fc3d_test14.c)

  SET(NSGS_TOL 1e-6)
  SET(NSGS_NB_IT 10000)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Example1_Fc3D.dat 1e-5 ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_NSN_AC 1e-18 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover1039.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover1040.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover1041.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover3865.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4144.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4396.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4493.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4516.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4609.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4613.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover4622.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover9770.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover11035.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSGS Rover11211.dat ${NSGS_TOL} ${NSGS_NB_IT} 
    SICONOS_FRICTION_3D_ONECONTACT_QUARTIC 1e-6 10)

  SET(test_fc3d_20_PROPERTIES WILL_FAIL TRUE)

  NEW_TEST(test_fc3d_20 fc3d_test20.c) 
  NEW_TEST(test_fc3d_21 fc3d_test21.c)
  NEW_TEST(test_fc3d_22 fc3d_test22.c)
  NEW_TEST(test_fc3d_23 fc3d_test23.c)
  
  NEW_TEST(test_fc3d_30 fc3d_test30.c)
  NEW_TEST(test_fc3d_31 fc3d_test31.c)
  NEW_TEST(test_fc3d_32 fc3d_test32.c)
  
  NEW_TEST(test_fc3d_33 fc3d_test33.c) # DSFP converges with specific rho

  NEW_TEST(test_fc3d_3400 fc3d_test3400.c) # EG 
  NEW_TEST(test_fc3d_35 fc3d_test35.c)     # EG

  NEW_TEST(test_fc3d_36 fc3d_test36.c) #TFP with NSGS and projection on cylinder

  #SET(test_fc3d_37_PROPERTIES WILL_FAIL TRUE)
  #NEW_TEST(test_fc3d_37 fc3d_test37.c) #TFP with ProjectedGradientOnCylinder is not working ...

  #NEW_TEST(test_fc3d_38 fc3d_test38.c) # HP is not converging

  NEW_TEST(test_fc3d_40 fc3d_test40.c) # VI_EG 
  NEW_TEST(test_fc3d_41 fc3d_test41.c) # VI_FPP
  NEW_TEST(test_fc3d_42 fc3d_test42.c) # VI_FPP
  NEW_TEST(test_fc3d_43 fc3d_test43.c) # DSFP converges with specific rho
  NEW_TEST(test_fc3d_44 fc3d_test44.c) # VI_EG
  NEW_TEST(test_fc3d_45 fc3d_test45.c) # VI_EG
  
  NEW_TEST(test_fc3d_46 fc3d_test46.c) # FPP
  NEW_TEST(test_fc3d_47 fc3d_test47.c) # EG

  SET(test_fc3d_50_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(test_fc3d_50 fc3d_test50.c)
  

  NEW_TEST(test_fc3d_60 fc3d_test60.c)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Example1_Fc3D_SBM.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Example1_Fc3D_SBM.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Example1_Fc3D_SBM.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC OneObject-i100000-499.hdf5.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB OneObject-i100000-499.hdf5.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM OneObject-i100000-499.hdf5.dat)
  
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Capsules-i100-1090.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Capsules-i100-1090.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Capsules-i100-1090.dat)
  
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Capsules-i101-404.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Capsules-i101-404.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Capsules-i101-404.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Capsules-i103-990.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Capsules-i103-990.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Capsules-i103-990.dat)
  
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Capsules-i122-1617.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Capsules-i122-1617.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Capsules-i122-1617.dat)
  
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Capsules-i100-889.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Capsules-i100-889.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Capsules-i100-889.dat)
  
  NEW_TEST(test_fc3d_70 fc3d_test70.c) # DSFP
  NEW_TEST(test_fc3d_71 fc3d_test71.c) # DSFP
  NEW_TEST(test_fc3d_72 fc3d_test72.c)	
  NEW_TEST(test_fc3d_73 fc3d_test73.c)
  NEW_TEST(test_fc3d_74 fc3d_test74.c)
  
  NEW_TEST(test_fc3d_80 fc3d_test80.c) # Proximal
  NEW_TEST(test_fc3d_81 fc3d_test81.c) # Proximal
  NEW_TEST(test_fc3d_83 fc3d_test83.c) # Proximal

  NEW_TEST(test_fc3d_82 fc3d_test82.c) # Proximal
  SET(test_fc3d_90_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(test_fc3d_90 fc3d_test90.c) # TFP with cycling

  NEW_TEST(FC3DNewFromFortranData fc3d_newFromFortranData.c)
  NEW_TEST(FC3DLmgcDriver1 fc3d_LmgcDriver_test1.c)
  NEW_TEST(FC3DLmgcDriver2 fc3d_LmgcDriver_test2.c)
  NEW_TEST(FC3DLmgcDriver3 fc3d_LmgcDriver_test3.c)

  NEW_TEST(FC3DLmgcDriver4 fc3d_LmgcDriver_test4.c)

  NEW_TEST(FC3DLmgcDriver5 fc3d_LmgcDriver_test5.c)

  NEW_TEST(test_fc3d_100 fc3d_test100.c)# TFP
  NEW_TEST(test_fc3d_110 fc3d_test110.c)# TFP

  NEW_TEST(test_fc3d_120 fc3d_test120.c) # TFP
  NEW_TEST(test_fc3d_121 fc3d_test121.c) # DSFP

  NEW_TEST(test_fc3d_122 fc3d_test122.c) # ACLMFP
  NEW_TEST(test_fc3d_123 fc3d_test123.c) # SOCLCP
  NEW_TEST(test_fc3d_124 fc3d_test124.c) # NSGS

  NEW_TEST(test_fc3d_125 fc3d_test125.c) # TFP with other strategy for internal solver
  NEW_TEST(test_fc3d_126 fc3d_test126.c) # ACLMFPwith other strategy for internal solver
 

  NEW_TEST(test_fc3d_130 fc3d_test130.c)

  NEW_TEST(test_fc3d_200 fc3d_test200.c)
  NEW_TEST(test_fc3d_210 fc3d_test210.c)

  NEW_TEST(test_fc3d_300 fc3d_test300.c)


  # Global Alart Curnier + Rover ok.
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Example1_Fc3D.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Example1_Fc3D.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Example1_Fc3D.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover1039.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover1039.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover1039.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover1040.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover1040.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover1040.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover1041.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover1041.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover1041.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover3865.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover3865.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover3865.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover4144.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover4144.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover4144.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover4396.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover4396.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover4396.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover4493.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover4493.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover4493.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover4516.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover4516.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover4516.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover4609.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover4609.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover4609.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover4613.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover4613.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover4613.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover4622.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover4622.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover4622.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover9770.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover9770.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover9770.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover11035.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover11035.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover11035.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC Rover11211.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB Rover11211.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM Rover11211.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC FrictionContact3D_1c.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC FrictionContact3D_RR_1c.dat)

  # Newton Euler 3D Spheres : some failures with  Alart Curnier
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC NESpheres_10_1.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB NESpheres_10_1.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM NESpheres_10_1.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_AC NESpheres_30_1.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_FB NESpheres_30_1.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_NSN_NM NESpheres_30_1.dat)

  ## test from the Rover Example
  NEW_TEST(test_fc3d_501 fc3d_Rover2.c) ## test ok without LineSearch
  NEW_TEST(test_fc3d_503 fc3d_Rover3.c) ## JeanMoreau ok, AC:no


  ## test 2D dense on two differents files
  NEW_TEST(test_fc2d_1 fc2d_test1.c)
  NEW_TEST(test_fc2d_2 fc2d_test2.c)
  NEW_TEST(test_fc2d_3 fc2d_test3.c)
  
  
  NEW_TEST(test_fc2d_10 fc2d_test10.c)
  NEW_TEST(test_fc2d_11 fc2d_test11.c)
  NEW_TEST(test_fc2d_12 fc2d_test12.c)
  

  ## test 2D sparse on 4 differents files
  NEW_TEST(test_fc2d_20 fc2d_test20.c)
  NEW_TEST(test_fc2d_21 fc2d_test21.c)
  NEW_TEST(test_fc2d_22 fc2d_test22.c)
  NEW_TEST(test_fc2d_23 fc2d_test23.c)

  ## test 2D dense with Lemke NSGS failed on it
  NEW_TEST(test_fc2d_30 fc2d_test30.c)
  NEW_TEST(test_fc2d_31 fc2d_test31.c)
  #SET(test_fc2d_32 WILL_FAIL TRUE)
  #NEW_TEST(test_fc2d_32 fc2d_test32.c)




  ## test 2D dense with Enum Lemke failed on 41 !!
  NEW_TEST(test_fc2d_40 fc2d_test40.c)
  NEW_TEST(test_fc2d_41 fc2d_test41.c)


  NEW_TEST(GFC3D_test1 gfc3d_test1.c)
  NEW_TEST(GFC3D_test2 gfc3d_test2.c)
  NEW_TEST(GFC3D_test3 gfc3d_test3.c)
  NEW_TEST(GFC3D_test4 gfc3d_test4.c)
  NEW_TEST(GFC3D_test5 gfc3d_test5.c)


  SET(GFC3D_test6_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(GFC3D_test6 gfc3d_test6.c)
  
  NEW_TEST(GFC3D_test7 gfc3d_test7.c)
  
  NEW_TEST(GFC3D_test10 gfc3d_test10.c)
  NEW_TEST(GFC3D_test11 gfc3d_test11.c)
  NEW_TEST(GFC3D_test12 gfc3d_test12.c)
  NEW_TEST(GFC3D_test13 gfc3d_test13.c)
  NEW_TEST(GFC3D_test14 gfc3d_test14.c)
  NEW_TEST(GFC3D_test14bis gfc3d_test14bis.c)
  NEW_TEST(GFC3D_test15 gfc3d_test15.c)
  NEW_TEST(GFC3D_test16 gfc3d_test16.c)

  
  ## Alart Curnier functions
  NEW_TEST(AlartCurnierFunctions_test fc3d_AlartCurnierFunctions_test.c)
  IF(WITH_FCLIB)
  NEW_TEST(FCLIB_test1 fc3d_writefclib_local_test.c)
  NEW_TEST(FCLIB_GFC3D_test1 gfc3d_fclib_cubeH8.c)
  NEW_TEST(FCLIB_GFC3D_test2 gfc3d_test20.c)
  ENDIF(WITH_FCLIB)
  SET(FC3D_DATA_SET
   "Capsules-i100-1090.dat;Capsules-i100-889.dat;Capsules-i101-404.dat;Capsules-i103-990.dat;Capsules-i122-1617.dat;Example1_Fc3D.dat;Example1_Fc3D_SBM.dat;FrictionContact3D_1c.dat;FrictionContact3D_RR_1c.dat;NESpheres_10_1.dat;NESpheres_30_1.dat;OneObject-i100000-499.hdf5.dat;Rover1039.dat;Rover1040.dat;Rover1041.dat;Rover11035.dat;Rover11211.dat;Rover3865.dat;Rover4144.dat;Rover4396.dat;Rover4493.dat;Rover4516.dat;Rover4609.dat;Rover4613.dat;Rover4622.dat;Rover9770.dat;KaplasTower-i1061-4.hdf5.dat;Confeti-ex13-Fc3D-SBM.dat;BoxesStack1-i100000-32.hdf5.dat;FrictionContactProblem00237.dat;FrictionContactProblem00727.dat")

  IF(HAVE_GAMS_C_API)
    FOREACH(_DAT ${FC3D_DATA_SET})
     NEW_FC_TEST(SICONOS_FRICTION_3D_GAMS_PATH ${_DAT})
     NEW_FC_TEST(SICONOS_FRICTION_3D_GAMS_LCP_PATH ${_DAT})
     IF(HAVE_GAMS_PATHVI)
       NEW_FC_TEST(SICONOS_FRICTION_3D_GAMS_PATHVI ${_DAT})
       NEW_FC_TEST(SICONOS_FRICTION_3D_GAMS_LCP_PATHVI ${_DAT})
      ENDIF(HAVE_GAMS_PATHVI)
    ENDFOREACH(_DAT ${FC3D_DATA_SET})
  ENDIF(HAVE_GAMS_C_API)

  END_TEST()





  
  BEGIN_TEST(src/GenericMechanical/test)

  # Warning. Only test GMP4 GMP5 have fc3d problem inside

  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP0.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP1.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP2.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP3.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP4.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP5.dat)
  SET(test-GMP-REDUCED0_3D_ONECONTACT_QUARTIC-GMP6_PROPERTIES WILL_FAIL TRUE)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP6.dat)

  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP0.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP1.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP2.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP3.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP4.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP5.dat)

  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC GMP4.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC GMP5.dat)

  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration GMP4.dat)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration GMP5.dat)
  
  IF(HAS_LAPACK_DGESVD)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP0.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP1.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP2.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP3.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP4.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP5.dat 0 0 0 0 0 1)
#    SET(test-GMP-REDUCED1_3D_QUARTIC-GMP6_PROPERTIES WILL_FAIL TRUE)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP6.dat 0 0 0 0 0 1)
  endif()
  
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP0.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP1.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP2.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP3.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP4.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP5.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP6.dat 0 0 0 0 0 2)
  
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP0.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP1.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP2.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP3.dat 0 0 0 0 0 3)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC GMP6.dat 0 0 0 0 0 3)

  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP4.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP5.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP GMP6.dat 0 0 0 0 0 2)

  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC GMP4.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC GMP5.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC GMP6.dat 0 0 0 0 0 2)

  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration GMP4.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration GMP5.dat 0 0 0 0 0 2)
  NEW_GMP_TEST(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration GMP6.dat 0 0 0 0 0 2)





  END_TEST()
  #BEGIN_TEST(src/GenericMechanical/test)
  #NEW_TEST(GMP_FAILED GenericMechanical_test1.c)
  #END_TEST()
  BEGIN_TEST(src/VI/test)
  NEW_TEST(VI_test0 VI_test.c)
  NEW_TEST(VI_test1 VI_test1.c)
  NEW_TEST(VI_test2 VI_test2.c)
  NEW_TEST(VI_test3 VI_test3.c)
  NEW_TEST(VI_test5 VI_test5.c)
  NEW_TEST(VI_testFC3D1 VI_testFC3D1.c)
  NEW_TEST(VI_testFC3D2 VI_testFC3D2.c)
  NEW_TEST(VI_testFC3D3 VI_testFC3D3.c)
  SET(VI_testFC3D3_PROPERTIES WILL_FAIL TRUE)


  SET(SICONOS_VI_SOLVERS "BOX_QI;BOX_AVI_LSA")
  IF(HAVE_PATHFERRIS)
   LIST(APPEND SICONOS_VI_SOLVERS "BOX_PATH")
  ENDIF(HAVE_PATHFERRIS)

  IF(DEV_MODE)
   SET(SICONOS_VI_TEST_PROBLEMS "VI_ZI1;VI_ZIT1")
  ENDIF(DEV_MODE)

  FOREACH(_PB ${SICONOS_VI_TEST_PROBLEMS})
   FOREACH(_SOLVER ${SICONOS_VI_SOLVERS})
    NEW_NCP_TEST(${_PB} SICONOS_VI_${_SOLVER})
   ENDFOREACH()
  ENDFOREACH()
 END_TEST()

  BEGIN_TEST(src/AVI/test)

  IF(HAS_ONE_LP_SOLVER)
   NEW_TEST(AVI_twisting implicit_twisting.c)
  ENDIF(HAS_ONE_LP_SOLVER)

  END_TEST(AVI/test)


  BEGIN_TEST(src/SOCP/test)
    NEW_TEST(SOCLCP_test1 soclcp_test1.c)
    NEW_TEST(SOCLCP_test2 soclcp_test2.c)
    NEW_TEST(SOCLCP_test3 soclcp_test3.c)
    # timeout on all machines, see
    # http://cdash-bipop.inrialpes.fr/testSummary.php?project=1&name=SOCLCP_test4&date=2015-09-03
    # Feel free to remove this once it is fixed --xhub
    #NEW_TEST(SOCLCP_test4 soclcp_test4.c)
    #NEW_TEST(SOCLCP_test5 soclcp_test5.c)
    NEW_TEST(SOCLCP_fc3d_to_soclcp  fc3d_to_soclcp.c)
  END_TEST(SOCP/test)

  add_library(numerics-test SHARED ${TEST_UTILS_SOURCES})
  set_target_properties("numerics-test" PROPERTIES
   LINKER_LANGUAGE "C")
  include(WindowsLibrarySetup)
  windows_library_extra_setup("numerics-test" "numerics-test")
  target_link_libraries(numerics-test ${PRIVATE} ${COMPONENT})
  target_link_libraries(numerics-test ${PRIVATE} ${${COMPONENT}_LINK_LIBRARIES})
endif()
