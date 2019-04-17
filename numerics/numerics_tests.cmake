include(tools4tests)

# wrapper are not needed
set(TEST_WRAP)

if(WITH_${COMPONENT}_TESTING)

  BEGIN_TEST(src/tools/test)

  NEW_TEST(tools_test_op3x3 test_op3x3.c)
  
  NEW_TEST(tools_test_timers_interf test_timers_interf.c)

  NEW_TEST(tools_test_blas_lapack test_blas_lapack.c)
  if(HAS_LAPACK_DGESVD)
    NEW_TEST(tools_test_pinv test_pinv.c)
  endif()
  NEW_TEST(tools_projection test_projection.c)
  NEW_TEST(tools_test_NumericsArrays NumericsArrays.c)

  #  tests for NumericsMatrix
  NEW_TEST(tools_test_NumericsMatrix NM_test.c)

  # MUMPS interface tests
  if(WITH_MUMPS)
    NEW_TEST(tools_test_MUMPS NM_MUMPS_test.c)
  endif()
  
  # Specific tests for SBM matrices
  NEW_TEST(tools_test_SBM SBM_test.c)
  NEW_TEST(tools_test_SBCM_to_SBM SBCM_to_SBM.c)
  
  # Specific tests for sparse matrices
  NEW_TEST(tools_test_SparseMatrix SparseMatrix_test.c)
  
  IF(HAS_ONE_LP_SOLVER)
    NEW_TEST(tools_test_Vertex_extraction vertex_problem.c)
  ENDIF(HAS_ONE_LP_SOLVER)
  END_TEST()

  BEGIN_TEST2(src/LCP/test)

  FILE(GLOB_RECURSE _DATA_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_D}
    data_collection*.c
    test_*.c)
  
  FOREACH(_F ${_DATA_FILES})
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_F} ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_F} COPYONLY)
  ENDFOREACH(_F ${_DATA_FILES})

  NEW_TEST(lcp_test_DefaultSolverOptions LinearComplementarity_DefaultSolverOptions_test.c)

  NEW_LCP_TEST_COLLECTION(TEST_LCP_COLLECTION_1)
  NEW_LCP_TEST_COLLECTION(TEST_LCP_COLLECTION_2)
  NEW_LCP_TEST_COLLECTION(TEST_LCP_COLLECTION_3)
  NEW_LCP_TEST_COLLECTION(TEST_LCP_COLLECTION_4)


  
  # MACRO(SET_LCP_TEST_AS_FAILED DATASET_LCP_DIAG FAILING_ALGO)
  #   FOREACH(_DS ${DATASET_LCP_DIAG})
  #     FOREACH(_SOLVER ${FAILING_ALGO})
  # 	SET(test-LCP_${_SOLVER}-lcp_${_DS}_PROPERTIES WILL_FAIL TRUE)
  #     ENDFOREACH()
  #   ENDFOREACH()
  # ENDMACRO()

  # SET(DATASET_LCP "lcp_mmc.dat;lcp_deudeu.dat;lcp_trivial.dat;lcp_ortiz.dat;lcp_enum_fails.dat")
  # LIST(APPEND DATASET_LCP
  #   "lcp_exp_murty.dat;lcp_exp_murty2.dat;lcp_CPS_1.dat;lcp_CPS_2.dat;lcp_CPS_3.dat;lcp_CPS_4.dat;lcp_CPS_4bis.dat;lcp_CPS_5.dat")
  # SET(DATASET_BLOCK_LCP "lcp_deudeu_block.dat")
  # # PSOR is not working :(

  # SET(SICONOS_LCP_SOLVERS
  #   "ENUM;LEMKE;CPG;PGS;RPGS;LATIN;LATIN_W;AVI_CAOFERRIS;NEWTONMIN;NEWTON_FBLSA;NEWTON_MINFBLSA;BARD;MURTY;PIVOT;PIVOT_LUMOD;PATHSEARCH;CONVEXQP_PG")
  # if(HAS_FORTRAN AND HAVE_QL0001)
  #   LIST(APPEND SICONOS_LCP_SOLVERS "QP;NSQP;")
  # endif()
  # IF(HAVE_PATHFERRIS)
  #   LIST(APPEND SICONOS_LCP_SOLVERS "PATH")
  # ENDIF()
  # IF(HAVE_GAMS_C_API)
  #   LIST(APPEND SICONOS_LCP_SOLVERS "GAMS")
  # ENDIF(HAVE_GAMS_C_API)
  # FOREACH(_DS ${DATASET_LCP})
  #   FOREACH(_SOLVER ${SICONOS_LCP_SOLVERS})
  #     NEW_LCP_TEST(SICONOS_LCP_${_SOLVER} ${_DS})
  #   ENDFOREACH()
  # ENDFOREACH()
  # FOREACH(_DS ${DATASET_BLOCK_LCP})
  #   FOREACH(_SOLVER ${SICONOS_LCP_SOLVERS})
  #     NEW_LCP_TEST(SICONOS_LCP_${_SOLVER} ${_DS} 1)
  #   ENDFOREACH()
  # ENDFOREACH()

  # # CPG does not work everywhere
  # SET(test-LCP_CPG-lcp_exp_murty_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_CPG-lcp_CPS_2_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_CPG-lcp_CPS_4_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_CPG-lcp_CPS_4bis_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_CPG-lcp_enum_fails_PROPERTIES WILL_FAIL TRUE)

  # # problem with Cholesky here
  # SET_LCP_TEST_AS_FAILED("exp_murty;exp_murty2" "LATIN;LATIN_W")
  # RM_TEST2(SICONOS_LCP_LATIN "lcp_ortiz.dat")
  # RM_TEST2(SICONOS_LCP_LATIN_W "lcp_ortiz.dat")

  # # QP reformulation does not always work when the matrix is not symmetric
  # # Use NSQP
  # SET_LCP_TEST_AS_FAILED("exp_murty;exp_murty2;ortiz;enum_fails;CPS_2;CPS_3;CPS_4;CPS_4bis" "QP")
  # SET_LCP_TEST_AS_FAILED("exp_murty;exp_murty2;" "CONVEXQP_PG")

  # # NEWTONMIN has no backup descent dir -> problem in DGESV -> GAME OVER !
  # SET(test-LCP_NEWTONMIN-lcp_CPS_1_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_NEWTONMIN-lcp_CPS_2_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_NEWTONMIN-lcp_CPS_5_PROPERTIES WILL_FAIL TRUE)



  # # NaN showing up in DGESV -> NEWTONMIN looks really buggy
  # SET(test-LCP_NEWTONMIN-lcp_CPS_4_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_NEWTONMIN-lcp_CPS_4bis_PROPERTIES WILL_FAIL TRUE)
  # SET(test-LCP_NEWTONMIN-lcp_enum_fails_PROPERTIES WILL_FAIL TRUE)

  # IF(NOT WITH_UNSTABLE_TEST)
  #   RM_TEST2(SICONOS_LCP_NEWTONMIN "lcp_mmc.dat")
  # ENDIF()


  # # those test cannot be solved with an algorithm that requires non-zero
  # # diagonal elements, that is PGS, BARD, MURTY, LATIN and LATIN_W
  # SET_LCP_TEST_AS_FAILED("enum_fails;CPS_2;CPS_3;CPS_4;CPS_4bis" "PGS;BARD;MURTY;LATIN;LATIN_W;CONVEXQP_PG")
  # # suprinsingly this works ...
  # SET(test-LCP_MURTY-lcp_enum_fails_PROPERTIES WILL_FAIL FALSE)

  # # those test cannot be solved with Lemke-based solvers (CPS_3 is for Lemke-Howson)
  # SET_LCP_TEST_AS_FAILED("CPS_3" "LEMKE;AVI_CAOFERRIS;PIVOT;PIVOT_LUMOD;PATHSEARCH")

  # # PSD matrices and those algo does not seem to be a good idea
  # SET_LCP_TEST_AS_FAILED("CPS_2;CPS_3" "NSQP;RPGS")

  # # lcp_mmc is of size 26, way too much for enum
  # RM_TEST2(SICONOS_LCP_ENUM "lcp_mmc.dat")
  # # this LCP was put here to show that enum does not work on every LCP, likely
  # # due to numerical problems, but works on some system ...
  # RM_TEST2(SICONOS_LCP_ENUM "lcp_enum_fails.dat")

  # # TODO backup path when GDESV fails
  # SET(test-LCP_NEWTON_FBLSA-lcp_CPS_1_PROPERTIES WILL_FAIL TRUE)

  # special tests
  # NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol.dat)
  # NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol_perturbed.dat)
  # SET(test-LCP_ENUM-lcp_Pang_isolated_sol_perturbed_PROPERTIES WILL_FAIL TRUE)
  # NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_inf_sol_perturbed.dat)

  # TODO refinment of solution
  # NEW_LCP_TEST(SICONOS_LCP_LEMKE lcp_tobenna.dat)
  # NEW_LCP_TEST(SICONOS_LCP_PIVOT lcp_tobenna.dat)
  #  NEW_LCP_TEST(SICONOS_LCP_PIVOT_LUMOD lcp_tobenna.dat)
  # LUMOD is not ready for prime time now
  # SET(test-LCP_PIVOT_LUMOD-lcp_mmc_PROPERTIES WILL_FAIL TRUE)
  # IF(DEV_MODE)
  #   SET(test-LCP_PIVOT-lcp_tobenna_PROPERTIES WILL_FAIL FALSE)
  #   #   SET(test-LCP_PIVOT_LUMOD-lcp_tobenna_PROPERTIES WILL_FAIL FALSE)
  # ENDIF(DEV_MODE)

  # IF(HAVE_PATHFERRIS)
  #   NEW_LCP_TEST(SICONOS_LCP_PATH lcp_tobenna.dat)
  # ENDIF(HAVE_PATHFERRIS)
  # IF(HAVE_GAMS_C_API)
  #   NEW_LCP_TEST(SICONOS_LCP_GAMS lcp_tobenna.dat)
  # ENDIF(HAVE_GAMS_C_API)


  END_TEST(LCP/test)
  
  BEGIN_TEST2(src/Relay/test)
  FILE(GLOB_RECURSE _DATA_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_D}
    data_collection*.c
    test_*.c)  
  FOREACH(_F ${_DATA_FILES})
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_F} ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_F} COPYONLY)
  ENDFOREACH(_F ${_DATA_FILES})

  NEW_RELAY_TEST_COLLECTION(TEST_RELAY_COLLECTION_1)

  NEW_TEST(relay_test_20 relay_test20.c)
  NEW_TEST(step_test_1 step_test1.c)
  NEW_TEST(step_test_2 step_test2.c)
  NEW_TEST(step_test_3 step_test3.c)
  NEW_TEST(step_test_4 step_test4.c)

  END_TEST()

  BEGIN_TEST(src/MLCP/test)
  IF(HAVE_SYSTIMES_H AND WITH_CXX)
    NEW_TEST(MLCP_test_0 main_mlcp.cpp)
  ENDIF(HAVE_SYSTIMES_H AND WITH_CXX)
  NEW_TEST(MLCP_test_read_write MixedLinearComplementarity_ReadWrite_test.c)
  END_TEST()

  BEGIN_TEST(src/MCP/test)
  NEW_TEST(MCP_test_0 MCP_test.c)
  NEW_TEST(MCP_test_1 MCP_test1.c)
  END_TEST()

  BEGIN_TEST(src/NCP/test)
  SET(SICONOS_NCP_SOLVERS "NEWTON_FBLSA;NEWTON_MINFBLSA;PATHSEARCH")
  IF(HAVE_PATHFERRIS)
    LIST(APPEND SICONOS_NCP_SOLVERS "PATH")
  ENDIF(HAVE_PATHFERRIS)


  IF(WITH_UNSTABLE_TEST)
    SET(SICONOS_NCP_TEST_PROBLEMS "NCP_ZI1")
    IF(DEV_MODE)
      LIST(APPEND SICONOS_NCP_TEST_PROBLEMS "NCP_ZIT1")
    ENDIF(DEV_MODE)
  ENDIF()

  FOREACH(_PB ${SICONOS_NCP_TEST_PROBLEMS})
    FOREACH(_SOLVER ${SICONOS_NCP_SOLVERS})
      NEW_NCP_TEST(${_PB} SICONOS_NCP_${_SOLVER})
    ENDFOREACH()
  ENDFOREACH()

  # Oliverie

  IF(NOT DEV_MODE)
    SET(NCP_NEWTON_FBLSA-NCP_ZI1_PROPERTIES WILL_FAIL TRUE)
  ENDIF(NOT DEV_MODE)
  SET(NCP_NEWTON_FBLSA-NCP_ZI1_TIMEOUT 60)

  END_TEST() # NCP

  BEGIN_TEST(src/FrictionContact/test)
  #===========================================
  # 3D Friction Contact tests
  #===========================================
  # (see FrictionContact/test/README for short details)
  # --> Must be uptodated!
  # Set name of input file used to generate c-files for tests
  FILE(GLOB_RECURSE _DATA_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_D}
    data_collection*.c
    test_*.c)
  
  FOREACH(_F ${_DATA_FILES})
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_F} ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_F} COPYONLY)
  ENDFOREACH(_F ${_DATA_FILES})
  
  NEW_FC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_1)
  NEW_FC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_2)
  NEW_FC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_3)
  NEW_FC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_4)
  NEW_FC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_5)

  
  NEW_FC_3D_TEST_COLLECTION(TEST_ADMM_COLLECTION_1)
  NEW_FC_3D_TEST_COLLECTION(TEST_ADMM_COLLECTION_2)
  
  NEW_FC_3D_TEST_COLLECTION(TEST_NSN_COLLECTION_1)
  NEW_FC_3D_TEST_COLLECTION(TEST_NSN_COLLECTION_2)
  NEW_FC_3D_TEST_COLLECTION(TEST_NSN_COLLECTION_3)

  NEW_FC_3D_TEST_COLLECTION(TEST_VI_BASED_COLLECTION_1)
  NEW_FC_3D_TEST_COLLECTION(TEST_FP_COLLECTION_0)
  NEW_FC_3D_TEST_COLLECTION(TEST_FP_COLLECTION_1)
  NEW_FC_3D_TEST_COLLECTION(TEST_FP_COLLECTION_2)
  NEW_FC_3D_TEST_COLLECTION(TEST_PROX_COLLECTION_1)
  
  NEW_FC_3D_TEST_COLLECTION(TEST_QUARTIC_COLLECTION_1)


  #  foreach(_DAT ${FC3D_DATA_SET})
  #    #MESSAGE(STATUS "Setting test for ${_DAT}")
  #    # --- GAMS Solvers ---
  #    if(HAVE_GAMS_C_API)
  #      NEW_FC_3D_TEST(${_DAT} SICONOS_FRICTION_3D_GAMS_PATH)
  #      NEW_FC_3D_TEST(${_DAT} SICONOS_FRICTION_3D_GAMS_LCP_PATH)
  #      if(HAVE_GAMS_PATHVI)
  #  	NEW_FC_3D_TEST(${_DAT} SICONOS_FRICTION_3D_GAMS_PATHVI)
  #  	NEW_FC_3D_TEST(${_DAT} SICONOS_FRICTION_3D_GAMS_LCP_PATHVI)
  #      endif()
  #    endif()
  #  endforeach()


  #  IF(WITH_UMFPACK)
  #    SET(fc3d__ADMM_Tol_1e-5_Max_10000_inTol_0_inMax_0_IPARAM_SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT_Confeti-ex13-Fc3D-SBM_PROPERTIES WILL_FAIL TRUE)
  #    SET(fc3d__ADMM_Tol_1e-5_Max_10000_inTol_0_inMax_0_IPARAM_SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY_Confeti-ex13-Fc3D-SBM_PROPERTIES WILL_FAIL TRUE)
  #    SET(fc3d__NSN_AC_Tol_5e-2_Max_1000_inTol_0_inMax_0_KaplasTower-i1061-4.hdf5_PROPERTIES WILL_FAIL TRUE)
  #  ENDIF()

  IF (WITH_UNSTABLE_TEST)
    NEW_FC_3D_TEST(BoxesStack1-i100000-32.hdf5.dat SICONOS_FRICTION_3D_NSN_AC 1e-5 5000
      0 0 0
      IPARAM 1 1)
    NEW_FC_3D_TEST(BoxesStack1-i100000-32.hdf5.dat SICONOS_FRICTION_3D_NSN_AC_TEST 1e-5 500
      0 0 0
      IPARAM 1 1
      DPARAM SICONOS_DPARAM_LSA_ALPHA_MIN 0.0) # alpha_min needs to be equal to zero for convergence
    NEW_FC_3D_TEST(BoxesStack1-i100000-32.hdf5.dat SICONOS_FRICTION_3D_NSN_AC_TEST 1e-3 1000
      0 0 0
      DPARAM SICONOS_DPARAM_LSA_ALPHA_MIN 0.0) # alpha_min needs to be equal to zero for convergence
  ENDIF()

  # --- LMGC driver ---

  NEW_TEST(FC3DNewFromFortranData fc3d_newFromFortranData.c)
  NEW_TEST(FC3DLmgcDriver1 fc3d_LmgcDriver_test1.c)
  NEW_TEST(FC3DLmgcDriver2 fc3d_LmgcDriver_test2.c)
  NEW_TEST(FC3DLmgcDriver3 fc3d_LmgcDriver_test3.c)

  NEW_TEST(FC3DLmgcDriver4 fc3d_LmgcDriver_test4.c)

  NEW_TEST(FC3DLmgcDriver5 fc3d_LmgcDriver_test5.c)

  
  # # --- Quartic ---
  # NEW_FC_3D_TEST(FrictionContact3D_1c.dat SICONOS_FRICTION_3D_ONECONTACT_QUARTIC)
  # NEW_FC_3D_TEST(FrictionContact3D_RR_1c.dat SICONOS_FRICTION_3D_ONECONTACT_QUARTIC)

  # ---------------------------------------------------
  # --- Global friction contact problem formulation ---
  # ---------------------------------------------------
  
  NEW_GFC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_1)
  NEW_GFC_3D_TEST_COLLECTION(TEST_WR_COLLECTION_1)
  NEW_GFC_3D_TEST_COLLECTION(TEST_NSN_COLLECTION_1)


  IF (WITH_UNSTABLE_TEST)
    NEW_GFC_3D_TEST(GFC3D_TwoRods1.dat SICONOS_GLOBAL_FRICTION_3D_NSN_AC 0 0
      0 0 0
      WILL_FAIL) # pass with mumps only
  ENDIF()
  # Alart Curnier functions
  NEW_TEST(AlartCurnierFunctions_test fc3d_AlartCurnierFunctions_test.c)

  #
  if(WITH_FCLIB)
    NEW_TEST(FCLIB_test1 fc3d_writefclib_local_test.c)
    
    NEW_FC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_6)
    NEW_FC_3D_TEST_COLLECTION(TEST_ADMM_COLLECTION_3)

    # NEW_FC_3D_TEST_HDF5(Capsules-i125-1213.hdf5 SICONOS_FRICTION_3D_NSGS)
    # NEW_FC_3D_TEST_HDF5(Capsules-i125-1213.hdf5 SICONOS_FRICTION_3D_ADMM 1e-10 0
    #   0 0 0
    #   IPARAM SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING)
    
    # NEW_FC_3D_TEST_HDF5(LMGC_100_PR_PerioBox-i00361-60-03000.hdf5 SICONOS_FRICTION_3D_ADMM 1e-08 100000
    #   0 0 0
    #   IPARAM SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING)

    # NEW_FC_3D_TEST_HDF5(LMGC_100_PR_PerioBox-i00361-60-03000.hdf5 SICONOS_FRICTION_3D_ADMM 1e-08 100000
    #   0 0 0
    #   IPARAM SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING
    #   IPARAM SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE  SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE)


    NEW_GFC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_2)
    
    NEW_GFC_3D_TEST_COLLECTION(TEST_WR_COLLECTION_2)
    
    NEW_GFC_3D_TEST_COLLECTION(TEST_NSN_COLLECTION_2)
    
  # ---------------------------------------------------
  # --- Rolling friction contact problem formulation ---
  # ---------------------------------------------------
  
  NEW_RFC_3D_TEST_COLLECTION(TEST_NSGS_COLLECTION_1)
  endif()

  #===========================================
  # 2D Friction Contact tests
  #===========================================
  ## test 2D dense on two differents files

  NEW_FC_2D_TEST_COLLECTION(TEST_FC2D_COLLECTION_1)

  END_TEST()


  BEGIN_TEST(src/GenericMechanical/test)

  NEW_GMP_TEST_COLLECTION(TEST_NSGS_COLLECTION_1)

  END_TEST()
  #BEGIN_TEST(src/GenericMechanical/test)
  #NEW_TEST(GMP_FAILED GenericMechanical_test1.c)
  #END_TEST()


  # === Variationnal inequalities tests ===

  BEGIN_TEST(src/VI/test)

  
  NEW_TEST(VI_test_collection VI_test_collection_1.c)
  NEW_TEST(VI_fc3d_test_collection VI_fc3d_test_collection_1.c)

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

  # === QP tests ===
  BEGIN_TEST(src/QP/test)

  NEW_TEST(ConvexQP_test_collection ConvexQP_test.c)

  NEW_TEST(ConvexQP_FC3D_test_collection ConvexQP_FC3D_test.c)
 
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
