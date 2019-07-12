include(tools4tests)

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
  NEW_TEST(MCP_test_2 MCP_test2.c)
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
  target_link_libraries(numerics-test PRIVATE ${COMPONENT})
  target_link_libraries(numerics-test PRIVATE ${${COMPONENT}_LINK_LIBRARIES})
endif()
