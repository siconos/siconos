include(tools4tests)

if(WITH_TESTING)

  # If WITH_SYSTEM_SUITESPARSE, suite sparse is an imported target
  # that must sometimes be taken into account by tests.
  if(WITH_SYSTEM_SUITESPARSE)
    set(suitesparse SuiteSparse::CXSparse)
  else()
    set(suitesparse)
  endif()
      

  begin_tests(src/tools/test)

  new_test(SOURCES test_op3x3.c)

  new_test(SOURCES test_timers_interf.c)

  new_test(SOURCES test_blas_lapack.c)

  #if(HAS_LAPACK_dgesvd) # Some lapack versions miss dgesvd
  new_test(SOURCES test_pinv.c)# DEPS "externals")
  #endif()
  
  new_test(NAME tools_projection SOURCES test_projection.c)

  new_test(SOURCES NumericsArrays.c)

  #  tests for NumericsMatrix
  new_test(SOURCES NM_test.c DEPS "${suitesparse}")

  #  tests for JordanAlgebra
  NEW_TEST(NAME tools_test_JordanAlgebra SOURCES JordanAlgebra_test.c)

  # MUMPS interface tests
  if(WITH_MUMPS)
    new_test(SOURCES NM_MUMPS_test.c)
  endif()
  
  # Specfic tests for SBM matrices 
  new_test(SOURCES SBM_test.c DEPS "${suitesparse}")
  new_test(SOURCES SBCM_to_SBM.c)

  # Specfic tests for sparse matrices 
  new_test(SOURCES SparseMatrix_test.c DEPS "${suitesparse}")

  if(HAS_ONE_LP_SOLVER)
    new_test(SOURCES vertex_problem.c)
  endif(HAS_ONE_LP_SOLVER)

  # ----------- LCP solvers tests -----------
  # Start tests for LCP dir.
  begin_tests(src/LCP/test)

  # Two kinds of tests :
  # * those with existing source file ('standards') and those where sources
  #  --> use new_test function (see details in cmake/tools4tests.cmake)
  #   new_test(NAME <test_name> SOURCES <source file name>)
  # * those generated from a global driver for a collection of solvers and data : run a test
  #  for each couple (data file, solver name).
  #  Usually, solvers list is defined in test_solvers_collection* files and data files list
  #  in data_collection* files.
  #  Use new_tests_collection function as below.
  
  new_test(NAME lcp_test_DefaultSolverOptions SOURCES LinearComplementarity_DefaultSolverOptions_test.c)

  new_tests_collection(
    DRIVER lcp_test_collection.c.in FORMULATION lcp COLLECTION TEST_LCP_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c)
  new_tests_collection(
    DRIVER lcp_test_collection.c.in FORMULATION lcp COLLECTION TEST_LCP_COLLECTION_2
    EXTRA_SOURCES data_collection_2.c)
  new_tests_collection(
    DRIVER lcp_test_collection.c.in FORMULATION lcp COLLECTION TEST_LCP_COLLECTION_3
    EXTRA_SOURCES data_collection_3.c)
  new_tests_collection(
    DRIVER lcp_test_collection.c.in FORMULATION lcp COLLECTION TEST_LCP_COLLECTION_4
    EXTRA_SOURCES data_collection_4.c)
  new_tests_collection(
    DRIVER lcp_test_collection.c.in FORMULATION lcp COLLECTION TEST_LCP_COLLECTION_5
    EXTRA_SOURCES data_collection_5.c)

  # ----------- Relay solvers tests -----------
  # Start tests for Relay dir.
  begin_tests(src/Relay/test)

  new_tests_collection(
    DRIVER relay_test_collection.c.in FORMULATION relay COLLECTION TEST_RELAY_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c)

  # ----------- MLCP solvers tests -----------
  begin_tests(src/MLCP/test)

  new_tests_collection(
    DRIVER mlcp_test_collection.c.in FORMULATION mlcp COLLECTION TEST_MLCP_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c)

  new_tests_collection(
    DRIVER mlcp_test_collection.c.in FORMULATION mlcp COLLECTION TEST_MLCP_COLLECTION_2
    EXTRA_SOURCES data_collection_2.c)

  new_tests_collection(
    DRIVER mlcp_test_collection.c.in FORMULATION mlcp COLLECTION TEST_MLCP_COLLECTION_3
    EXTRA_SOURCES data_collection_3.c)

  # new_tests_collection(
  #   DRIVER mlcp_test_collection.c.in FORMULATION mlcp COLLECTION TEST_MLCP_COLLECTION_4
  #   EXTRA_SOURCES data_collection_4.c)


  if(HAVE_SYSTIMES_H AND WITH_CXX)
    new_test(NAME MLCPtest SOURCES main_mlcp.cpp)
  endif()
  new_test(SOURCES MixedLinearComplementarity_ReadWrite_test.c)

  # ----------- MCP solvers tests -----------
  begin_tests(src/MCP/test)
  new_test(SOURCES MCP_test.c)
  new_test(SOURCES MCP_test1.c)
  new_test(SOURCES MCP_test2.c)
  new_test(SOURCES MCP_2_test1.c)

  # ----------- NCP solvers tests -----------
  begin_tests(src/NCP/test)
  # Generated tests
  # -- List of data files used in NCP tests --
  #  -- List of solvers for NCPs --
  set(SICONOS_NCP_SOLVERS
    SICONOS_NCP_NEWTON_FB_FBLSA
    SICONOS_NCP_NEWTON_MIN_FBLSA
    SICONOS_NCP_PATHSEARCH
    )
  if(HAVE_PATHFERRIS)
    list(APPEND SICONOS_NCP_SOLVERS "SICONOS_NCP_PATH")
  endif()
  
  if(WITH_UNSTABLE_TEST)
  
    # -- Declare the tests --
    foreach(SOLVER IN LISTS SICONOS_NCP_SOLVERS)
      new_tests_collection(DRIVER NCP_ZI1.c.in FORMULATION NCP COLLECTION ${SOLVER})
      if(WARNINGS_LEVEL GREATER 0)
        new_tests_collection(DRIVER NCP_ZIT1.c.in FORMULATION NCP COLLECTION ${SOLVER} SUFFIX _UNSTABLE)
      endif()
    endforeach()
    
  endif()

  #===========================================
  # 3D Friction Contact tests
  #===========================================

  begin_tests(src/FrictionContact/test DEPS "${suitesparse}")
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSGS_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c test_nsgs_1.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSGS_COLLECTION_2
    EXTRA_SOURCES data_collection_2.c test_nsgs_1.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSGS_COLLECTION_FREEZE_2
    EXTRA_SOURCES data_collection_2.c test_nsgs_freeze_1.c)

  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSGS_COLLECTION_3
    EXTRA_SOURCES data_collection_3.c test_nsgs_3.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSGS_COLLECTION_QUARTIC
    EXTRA_SOURCES rover_collection.c test_nsgs_quartic.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSGS_COLLECTION_5
    EXTRA_SOURCES data_collection_3.c test_nsgs_5.c)
  
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_ADMM_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c test_admm_1.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_ADMM_COLLECTION_2
    EXTRA_SOURCES data_collection_2.c test_admm_2.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSN_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c test_nsn_1.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSN_COLLECTION_2
    EXTRA_SOURCES data_collection_3.c test_nsn_2.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSN_COLLECTION_3
    EXTRA_SOURCES data_collection_3.c test_nsn_3.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_VI_BASED_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c test_vi_based_1.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_DSFP_COLLECTION
    EXTRA_SOURCES data_collection_3.c test_dsfp.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_FP_COLLECTION_1
    EXTRA_SOURCES data_collection_3.c test_fp_1.c)
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_FP_COLLECTION_2
    EXTRA_SOURCES data_collection_3.c test_fp_2.c)

  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_PROX_COLLECTION_1
    EXTRA_SOURCES data_collection_3.c test_prox_1.c)

  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_QUARTIC_COLLECTION_1
    EXTRA_SOURCES data_collection_5.c test_quartic_1.c)
    
  # --- LMGC driver ---
  new_test(SOURCES fc3d_newFromFortranData.c)
  new_test(SOURCES fc3d_LmgcDriver_test1.c)
  new_test(SOURCES fc3d_LmgcDriver_test2.c)
  new_test(SOURCES fc3d_LmgcDriver_test3.c)
  new_test(SOURCES fc3d_LmgcDriver_test4.c)
  new_test(SOURCES fc3d_LmgcDriver_test5.c)

  # ---------------------------------------------------
  # --- Global friction contact problem formulation ---
  # ---------------------------------------------------

  new_tests_collection(
    DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_FIRST_ORDER_COLLECTION_1
    EXTRA_SOURCES data_collection_gfc3d_1.c test_first_order_gfc3d_1.c)
  new_tests_collection(
    DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_WR_COLLECTION_1
    EXTRA_SOURCES data_collection_gfc3d_1.c test_solvers_wr_gfc3d_1.c)
  new_tests_collection(
    DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_NSN_COLLECTION_1
    EXTRA_SOURCES data_collection_gfc3d_1.c test_nsn_gfc3d_1.c)
  new_tests_collection(
    DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_IPM_COLLECTION_1
    EXTRA_SOURCES data_collection_gfc3d_1.c test_ipm_gfc3d_1.c )
  new_tests_collection(
    DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_ADMM_COLLECTION_1
    EXTRA_SOURCES data_collection_gfc3d_1.c test_admm_gfc3d_1.c )

  
  new_tests_collection(
    DRIVER gfc2d_test_collection.c.in FORMULATION gfc2d COLLECTION TEST_FIRST_ORDER_COLLECTION_1
    EXTRA_SOURCES data_collection_gfc2d_1.c test_first_order_gfc2d_1.c )

  # Alart Curnier functions
  new_test(NAME AlartCurnierFunctions_test SOURCES fc3d_AlartCurnierFunctions_test.c)

  # ---------------------------------------------------
  # --- Rolling friction contact problem formulation ---
  # ---------------------------------------------------
  
  new_tests_collection(
    DRIVER rfc3d_test_collection.c.in  FORMULATION rolling_fc3d COLLECTION TEST_FIRST_ORDER_COLLECTION
    EXTRA_SOURCES data_collection_rfc3d.c test_first_order_rfc3d_1.c )
  new_tests_collection(
    DRIVER grfc3d_test_collection.c.in  FORMULATION grfc3d COLLECTION TEST_IPM_COLLECTION_1
    EXTRA_SOURCES data_collection_grfc3d.c test_ipm_grfc3d_1.c )
      
  if(WITH_FCLIB)

    new_test(NAME FCLIB_test1 SOURCES fc3d_writefclib_local_test.c DEPS FCLIB::fclib)

    new_tests_collection(
      DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_NSGS_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_fclib.c test_nsgs_1.c DEPS FCLIB::fclib
      HDF5 ON
      )

    new_tests_collection(
      DRIVER fc_test_collection.c.in FORMULATION fc3d COLLECTION TEST_ADMM_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_fclib.c test_admm_1.c DEPS FCLIB::fclib
      HDF5 ON
      )

    new_tests_collection(
      DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_FIRST_ORDER_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_gfc3d_fclib.c test_first_order_gfc3d_1.c DEPS FCLIB::fclib
      HDF5 ON
      )
    
    new_tests_collection(
      DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_ADMM_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_gfc3d_3.c test_admm_gfc3d_1.c DEPS FCLIB::fclib
      HDF5 ON
      )
    new_tests_collection(
      DRIVER gfc3d_test_collection.c.in FORMULATION gfc3d COLLECTION TEST_ADMM_COLLECTION_FCLIB_PA
      EXTRA_SOURCES data_collection_gfc3d_fclib.c test_admm_gfc3d_1.c DEPS FCLIB::fclib
      HDF5 ON
      )
 
    new_tests_collection(
      DRIVER gfc3d_test_collection.c.in  FORMULATION gfc3d COLLECTION TEST_WR_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_gfc3d_fclib_1.c test_solvers_wr_gfc3d_fclib.c DEPS FCLIB::fclib
      HDF5 ON
      )

    new_tests_collection(
      DRIVER gfc3d_test_collection.c.in  FORMULATION gfc3d COLLECTION TEST_NSN_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_gfc3d_fclib.c test_nsn_gfc3d_1.c DEPS FCLIB::fclib
      HDF5 ON
      )
    new_tests_collection(
      DRIVER gfc3d_test_collection.c.in  FORMULATION gfc3d COLLECTION TEST_IPM_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_gfc3d_fclib.c test_ipm_gfc3d_1.c DEPS FCLIB::fclib
      HDF5 ON
      )
    new_tests_collection(
      DRIVER gfc3d_test_collection.c.in  FORMULATION gfc3d_nonsmooth COLLECTION TEST_IPM_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_gfc3d_fclib.c test_ipm_gfc3d_1.c DEPS FCLIB::fclib
      HDF5 ON
      )

    # ---------------------------------------------------
    # --- Rolling friction contact problem formulation ---
    # ---------------------------------------------------
    
    new_tests_collection(
      DRIVER rfc3d_test_collection.c.in  FORMULATION rolling_fc3d COLLECTION TEST_FIRST_ORDER_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_rfc3d_fclib.c test_first_order_rfc3d_1.c DEPS FCLIB::fclib)

    new_tests_collection(
      DRIVER grfc3d_test_collection.c.in  FORMULATION grfc3d COLLECTION TEST_IPM_COLLECTION_FCLIB
      EXTRA_SOURCES data_collection_grfc3d_fclib.c test_ipm_grfc3d_1.c DEPS FCLIB::fclib
      HDF5 ON
      )
    
  endif()

  #===========================================
  # 2D Friction Contact tests
  #===========================================
  # test 2D dense, all solvers but enum
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc2d COLLECTION TEST_FC2D_COLLECTION_1
    EXTRA_SOURCES data_collection_fc2d_1.c test_fc2d_1.c)

  
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc2d COLLECTION TEST_FC2D_COLLECTION_2
    EXTRA_SOURCES data_collection_fc2d_2.c test_fc2d_2.c)
  # test fc2d dense, for enum solvers. Files from data_collection_fc2d_1.c lead to timeout.
  
  
  new_tests_collection(
    DRIVER fc_test_collection.c.in FORMULATION fc2d COLLECTION TEST_FC2D_COLLECTION_ENUM
    EXTRA_SOURCES data_collection_fc2d_enum.c test_fc2d_enum.c)

  #===========================================
  # Generic mechanical tests
  #===========================================
  begin_tests(src/GenericMechanical/test)

  new_tests_collection(
    DRIVER gmp_test_collection.c.in FORMULATION gmp COLLECTION TEST_NSGS_COLLECTION_1
    EXTRA_SOURCES data_collection_1.c test_solvers_1.c)

  # ----------- Variationnal inequalities solvers tests -----------
  begin_tests(src/VI/test)

  new_test(SOURCES VI_test_collection_1.c)
  new_test(SOURCES VI_fc3d_test_collection_1.c)

  set(SICONOS_VI_SOLVERS
    SICONOS_VI_BOX_QI
    SICONOS_VI_BOX_AVI_LSA
    )
  if(HAVE_PATHFERRIS)
    list(APPEND SICONOS_VI_SOLVERS "SICONOS_VI_BOX_PATH")
  endif()
  
  if(WARNINGS_LEVEL GREATER 0)
    foreach(SOLVER IN LISTS SICONOS_VI_SOLVERS)
      new_tests_collection(DRIVER VI_ZI1.c.in FORMULATION vi COLLECTION ${SOLVER} SUFFIX I1 )
      new_tests_collection(DRIVER VI_ZIT1.c.in FORMULATION vi COLLECTION ${SOLVER} SUFFIX IT1 )
    endforeach()
  endif()

  # ----------- QP solvers tests -----------
  begin_tests(src/QP/test)

  new_test(NAME ConvexQP_test_collection SOURCES ConvexQP_test.c)
  new_test(NAME ConvexQP_FC3D_test_collection SOURCES  ConvexQP_FC3D_test.c)
  
  # ----------- AVI solvers tests -----------
  begin_tests(src/AVI/test)

  if(HAS_ONE_LP_SOLVER)
    new_test(NAME AVI_twisting SOURCES implicit_twisting.c)
  endif()

  # ----------- SOCP solvers tests -----------
  begin_tests(src/SOCP/test)
  new_test(SOURCES soclcp_test1.c)
  new_test(SOURCES soclcp_test2.c)
  new_test(SOURCES soclcp_test3.c)
  # timeout on all machines, see
  # http://cdash-bipop.inrialpes.fr/testSummary.php?project=1&name=SOCLCP_test4&date=2015-09-03
  # Feel free to remove this once it is fixed --xhub
  #new_test(SOURCES soclcp_test4.c)
  #new_test(SOURCES soclcp_test5.c)
  new_test(SOURCES fc3d_to_soclcp.c)

  # ---- Extra conf for ${COMPONENT}-test library ---
  if(TARGET numerics-test)
    if(MSVC)
      # This part should be reviewed by a windows expert ...
      include(WindowsLibrarySetup)
      windows_library_extra_setup("numerics-test" "numerics-test")
    endif()
  endif()

  # For SuiteSparse and SiconosLapack.h 
  target_link_libraries(numerics-test PUBLIC externals)
  target_link_libraries(numerics-test PUBLIC LAPACK::LAPACK)
  #target_include_directories(numerics-test PUBLIC ${CMAKE_SOURCE_DIR}/externals/blas_lapack)

endif()
