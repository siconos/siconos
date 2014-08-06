#
# Here comes Numerics
#
PROJECT(SiconosNumerics C CXX Fortran)

include(FortranCInterface)
# Set module files directory (i.e. where .mod will be created)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/Modules)
# Define a Fortran interface file (FCMangle.hpp)
#FortranCInterface_VERIFY(QUIET)
if(NOT FortranCInterface_VERIFIED_C)
  configure_file(../DefaultMangle.h.in ${CMAKE_BINARY_DIR}/FCMangle.h)
else()
  FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/FCMangle.h
    MACRO_NAMESPACE "myF2C" 
    SYMBOL_NAMESPACE "myF2C"
    )
endif()
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
  hairer
  netlib/odepack
  netlib/dftemplates
  utils
  utils/BlasLapackInterfaces
  NonSmoothTools
  NonSmoothTools/test-utils
  NonSmoothSolvers
  NonSmoothSolvers/AVI
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
  NonSmoothSolvers/VI
)
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
INSTALL(FILES ${CMAKE_BINARY_DIR}/FCMangle.h DESTINATION include/Siconos/Numerics)

#
# compilation setup 
#

compile_with(BLAS REQUIRED)
compile_with(LAPACK REQUIRED)

include_directories(utils/test)

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
  #  APPEND_C_FLAGS("-DWITH_FCLIB")	
ENDIF()

# set extra CFLAGS

IF(C_HAVE_C99 AND NOT MSVC)
  APPEND_C_FLAGS("-std=c99")
ELSEIF(C_HAVE_C99 AND NOT MSVC)
  IF(C_HAVE_XC99 AND NOT MSVC)
    APPEND_C_FLAGS("-xc99")
  ENDIF(C_HAVE_XC99 AND NOT MSVC)
ENDIF(C_HAVE_C99 AND NOT MSVC)

IF(C_HAVE_WALL)
  APPEND_C_FLAGS("-Wall")
ENDIF(C_HAVE_WALL)

IF(C_HAVE_WEXTRA)
 APPEND_C_FLAGS("-Wextra -Wno-unused-parameter")
ENDIF(C_HAVE_WEXTRA)

IF(C_HAVE_IMPL)
  APPEND_C_FLAGS("-Werror=implicit-function-declaration")
ENDIF(C_HAVE_IMPL)

IF(C_HAVE_UNREACH AND NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  APPEND_C_FLAGS("-Werror=unreachable-code")
ENDIF(C_HAVE_UNREACH AND NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

IF(C_HAVE_CONV)
   APPEND_C_FLAGS("-Werror=conversion -Wno-sign-conversion")
   APPEND_C_FLAGS("-Wno-error=sign-conversion -Wno-error=shorten-64-to-32")
ENDIF(C_HAVE_CONV)

# too many errors right now ...
#IF(C_HAVE_MISS)
#  APPEND_C_FLAGS("-Wmissing-prototypes")
#ENDIF(C_HAVE_MISS)

IF(CMAKE_CXX_COMPILER)

 IF(NOT CMAKE_COMPILER_IS_GNUCXX)
  IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
   # Disable warnings with intel compiler due (mainly) to visitors visit function overloading
   IF(CXX_HAVE_DIAG_DISABLE_654)
    APPEND_CXX_FLAGS("-diag-disable 654")
   ENDIF(CXX_HAVE_DIAG_DISABLE_654)
   # Error on intel compiler, see: http://software.intel.com/en-us/forums/showthread.php?t=65041
   # This issue have been solved with ICC >= 12.1
   if(NOT ICCOK)
    if(CXX_HAVE_D__ALIGNED__IGNORED)
     APPEND_CXX_FLAGS("-D__aligned__=ignored")
    endif(CXX_HAVE_D__ALIGNED__IGNORED)
   endif(NOT ICCOK)
  endif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
 endif(NOT CMAKE_COMPILER_IS_GNUCXX)
 # way too verbose with MSVC
 IF(CXX_HAVE_WALL AND NOT MSVC)
  APPEND_CXX_FLAGS("-Wall")
 ENDIF(CXX_HAVE_WALL AND NOT MSVC)
 IF(CXX_HAVE_WEXTRA)
  APPEND_CXX_FLAGS("-Wextra -Wno-unused-parameter")
 ENDIF(CXX_HAVE_WEXTRA)
 IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   # stupid warning
  APPEND_CXX_FLAGS("-Wno-string-plus-int")
 ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

 IF(CXX_HAVE_WOVERLOADED_VIRTUAL)
   APPEND_CXX_FLAGS("-Woverloaded-virtual")
 ENDIF(CXX_HAVE_WOVERLOADED_VIRTUAL)

 IF(CXX_HAVE_IMPL AND NOT MSVC)
  APPEND_CXX_FLAGS("-Werror=implicit-function-declaration")
 ENDIF(CXX_HAVE_IMPL AND NOT MSVC)

 IF(CXX_HAVE_UNREACH AND NOT MSVC  AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    APPEND_CXX_FLAGS("-Werror=unreachable-code")
 ENDIF(CXX_HAVE_UNREACH AND NOT MSVC  AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

 IF(CXX_HAVE_CONV AND NOT MSVC)
   APPEND_CXX_FLAGS("-Werror=conversion -Wno-sign-conversion")
   APPEND_CXX_FLAGS("-Wno-error=sign-conversion -Wno-error=shorten-64-to-32")
 ENDIF(CXX_HAVE_CONV AND NOT MSVC)

 IF(CXX_HAVE_MISS AND NOT MSVC)
    APPEND_CXX_FLAGS("-Wmissing-declarations")
 ENDIF(CXX_HAVE_MISS AND NOT MSVC)

ENDIF(CMAKE_CXX_COMPILER)


#
# Fedora13 https://fedoraproject.org/wiki/UnderstandingDSOLinkChange
#
IF(UNIX)
 IF(C_HAVE_LINKER_M)
  FIND_LIBRARY(M_FOUND "m")
  SET(M_LIBRARIES ${M_FOUND})
  COMPILE_WITH(M STANDARD)
 ENDIF(C_HAVE_LINKER_M)  
ENDIF(UNIX)
APPEND_CXX_FLAGS("-D_NUMERICS_INTERNAL_CXX_")

IF(LLVM_ANALYSE)
  SET(CMAKE_C_OUTPUT_EXTENSION ".bc")
  SET(CMAKE_C_FLAGS "-emit-llvm")
ENDIF()

# CXX may be needed with MlcpSimplex
SET(SiconosNumerics_LINKER_LANGUAGE C)
IF(MLCPSIMPLEX_FOUND) 
  SET(SiconosNumerics_LINKER_LANGUAGE CXX)
ENDIF(MLCPSIMPLEX_FOUND)

# for dlopen or equivalent (used for plugins)
SET(SiconosNumerics_LINK_LIBRARIES ${SiconosNumerics_LINK_LIBRARIES}
 ${CMAKE_DL_LIBS})

IF(MSVC AND "${CMAKE_Fortran_COMPILER}" MATCHES "gfortran")
  SET(GCC_A "")
  EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} -print-libgcc-file-name
  OUTPUT_VARIABLE GCC_A)
  STRING(REGEX REPLACE "\n.*$" "" GCC_A ${GCC_A})
  GET_FILENAME_COMPONENT(GFORTRAN_DLL_A_DIR ${GCC_A} PATH)
  MESSAGE(STATUS "libgfortran.dll.a dir :: ${GFORTRAN_DLL_A_DIR}")
  STRING(REGEX REPLACE "\n.*$" "" GFORTRAN_DLL_A_DIR ${GFORTRAN_DLL_A_DIR})
  find_library(LIBGFORTRAN_FOUND gfortran PATHS ${GFORTRAN_DLL_A_DIR} NO_DEFAULT_PATH)
  SET(SiconosNumerics_LINK_LIBRARIES ${SiconosNumerics_LINK_LIBRARIES} ${LIBGFORTRAN_FOUND})
  GET_FILENAME_COMPONENT(GFORTRAN_DIR ${CMAKE_Fortran_COMPILER} PATH)
  SET(CMAKE_AR "${GFORTRAN_DIR}/ar")
  APPEND_Fortran_FLAGS("-static -static-libgcc -static-libgfortran") # XXX No test :( -- xhub
ENDIF()

# for LAPACKE, see http://icl.cs.utk.edu/lapack-for-windows/lapack/
IF(MSVC)
  APPEND_C_FLAGS("/DHAVE_LAPACK_CONFIG_H")
  APPEND_C_FLAGS("/DLAPACK_COMPLEX_STRUCTURE")
  APPEND_C_FLAGS("/DADD_")
  APPEND_CXX_FLAGS("/DHAVE_LAPACK_CONFIG_H")
  APPEND_CXX_FLAGS("/DLAPACK_COMPLEX_STRUCTURE")
  APPEND_CXX_FLAGS("/DADD_")
ENDIF(MSVC)

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  STRING(REGEX REPLACE "\\\\" "/" ENV_PATH "$ENV{PATH}")
  STRING(REGEX REPLACE "\;" "\\\;" ENV_PATH "${ENV_PATH}")
  STRING(REGEX REPLACE "\"" "" ENV_PATH "${ENV_PATH}")
ENDIF()


LIBRARY_PROJECT_SETUP()


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
  NEW_TEST(dr1_radau5 dr1_radau5.f)
  NEW_TEST(dr2_radau5 dr2_radau5.f)
  NEW_TEST(dr_radau dr_radau.f)
  NEW_TEST(dr_radaup dr_radaup.f)
  NEW_TEST(dr_rodas dr_rodas.f)
  NEW_TEST(dr_seulex dr_seulex.f)
  END_TEST()
  
  BEGIN_TEST(utils/test)
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
  NEW_TEST(SBMTest4 SBM_test4.c)
  NEW_TEST(SBMTest5 SBM_test5.c)
  NEW_TEST(SparseMatrix0 SparseMatrix_test0.c)
 END_TEST()

 BEGIN_TEST2(NonSmoothSolvers/LCP/test)

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
   "ENUM;LEMKE;CPG;PGS;RPGS;LATIN;LATIN_W;QP;NSQP;AVI_CAOFERRIS;NEWTONMIN;NEWTONFB;NEWTONMINFB;BARD;MURTY;PIVOT")
  IF(PATHFERRIS_FOUND)
   LIST(APPEND SICONOS_LCP_SOLVERS "PATH")
  ENDIF()
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
  SET_LCP_TEST_AS_FAILED("CPS_3" "LEMKE;AVI_CAOFERRIS;PIVOT")

  # PSD matrices and those algo seems not to be a good mix
  SET_LCP_TEST_AS_FAILED("CPS_2;CPS_3" "NSQP;RPGS")

  # lcp_mmc is of size 26, way too much for enum
  RM_TEST2(SICONOS_LCP_ENUM "lcp_mmc.dat")
  # this LCP was put here to show that enum does not work on every LCP, likely
  # due to numerical problems, but works on some system ...
  RM_TEST2(SICONOS_LCP_ENUM "lcp_enum_fails.dat")

  # special tests
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol.dat)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_Pang_isolated_sol_perturbed.dat)
  SET(test-LCP_ENUM-lcp_Pang_isolated_sol_perturbed_PROPERTIES WILL_FAIL TRUE)
  NEW_LCP_TEST(SICONOS_LCP_ENUM lcp_inf_sol_perturbed.dat)

  NEW_TEST(LCP_DefaultSolverOptionstest LinearComplementarity_DefaultSolverOptions_test.c)

  END_TEST(NonSmoothSolvers/LCP/test)





  BEGIN_TEST(NonSmoothSolvers/LinearSystem/test)

  NEW_LS_TEST(SICONOS_LS_0 ls_trivial.dat)
  SET(test-LS_0-ls_inf_sol_perturbed_PROPERTIES WILL_FAIL TRUE)
  NEW_LS_TEST(SICONOS_LS_0 ls_inf_sol_perturbed.dat)
  END_TEST()

  BEGIN_TEST2(NonSmoothSolvers/Relay/test)

  SET(DATA_SET "relay1.dat;relay_2x2.dat;relay_4x4.dat;relay_simple2.dat;step_1x1.dat;step_2x2.dat;step_4x4.dat")
  SET(SICONOS_RELAY_SOLVERS "ENUM;LEMKE;PGS;AVI_CAOFERRIS")
  IF(PATHFERRIS_FOUND)
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

  IF(PATHFERRIS_FOUND)
    NEW_TEST(Relaytest3 relay_test3.c)
  ENDIF(PATHFERRIS_FOUND)

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


  BEGIN_TEST(NonSmoothSolvers/MLCP/test)
  IF(HAVE_SYSTIMES_H)
    NEW_TEST(MLCPtest main_mlcp.cpp)
  ENDIF(HAVE_SYSTIMES_H)
  NEW_TEST(ReadWrite_MLCPtest MixedLinearComplementarity_ReadWrite_test.c)
  END_TEST()

  BEGIN_TEST(NonSmoothSolvers/MCP/test)
  NEW_TEST(MCPtest MCP_test.c)
  NEW_TEST(MCPtest1 MCP_test1.c)
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
  
  IF(PATHFERRIS_FOUND)
    NEW_TEST(FC3Dtest6 FrictionContact3D_test6.c)
    SET(FC3Dtest6_PROPERTIES WILL_FAIL TRUE)
  ENDIF(PATHFERRIS_FOUND)
  
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
  NEW_TEST(FC3Dtest22 FrictionContact3D_test22.c)
  NEW_TEST(FC3Dtest23 FrictionContact3D_test23.c)
  
  NEW_TEST(FC3Dtest30 FrictionContact3D_test30.c)
  NEW_TEST(FC3Dtest31 FrictionContact3D_test31.c)
  NEW_TEST(FC3Dtest32 FrictionContact3D_test32.c)
  
  NEW_TEST(FC3Dtest33 FrictionContact3D_test33.c) # DSFP converges with specific rho

  NEW_TEST(FC3Dtest3400 FrictionContact3D_test3400.c) # EG 
  NEW_TEST(FC3Dtest35 FrictionContact3D_test35.c)     # EG

  NEW_TEST(FC3Dtest36 FrictionContact3D_test36.c) #TFP with NSGS and projection on cylinder

  #NEW_TEST(FC3Dtest37 FrictionContact3D_test37.c) #TFP with ProjectedGradientOnCylinder is not working ...

  #NEW_TEST(FC3Dtest38 FrictionContact3D_test38.c) # HP is not converging

  NEW_TEST(FC3Dtest40 FrictionContact3D_test40.c) # VI_EG 
  NEW_TEST(FC3Dtest41 FrictionContact3D_test41.c) # VI_FPP
  NEW_TEST(FC3Dtest42 FrictionContact3D_test42.c) # VI_FPP
  NEW_TEST(FC3Dtest43 FrictionContact3D_test43.c) # DSFP converges with specific rho
  NEW_TEST(FC3Dtest44 FrictionContact3D_test44.c) # VI_EG
  NEW_TEST(FC3Dtest45 FrictionContact3D_test45.c) # VI_EG
  
  NEW_TEST(FC3Dtest46 FrictionContact3D_test46.c) # FPP
  NEW_TEST(FC3Dtest47 FrictionContact3D_test47.c) # EG

  SET(FC3Dtest50_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(FC3Dtest50 FrictionContact3D_test50.c)
  

  IF(WITH_MUMPS) # uses sparse Alart & Curnier
    NEW_TEST(FC3Dtest60 FrictionContact3D_test60.c)
  ENDIF(WITH_MUMPS)

  NEW_TEST(FC3Dtest70 FrictionContact3D_test70.c) # DSFP
  NEW_TEST(FC3Dtest71 FrictionContact3D_test71.c) # DSFP
  NEW_TEST(FC3Dtest72 FrictionContact3D_test72.c)	
  NEW_TEST(FC3Dtest73 FrictionContact3D_test73.c)
  NEW_TEST(FC3Dtest74 FrictionContact3D_test74.c)

  IF(WITH_MUMPS) # uses sparse Alart & Curnier
    NEW_TEST(FC3Dtest80 FrictionContact3D_test80.c) # Proximal
    NEW_TEST(FC3Dtest81 FrictionContact3D_test81.c) # Proximal
    NEW_TEST(FC3Dtest83 FrictionContact3D_test83.c) # Proximal
  ENDIF(WITH_MUMPS)


 # NEW_TEST(FC3Dtest82 FrictionContact3D_test82.c) # Proximal

 # NEW_TEST(FC3Dtest90 FrictionContact3D_test90.c) # TFP



  NEW_TEST(FC3DNewFromFortranData FrictionContact3D_newFromFortranData.c)
  NEW_TEST(FC3DLmgcDriver1 FrictionContact3D_LmgcDriver_test1.c)
  NEW_TEST(FC3DLmgcDriver2 FrictionContact3D_LmgcDriver_test2.c)
  NEW_TEST(FC3DLmgcDriver3 FrictionContact3D_LmgcDriver_test3.c)
IF(MUMPS_FOUND)
  NEW_TEST(FC3DLmgcDriver4 FrictionContact3D_LmgcDriver_test4.c)
ENDIF(MUMPS_FOUND)
  NEW_TEST(FC3DLmgcDriver5 FrictionContact3D_LmgcDriver_test5.c)

  NEW_TEST(FC3Dtest100 FrictionContact3D_test100.c)
  NEW_TEST(FC3Dtest110 FrictionContact3D_test110.c)

  SET(FC3Dtest120_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(FC3Dtest120 FrictionContact3D_test120.c) # TFP
  NEW_TEST(FC3Dtest121 FrictionContact3D_test121.c) # DSFP


  NEW_TEST(FC3Dtest130 FrictionContact3D_test130.c)

  NEW_TEST(FC3Dtest200 FrictionContact3D_test200.c)
  NEW_TEST(FC3Dtest210 FrictionContact3D_test210.c)

  # Global Alart Curnier + Rover ok.
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Example1_Fc3D.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover1039.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover1040.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover1041.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover3865.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover4144.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover4396.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover4493.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover4516.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover4609.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover4613.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover4622.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover9770.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover11035.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC Rover11211.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_QUARTIC FrictionContact3D_1c.dat)
  NEW_FC_TEST(SICONOS_FRICTION_3D_QUARTIC FrictionContact3D_RR_1c.dat)

  # Newton Euler 3D Spheres : some failures with global Alart Curnier
  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC NESpheres_10_1.dat)

  NEW_FC_TEST(SICONOS_FRICTION_3D_LOCALAC NESpheres_30_1.dat)

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


  NEW_TEST(GFC3D_test1 GlobalFrictionContact3D_test1.c)
  NEW_TEST(GFC3D_test2 GlobalFrictionContact3D_test2.c)
  NEW_TEST(GFC3D_test3 GlobalFrictionContact3D_test3.c)
  NEW_TEST(GFC3D_test4 GlobalFrictionContact3D_test4.c)
  NEW_TEST(GFC3D_test5 GlobalFrictionContact3D_test5.c)


  SET(GFC3D_test6_PROPERTIES WILL_FAIL TRUE)
  NEW_TEST(GFC3D_test6 GlobalFrictionContact3D_test6.c)
  
  NEW_TEST(GFC3D_test7 GlobalFrictionContact3D_test7.c)
  
  NEW_TEST(GFC3D_test10 GlobalFrictionContact3D_test10.c)
  NEW_TEST(GFC3D_test11 GlobalFrictionContact3D_test11.c)
  NEW_TEST(GFC3D_test12 GlobalFrictionContact3D_test12.c)
  NEW_TEST(GFC3D_test13 GlobalFrictionContact3D_test13.c)
  NEW_TEST(GFC3D_test14 GlobalFrictionContact3D_test14.c)
  NEW_TEST(GFC3D_test14bis GlobalFrictionContact3D_test14bis.c)
  NEW_TEST(GFC3D_test15 GlobalFrictionContact3D_test15.c)
  NEW_TEST(GFC3D_test16 GlobalFrictionContact3D_test16.c)

  
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
  
  IF(HAS_LAPACK_DGESVD)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP0.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP1.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP2.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP3.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP4.dat 0 0 0 0 0 1)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP5.dat 0 0 0 0 0 1)
#    SET(test-GMP-REDUCED1_3D_QUARTIC-GMP6_PROPERTIES WILL_FAIL TRUE)
    NEW_GMP_TEST(SICONOS_FRICTION_3D_QUARTIC GMP6.dat 0 0 0 0 0 1)
  endif()
  
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
  BEGIN_TEST(NonSmoothSolvers/VI/test)
  NEW_TEST(VI_test0 VI_test.c)
  NEW_TEST(VI_test1 VI_test1.c)
  NEW_TEST(VI_testFC3D1 VI_testFC3D1.c)
  NEW_TEST(VI_testFC3D2 VI_testFC3D2.c)
  NEW_TEST(VI_testFC3D3 VI_testFC3D3.c)
  
  END_TEST()
ENDIF(WITH_TESTING)

IF(MSVC)
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX)
ENDIF(MSVC)

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)

IF(PRINT_VAR)
  MESSAGE(STATUS)
  PRINT_VAR(${PRINT_VAR})
  MESSAGE(STATUS)
ENDIF(PRINT_VAR)

# print all features
FEATURE_SUMMARY(WHAT ALL)

WRITE_NOTES()
