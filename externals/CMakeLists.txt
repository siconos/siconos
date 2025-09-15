#================================================================
# cmake utilities to build external component
#================================================================

# Component name (i.e. target name)
message("-- Set up for ${PROJECT_NAME}_externals library ...\n")

include(ComponentSetup)

if(WITH_FORTRAN)
  ### Fortran part ###
  create_siconos_component(externals_fortran)
  
  # - netlib
  target_sources(externals_fortran PRIVATE
    netlib/odepack/opkda1.f
    netlib/odepack/opkda2.f
    netlib/odepack/opkdmain.f
    )

  # add fallow-argument-mismatch' option to compile fortran when compiler version is 10 or greater
  # See https://scivision.co/gfortran-type-mismatch-error/

  set_source_files_properties(netlib/odepack/opkda1.f netlib/odepack/opkdmain.f
    PROPERTIES COMPILE_OPTIONS "-fallow-argument-mismatch")
  # target_include_directories(externals_fortran INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/netlib>) # For SiconosNetlib.h
  target_include_directories(externals_fortran INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tools>) # For SiconosFortran.h
 
  # - hairer
  # https://www.unige.ch/~hairer/software.html
  target_sources(externals_fortran PRIVATE
    hairer/linsp.f
    hairer/decsol.f
    hairer/hem5.f)
  set_source_files_properties(hairer/hem5.f PROPERTIES COMPILE_OPTIONS "-fallow-argument-mismatch")
  # HEM5 calls DGETRF and DGETRS from LAPACK
  target_link_libraries(externals_fortran PRIVATE $<BUILD_INTERFACE:LAPACK::LAPACK>)

  #target_compile_options(externals PRIVATE "-Wno-unused-local-typedef") # Else boost bindings lead to many warnings
  #target_compile_options(externals_fortran PRIVATE "-Wno-unused-variable")
  #target_compile_options(externals_fortran PRIVATE "-Wno-unused-parameter")
  #target_compile_options(externals_fortran PRIVATE "-fallow-argument-mismatch")
  # --- ql0001 ---
  #if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/optim_misc/ql0001/ql0001.f")
  target_sources(externals_fortran PRIVATE optim_misc/ql0001/ql0001.f optim_misc/qnb.f)
  set(HAVE_QL0001 TRUE CACHE INTERNAL "True if ql0001 is available in externals.")


  # Fortran module used to defines iso_c_bindings API to all f77 routines
  target_sources(externals_fortran PRIVATE siconos_fortran2c.f90)
  
  set_property(TARGET externals_fortran PROPERTY LINKER_LANGUAGE Fortran)
  target_compile_options(externals_fortran PRIVATE "-w")

endif()

# Windows stuff --> should probably be reviewed and updated
#include(WindowsExternalsSetup)

# -- create/setup component target --
create_siconos_component(externals)
target_compile_options(externals PRIVATE "-w")
# - lumod part -
target_sources(externals PRIVATE
  SOL/lumod-c/commonlib.c SOL/lumod-c/lumod_dense.c)
# For numerics: 
target_include_directories(externals INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SOL/lumod-c>)

# - sort -
# --- Sort ---
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/sort/sort.h")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/sort/sort_common.h")
    set(HAVE_SORT TRUE CACHE INTERNAL "True if sort is available in externals.")
    target_include_directories(externals INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sort>)
    # Required by NumericsSparseMatrix
  endif()
endif()

# --- MA57  ---
#SET(WITH_MA57 1)
# - Extra setup for the component -
if(WITH_MA57) 
  target_sources(externals PRIVATE
    lbl/src/lbl_lib.c
    lbl/src/lbl_mem.c
    lbl/src/ma27_lib.c
    lbl/src/ma57_lib.c)
  get_sources(metisMa57 DIRS lbl/metis4)
  # metis for ma57
  target_sources(externals PRIVATE ${metisMa57_SRCS})
  # find_package(Metis) Metis 5  is not compatible with MA57. We include metis4
  # if(METIS_FOUND)
  #   target_link_libraries(externals PUBLIC Metis::Metis)
  # endif()
  # --> todo: write Find metis for cmake
  target_include_directories(externals PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lbl/include>)
endif()
 
# - Extra setup for the component -
if(GAMSCAPI_FOUND AND CMAKE_DL_LIBS) # needed by GAMS (?). Here or in numerics ?
  target_link_libraries(externals PRIVATE ${CMAKE_DL_LIBS})
endif()


if(WITH_FORTRAN)
  target_link_libraries(externals PUBLIC externals_fortran)
  target_compile_definitions(externals_fortran INTERFACE HAS_FORTRAN)
endif()

set_target_properties(externals PROPERTIES LINKER_LANGUAGE C)
target_compile_features(externals PUBLIC c_std_11)

# externals lib and tests need blas and lapack
target_link_libraries(externals PUBLIC $<BUILD_INTERFACE:BLAS::BLAS>)
target_include_directories(externals PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack>)

# ---- Installation ----
if(WITH_FORTRAN)
  siconos_component_install_setup(externals_fortran NODOC)
endif()
siconos_component_install_setup(externals NODOC)
install(FILES tools/SiconosFortran.h DESTINATION include/siconos)

# - Optional OCE Renderer setup (install only) -
if(WITH_OpenCASCADE)
  include(occ_renderer)
endif()

if(WITH_SERIALIZATION AND INSTALL_EXTERNAL_HEADERS AND NOT WITH_SYSTEM_BOOST_SERIALIZATION)
  install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/boost/serialization/unordered_collections_load_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/boost/serialization/unordered_collections_save_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/boost/serialization/unordered_map.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_serialization/boost/serialization/unordered_set.hpp
    DESTINATION include/siconos/boost/serialization)
endif()

# --- tests ---
include(externals_tests.cmake)

