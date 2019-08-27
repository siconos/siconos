#================================================================
# cmake utilities to build external component
#================================================================

set(COMPONENT externals)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# -- source directories for current component --

# What is needed by component to compile ?
# List here all directories that contain sources files
# for current component.
set(${COMPONENT}_DIRS
  hairer
  netlib/dftemplates
  netlib/odepack
  optim_misc
  optim_misc/ql0001
  SOL/lumod-c
  sort
  tools
)

if (NOT USE_SYSTEM_SUITESPARSE)
  LIST(APPEND ${COMPONENT}_DIRS SuiteSparse/CXSparse)
endif()

if(WITH_CXX)
  list(APPEND ${COMPONENT}_DIRS numeric_bindings boost_contribs)
  if(NOT WITH_SYSTEM_BOOST_SERIALIZATION)
    list(APPEND ${COMPONENT}_DIRS boost_serialization)
  endif()
endif()

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS)

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR PATH_SDK/include SOL/lumod-c)
if(NOT INSTALL_EXTERNAL_HEADERS)
  list(APPEND ${COMPONENT}_HDRS_EXCLUDE_DIR
    # blas_lapack
    hairer
    netlib/dftemplates
    netlib/odepack
    optim_misc
    optim_misc/ql0001
    sort
    boost_contribs
    tools)
endif()

# List of specific headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE tools/SiconosCompat.h)

if (USE_SYSTEM_SUITESPARSE)
  list(APPEND ${COMPONENT}_HDRS_EXCLUDE_DIR
    SuiteSparse/CXSparse)
endif()

# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
# No docs for externals.
set(${COMPONENT}_EXCLUDE_DOXY
  # .
  ${${COMPONENT}_DIRS}
  )

SET(${COMPONENT}_EXCLUDE_SRCS "SOL/lumod-c/sparselib.c" "SOL/lumod-c/lumod_sparse.c")

set_source_files_properties(SOL/lumod-c/lumod_dense.c PROPERTIES COMPILE_FLAGS -Wno-unused-variable)
set_source_files_properties(SOL/lumod-c/commonlib.c PROPERTIES COMPILE_FLAGS -Wno-unused-variable)
if (NOT USE_SYSTEM_SUITESPARSE)
  set_source_files_properties(SuiteSparse/CSparse/csparse.c PROPERTIES COMPILE_FLAGS -Wno-unused)
endif()

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE C)
# but CXX may be needed with MlcpSimplex
if(MLCPSIMPLEX_FOUND) 
  set(${COMPONENT}_LINKER_LANGUAGE CXX)
endif()
list(APPEND ${COMPONENT}_LINK_LIBRARIES ${CMAKE_DL_LIBS}) # needed by GAMS
list(APPEND ${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})

# Windows stuff
include(WindowsExternalsSetup)

if(BUILD_AS_CPP)
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX
    COMPILE_FLAGS -fpermissive)
  set(${COMPONENT}_LINKER_LANGUAGE CXX)
endif()

check_c_compiler_flag("-Wno-error=float-conversion" W_no_error_float_conversion_flag)
check_c_compiler_flag("-Wno-error=conversion" W_no_error_conversion_flag)
check_c_compiler_flag("-Wno-conversion" W_no_conversion_flag)
check_c_compiler_flag("-fpermissive" fpermissive_flag)


if(NOT USE_SYSTEM_SUITESPARSE)
  file(GLOB_RECURSE CXSPARSE_FILES ${CMAKE_CURRENT_SOURCE_DIR}  cxsparse_*.c)
  if(BUILD_AS_CPP)
    set(_flags)
    if(W_no_error_float_conversion_flag)
      set(_flags "-Wno-error=float-conversion")
    endif()
    if(W_no_error_conversion_flag)
      set(_flags "${_flags} -Wno-error=conversion")
    endif()
    if(W_no_conversion_flag)
      set(_flags "${_flags} -Wno-conversion")
    endif()
    if(fpermissive_flag)
      set(_flags "${_flags} -fpermissive")
    endif()
    if(_flags)
      set_source_files_properties(${CXSPARSE_FILES} PROPERTIES COMPILE_FLAGS ${_flags})
    endif()
  else()
    set(_flags)
    if(W_no_error_conversion_flag)
      set(_flags "-Wno-error=conversion")
    endif()
    if(W_no_error_float_conversion_flag)
      set(_flags "${_flags} -Wno-error=float-conversion")
    endif()
    if(W_no_conversion_flag)
      set(_flags "${_flags} -Wno-conversion")
    endif()
    if(_flags)
      set_source_files_properties(${CXSPARSE_FILES} PROPERTIES COMPILE_FLAGS ${_flags})
    endif()
  endif()
endif()

include(LibraryProjectSetup)
library_project_setup()

# If dtrtrs is not defined in blas/lapack (which is the case in atlas) add it to the sources
if(NOT HAS_LAPACK_dtrtrs)
  target_sources(externals PRIVATE blas_lapack/clapack_dtrtrs.c)
endif()

# externals needs blas
# externals tests need blas AND lapack
# To be tested later : link externals with blas only, private mode
# link tests with lapack.
target_link_libraries(externals PUBLIC LAPACK::LAPACK) # Later : turn this to private (to be tested ...).
target_include_directories(externals PRIVATE blas_lapack)


# --- specific install ----
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FCMangle.h DESTINATION include/${PROJECT_NAME})

# - Optional OCE Renderer setup (install only) -
if(WITH_OCE)
  include(oce_renderer)
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
include(${COMPONENT}_tests.cmake)


# We don't want warnings from externals libraries.
target_compile_options(externals PRIVATE "-w")
