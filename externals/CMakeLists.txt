#================================================================
# cmake utilities to build external component
#================================================================
#
# The objective is to call component_setup to create the target <COMPONENT>.
# Before, it's necessary to set:
# 
# - COMPONENT component name
# - <COMPONENT>_DIRS: the list of paths (relative to CMAKE_CURRENT_SOURCE_DIR) that
#   contain source files
# - <COMPONENT>_EXCLUDE_SRCS: the list of files, in <COMPONENT>_DIRS, that must be excluded
#   from build.
# - <COMPONENT>_INTERFACE_INCLUDE_DIRECTORIES: a list of directories
#   to populate the interface of the target for include directories at build time


# Component name (i.e. target name)
set(COMPONENT externals)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# ------ source directories for current component ------
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

# boost extras, when c++ is on only
if(WITH_CXX)
  list(APPEND ${COMPONENT}_DIRS numeric_bindings boost_contribs)
  if(NOT WITH_SYSTEM_BOOST_SERIALIZATION)
    list(APPEND ${COMPONENT}_DIRS boost_serialization)
  endif()
endif()

# -- Source files to be excluded from build --
set(${COMPONENT}_EXCLUDE_SRCS "SOL/lumod-c/sparselib.c" "SOL/lumod-c/lumod_sparse.c")

# --- Sort ---
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/sort/sort.h")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/sort/sort_common.h")
    set(HAVE_SORT TRUE CACHE INTERNAL "True if sort is available in externals.")
  endif()
endif()

# --- ql0001 ---
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/optim_misc/ql0001/ql0001.f")
  set(HAVE_QL0001 TRUE CACHE INTERNAL "True if ql0001 is available in externals.")
endif()

# -- Documentation --
# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
# No docs for externals.
set(${COMPONENT}_EXCLUDE_DOXY ${${COMPONENT}_DIRS})


# TEST TEST TEST
# ------ include interface ------
# What is needed at build time
# by other targets (numerics etc) to compile with externals.
# It means that a call to
#  target_link_libraries(truc PRIVATE externals)
# will imply -I<dirs> with dirs listed in
# ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES.
set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES
  netlib
   netlib/odepack  # For odepack.h in tests and kernel
   tools           # For SiconosCompat.h in numerics 
   blas_lapack     # For SiconosLapack.h in numerics
   SOL/lumod-c # For lumod_dense.h in numerics
   hairer # For kernel
   )

if(WITH_CXX)
  list(APPEND ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES
    numeric_bindings boost_contribs) # SiconosAlgebra needs this
endif()

if(HAVE_SORT)
  list(APPEND ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES sort) # NumericsSparseMatrix needs this
endif()

if(NOT WITH_SYSTEM_SUITESPARSE)
  # else suite sparse is an external dependency of numerics and handled in numerics
  # CMakeLists.txt.
  list(APPEND ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES SuiteSparse/CXSparse)
endif()
# TEST TEST TEST

# ---- Final setup for the library ----

# Windows stuff --> should probably be reviewed and updated
include(WindowsExternalsSetup)

# -- create/setup component target --
include(ComponentSetup)
create_siconos_component(${COMPONENT})

# --- python bindings ---
if(WITH_${COMPONENT}_PYTHON_WRAPPER)
  add_subdirectory(swig)
endif()
  
if(WITH_PYTHON_WRAPPER)
  add_dependencies(${COMPONENT} ${COMPONENT}_docstrings)
endif()  

# - Extra setup for the component -
if(GAMSCAPI_FOUND AND CMAKE_DL_LIBS) # needed by GAMS (?). Here or in numerics ?
  target_link_libraries(externals PRIVATE ${CMAKE_DL_LIBS})
endif()

# If dtrtrs is not defined in blas/lapack (which is the case in atlas) add it to the sources
if(NOT HAS_LAPACK_dtrtrs)
  target_sources(externals PRIVATE blas_lapack/clapack_dtrtrs.c)
endif()

if(BUILD_AS_CPP) 
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX
    COMPILE_FLAGS -fpermissive)
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE CXX)
  # We require C++ 11. PUBLIC to propagate to numerics.
  target_compile_features(externals PUBLIC cxx_std_11)
else()
  set_target_properties(${COMPONENT} PROPERTIES LINKER_LANGUAGE C)
  target_compile_features(externals PUBLIC c_std_11)
endif()

# externals lib needs blas while tests need Blas AND lapack.
# Link only at build time.
target_link_libraries(externals PUBLIC $<BUILD_INTERFACE:BLAS::BLAS>)
target_include_directories(externals PRIVATE blas_lapack)

# We don't want warnings from externals libraries.
target_compile_options(externals PRIVATE "-w")
target_compile_options(externals PRIVATE "-Wno-unused-variable")

# - Extras, optional -
# --- SuiteSparse ---
# It is allowed to switch between system or local suitesparse,
# depending on the chosen user option.
# 'local' : use WITH_SYSTEM_SUITESPARSE=OFF, 'system' WITH_SYSTEM_SUITESPARSE=ON (default).
# For 'local' suite sparse, some extra setup is required in externals,
# while system suitesparse is searched for and configured
# during numerics component setup.
if(NOT WITH_SYSTEM_SUITESPARSE)
  include(suite_sparse_setup)
endif()


# ---- Installation ----
# Call siconos_component_install_setup(<COMPONENT>)
# to prepare installation of the current target.
# Before, it's necessary to set:
# 
# - <COMPONENT>_INSTALL_INTERFACE_INCLUDE_DIRECTORIES with all directories
#    that contain headers files that must be installed.
# 
# common install
# No headers installed for externals.
set(${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES)
# set(${COMPONENT}_HDRS_EXCLUDE tools/SiconosCompat.h)
siconos_component_install_setup(${COMPONENT})

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


# # List of directories of headers not to be installed
# set(${COMPONENT}_HDRS_EXCLUDE_DIR PATH_SDK/include SOL/lumod-c)
# if(NOT INSTALL_EXTERNAL_HEADERS)
#   list(APPEND ${COMPONENT}_HDRS_EXCLUDE_DIR
#     blas_lapack
#     hairer
#     netlib/dftemplates
#     netlib/odepack
#     optim_misc
#     optim_misc/ql0001
#     sort
#     boost_contribs
#     tools)
# endif()

# # List of specific headers not to be installed
# set(${COMPONENT}_HDRS_EXCLUDE tools/SiconosCompat.h)

