# ===================================================
# Set main parameters for the Siconos cmake project
#
# ===================================================
include(SiconosTools)

# =========== Windows stuff ... ===========
include(WindowsSiconosSetup)

# --------- CMake project internal variables ---------
# Siconos current version
include(SiconosVersion)

# File used to print tests setup messages.
set(TESTS_LOGFILE ${CMAKE_BINARY_DIR}/tests.log)

# -- Set include directories that are required by ALL components
# and only those!
# Other includes must be specified to individual targets only.
# Current binary dir, for generated headers. Only at build time!
include_directories($<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)

set(tests_timeout 120 CACHE INTERNAL "Limit time for tests (in seconds)")

if(WITH_GIT) # User defined option, default = off
  # Check if git is available
  # and get last commit id (long and short).
  # Saved in SOURCE_ABBREV_GIT_SHA1 and SOURCE_GIT_SHA1
  # These vars are useful for tests logs and 'write_notes' macro.
  find_package(Git)
  if(GIT_FOUND)
    set(CTEST_GIT_COMMAND "${GIT_EXECUTABLE}" )     
    execute_process(COMMAND 
      ${GIT_EXECUTABLE} log -n 1 --pretty=format:%h 
      OUTPUT_VARIABLE SOURCE_ABBREV_GIT_SHA1
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    
    execute_process(COMMAND 
      ${GIT_EXECUTABLE} log -n 1 --pretty=format:%H
      OUTPUT_VARIABLE SOURCE_GIT_SHA1
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endif()

# Save date/time into BUILD_TIMESTAMP var
string(TIMESTAMP BUILD_TIMESTAMP)

# ---- PYTHON SETUP ----
# Python interpreter is always required.
# We force Python3!
# In addition, when WITH_PYTHON_WRAPPER is ON,
# we need python libraries and numpy.
if(${CMAKE_VERSION} VERSION_LESS "3.14")
  # Our FindPython3 is just a copy of the one distributed
  # with cmake = 3.14
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/extras)
endif()

if(WITH_PYTHON_WRAPPER)
  find_package(Python3 COMPONENTS Development Interpreter NumPy REQUIRED)
else()
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
endif()

# For backward compat ...
set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})

get_filename_component(PYTHON_EXE_NAME ${PYTHON_EXECUTABLE} NAME)
if(WITH_PYTHON_WRAPPER OR WITH_DOCUMENTATION)
  include(FindPythonModule)
  # --- xml schema. Used in tests. ---
  if(WITH_XML)
    set(SICONOS_XML_SCHEMA "${CMAKE_SOURCE_DIR}/kernel/swig/SiconosModelSchema-V3.7.xsd")
    if(NOT NO_RUNTIME_BUILD_DEP)
      find_python_module(lxml REQUIRED)
    endif()
  endif()
endif()

# --- End of python conf ---

# Choice of CSparse/CXSparse integer size
# Note FP :  this should be in externals isn't it?
IF(NOT DEFINED SICONOS_INT64)
  IF(NOT SIZE_OF_CSI)
    INCLUDE(CheckTypeSize)
    CHECK_TYPE_SIZE("size_t" SIZE_OF_CSI)
    IF(NOT SIZE_OF_CSI)
      message(FATAL_ERROR
        "Could not get size of size_t, please specify SICONOS_INT64.")
    ENDIF(NOT SIZE_OF_CSI)
  ENDIF(NOT SIZE_OF_CSI)

  IF ("${SIZE_OF_CSI}" EQUAL 8)
    SET(SICONOS_INT64 TRUE)
  ELSE ("${SIZE_OF_CSI}" EQUAL 8)
    SET(SICONOS_INT64 FALSE)
  ENDIF ("${SIZE_OF_CSI}" EQUAL 8)
ENDIF()

# =========== install setup ===========

# Set directory used to save cmake config files
# required to use Siconos (e.g. to call find_package(siconos) )
set(ConfigPackageLocation lib/cmake/siconos-${SICONOS_VERSION})

# Provides install directory variables as defined by the GNU Coding Standards.
include(GNUInstallDirs)  # It defines CMAKE_INSTALL_LIBDIR

# --- RPATH stuff ---
# Warning: RPATH settings must be defined before install(...) settings.
# Source : https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling

if(FORCE_SKIP_RPATH) # Build with no RPATH. Do we really need this option??
  set(CMAKE_SKIP_BUILD_RPATH TRUE)
else()
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
endif()

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# when building a binary package, it makes no sense to add this rpath
if(NOT FORCE_SKIP_RPATH)
  # the RPATH to be used when installing
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif(NOT FORCE_SKIP_RPATH)

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()

# init all common options for enabled components
set(common_options DOCUMENTATION TESTING UNSTABLE PYTHON_WRAPPER
  DOXYGEN_WARNINGS DOXY2SWIG)
foreach(opt ${common_options})
  init_to_default_option(${opt})
endforeach()

# ========= Documentation =========
if(WITH_DOCUMENTATION OR WITH_DOXY2SWIG OR WITH_DOXYGEN_WARNINGS OR WITH_GENERATION)
  set(USE_DOXYGEN TRUE)
endif()

# ----- Required dependencies (whatever Siconos components are) -----

# =========== Blas/Lapack ===========
# Find package stuff provided by cmake deals
# only with blas/lapack libraries.
# Since headers are also required for Siconos
# we use our own cmake "find package" stuff.
find_package(LAPACKDEV REQUIRED)

# =========== Boost ===========
# check https://cmake.org/cmake/help/latest/module/FindBoost.html?highlight=boost
if(WITH_CXX)
  find_package(Boost 1.61 REQUIRED)
endif()
#
# -- Python bindings --
if(WITH_PYTHON_WRAPPER)
  find_package(SWIG 3.0 REQUIRED)
  include(${SWIG_USE_FILE})
endif()

# ---- Extra 'developers only' options ----
# -- The options below should be used with caution
# and preferably by developers or advanced users --

# =========== OpenMP ==========
option(WITH_OPENMP "Use OpenMP" OFF)

# =========== use ccache if available ===========
option(WITH_CCACHE "Use ccache" OFF)
if(WITH_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif()
endif()

# =========== MPI ==========
option(WITH_MPI "Use MPI" OFF)
if(WITH_MPI)
  find_package(MPI REQUIRED)
  # https://cmake.org/cmake/help/v3.10/module/FindMPI.html
  if(MPI_CXX_FOUND)
    print_mpi_info(CXX)
    set(SICONOS_HAS_MPI TRUE) # for config.h
  endif()
  # if(MPI_Fortran_FOUND) # Do we need mpi fortran ?
  #   print_mpi_info(Fortran)
  # endif()
endif()

# ==== Python symlinks ===
# Useful during io installation, to use links for python scripts rather
# that copying files.
# FP : in which case do we need this ? If anyone knows, please document it ...
option(INSTALL_PYTHON_SYMLINKS "Install Python .py files as symlinks" OFF)

# For SiconosConfig.h
option(SICONOS_USE_MAP_FOR_HASH "Prefer std::map to std::unordered_map even if C++xy is enabled" ON)

