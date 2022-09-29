#== == == == == == == == == == == == == == == == == == == == == == == == == =
#Set main parameters for the Siconos cmake project
#
#== == == == == == == == == == == == == == == == == == == == == == == == == =
include(SiconosTools)

#== == == == == = Windows stuff... == == == == == =
include(WindowsSiconosSetup)

#-- -- -- -- - CMake project internal variables -- -- -- -- -
#Siconos current version
include(SiconosVersion)

include(FetchContent)

#-- Set include directories that are required by ALL components
#and only those !
#Other includes must be specified to individual targets only.
#Current binary dir, for generated headers.Only at build time !
include_directories($<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)

#Add extra logs about git references(branch name, commit number...)
#Useful for documentation and continuous integration
option(WITH_GIT "Consider sources are under GIT" OFF)

if(WITH_GIT) # User defined option, default = off
#Check if git is available
#and get last commit id(long and short).
#Saved in SOURCE_ABBREV_GIT_SHA1 and SOURCE_GIT_SHA1
#These vars are useful for tests logs and 'write_notes' macro.
  find_package(Git)
  if(GIT_FOUND)
    set(CTEST_GIT_COMMAND "${GIT_EXECUTABLE}" )     
    execute_process(COMMAND 
      ${GIT_EXECUTABLE} rev-parse --short HEAD
      OUTPUT_VARIABLE SOURCE_ABBREV_GIT_SHA1
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    
    execute_process(COMMAND 
      ${GIT_EXECUTABLE} rev-parse HEAD
      OUTPUT_VARIABLE SOURCE_GIT_SHA1
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endif()


#Save date / time into BUILD_TIMESTAMP var
string(TIMESTAMP BUILD_TIMESTAMP)

#-- -- PYTHON SETUP -- --
#Python interpreter is always required.
#We force Python3 !
#In addition, when WITH_PYTHON_WRAPPER is ON,
#we need python libraries and numpy.
if(${CMAKE_VERSION} VERSION_LESS "3.14")
#Our FindPython3 is just a copy of the one distributed
#with cmake = 3.14
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/extras)
endif()

#determine the order of preference between Apple - style and unix - style package components
#--> look for python framework when all other possibilities failed.
set(Python3_FIND_FRAMEWORK LAST)
if(WITH_PYTHON_WRAPPER)
  find_package(Python3 COMPONENTS Development Interpreter NumPy REQUIRED)
else()
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
endif()
#For backward compat...
set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
include(FindPythonModule)
find_python_module(packaging REQUIRED) # for siconos runtime
find_python_module(wheel REQUIRED) # for siconos runtime

get_filename_component(PYTHON_EXE_NAME ${PYTHON_EXECUTABLE} NAME)
if(WITH_PYTHON_WRAPPER OR WITH_DOCUMENTATION)
#-- - xml schema.Used in tests.-- -
  if(WITH_XML)
    set(SICONOS_XML_SCHEMA "${CMAKE_SOURCE_DIR}/kernel/swig/SiconosModelSchema-V3.7.xsd")
    if(NOT NO_RUNTIME_BUILD_DEP)
      find_python_module(lxml REQUIRED)
    endif()
  endif()
endif()

#-- - End of python conf -- -

#== == == == == = install setup == == == == == =

# --- A specific option to install everything in a completely isolated path (use-case: guix), lib, binaries AND python packages.
# In that case:
# - CMAKE_INSTALL_PREFIX is ignored and its content overwritten with ISOLATED_INSTALL value.
# - eveything is installed as usual (binaries in bin, libraries in lib ...) with ISOLATED_INSTALL as root dir
# - python packages are installed by pip with the option --prefix=${ISOLATED_INSTALL}.

if(ISOLATED_INSTALL)
  message(WARNING "You asked for a fully isolated installation in ${ISOLATED_INSTALL}.")
  if(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # if CMAKE_INSTALL_PREFIX has been explicitely set
    if(NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt) # if it's the first cmake run
      message(FATAL_ERROR "You can not set both ISOLATED_INSTALL and CMAKE_INSTALL_PREFIX. CMAKE_INSTALL_PREFIX value will be ignored.")
    endif()
    # Overwrite CMAKE_INSTALL_PREFIX with ISOLATED_INSTALL value
  endif()
  set(CMAKE_INSTALL_PREFIX ${ISOLATED_INSTALL} CACHE PATH "Install root directory." FORCE)
  set(sicopy_install_mode isolated CACHE STRING "Siconos Python packages nstall mode." FORCE)
else()
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(sicopy_install_mode standard CACHE STRING "Siconos Python packages install mode." FORCE)
    set(install_fixed ON CACHE INTERNAL "internal var: true after first run of cmake."  FORCE)
  # elseif(${install_fixed}) # second run (or more) of cmake, with standard install. We don't want any change.
  #   # this is just a double check to be sure to avoid install mess when users are playing with cache variables ...
  #   set(sicopy_install_mode standard CACHE STRING "Siconos Python packages install mode." FORCE)     
  else()
    set(sicopy_install_mode user CACHE STRING "Siconos Python packages install mode.")
  endif()
endif()
  
if(WITH_PYTHON_WRAPPER)
  if(siconos_python_install_mode)
    message(WARNING "You explicitely set siconos_python_install_mode to ${siconos_python_install_mode}! Be sure to know what you're doing and check the install paths printed after cmake run ... Or remove the whole build dir and start again!")
    set(sicopy_install_mode ${siconos_python_install_mode} CACHE STRING "Siconos Python packages install mode." FORCE)
  endif()
endif()

#Set directory used to save cmake config files
#required to use Siconos(e.g.to call find_package(siconos))
set(ConfigPackageLocation lib/cmake/siconos-${SICONOS_VERSION})

#Provides install directory variables as defined by the GNU Coding Standards.
include(GNUInstallDirs)  # It defines CMAKE_INSTALL_LIBDIR

#-- - RPATH stuff -- -
#Warning : RPATH settings must be defined before install(...) settings.
#Source : https: // gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling

if(FORCE_SKIP_RPATH) # Build with no RPATH. Do we really need this option??
  set(CMAKE_SKIP_BUILD_RPATH TRUE)
else()
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
endif()

#when building, don't use the install RPATH already
#(but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

#when building a binary package, it makes no sense to add this rpath
if(NOT FORCE_SKIP_RPATH)
#the RPATH to be used when installing
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif(NOT FORCE_SKIP_RPATH)

#don't add the automatically determined parts of the RPATH
#which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()

#init all common options for enabled components
set(common_options DOCUMENTATION TESTING UNSTABLE PYTHON_WRAPPER
  DOXYGEN_WARNINGS DOXY2SWIG)
foreach(opt ${common_options})
  init_to_default_option(${opt})
endforeach()

#== == == == = Documentation == == == == =
if(WITH_DOCUMENTATION OR WITH_DOXY2SWIG OR WITH_DOXYGEN_WARNINGS OR WITH_GENERATION)
  set(USE_DOXYGEN TRUE)
endif()

#-- -- - Required dependencies(whatever Siconos components are) -- -- -
#-- Python bindings --
if(WITH_PYTHON_WRAPPER)

  # Name of the generated Python package
  set(SICONOS_PYTHON_PACKAGE siconos CACHE INTERNAL "Name of the Siconos python package.")
  # --------------- Python install setup ---------------
  # Set path for siconos-python installation (SICONOS_PYTHON_INSTALL_DIR)
  # and get pip install options (PIP_INSTALL_OPTIONS).
  include(PythonInstallSetup)
  set_python_install_path()

  # -- swig stuff --
  include(swig_setup)
  
  #== == == Create(and setup) build / install target == == ==
  add_custom_target(python-install
    COMMAND ${PYTHON_EXECUTABLE} -m pip install -U ${CMAKE_BINARY_DIR}/wrap ${PIP_INSTALL_OPTIONS} -v 
    VERBATIM USES_TERMINAL
    COMMAND_EXPAND_LISTS
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "build/install siconos python package")
  #execute python - install when target install is called
  install(CODE "execute_process(COMMAND ${CMAKE_MAKE_PROGRAM} python-install WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")")

endif()

#== == == == == = Blas / Lapack == == == == == =
#Find package stuff provided by cmake deals
#only with blas / lapack libraries.
#Since headers are also required for Siconos
#we use our own cmake "find package" stuff.
find_package(LAPACKDEV REQUIRED)

#== == == == == = Boost == == == == == =
#check https: // cmake.org/cmake/help/latest/module/FindBoost.html?highlight=boost
if(WITH_CXX)

#From boost 1.71, something is wrong in cmake and boost support for multithread
#https: // github.com/boostorg/boost_install/issues/13
#https: // gitlab.kitware.com/cmake/cmake/issues/19714
#set(Boost_USE_MULTITHREADED ON)
  set(Boost_NO_BOOST_CMAKE 1)
  set(boost_min_version 1.61)
#Set the list of required boost components
  if(WITH_SERIALIZATION)
    list(APPEND boost_required_components serialization filesystem)
  endif()
  if(boost_required_components)
    set(boost_opts COMPONENTS ${boost_required_components})
  endif()

#Search boost...
  find_package(Boost ${boost_min_version} ${boost_opts} REQUIRED)

  if(WITH_SERIALIZATION)
    set(WITH_SYSTEM_BOOST_SERIALIZATION ON CACHE INTERNAL "Siconos uses boost serialization lib.")
  endif()
endif()

#-- Open Cascade Community Edition --
if(WITH_OCE)
  include(oce_setup)
endif()

#-- -- Extra 'developers only' options -- --
#-- The options below should be used with caution
#and preferably by developers or advanced users--

#== == == == == = OpenMP == == == == ==
option(WITH_OPENMP "Use OpenMP" OFF)

#== == == == == = use ccache if available == == == == == =
option(WITH_CCACHE "Use ccache" OFF)
if(WITH_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif()
endif()

#== == == == == = MPI == == == == ==
option(WITH_MPI "Use MPI" OFF)
if(WITH_MPI)
  find_package(MPI REQUIRED)
#https: // cmake.org/cmake/help/v3.10/module/FindMPI.html
  if(MPI_CXX_FOUND)
    print_mpi_info(CXX)
    set(SICONOS_HAS_MPI TRUE) # for config.h
  endif()
#if (MPI_Fortran_FOUND) #Do we need mpi fortran ?
#print_mpi_info(Fortran)
#endif()
endif()

#== == == == == = Tests env == == == == ==
if(WITH_TESTING)
#File used to print tests setup messages.
  set(TESTS_LOGFILE ${CMAKE_BINARY_DIR}/tests.log)

  if(CMAKE_BUILD_TYPE STREQUAL Debug)
#mlcp enum and dr_iso1(hairer) tests are quite long in debug mode...
    set(tests_timeout 700 CACHE INTERNAL "Limit time for tests (in seconds)")
  else()
    set(tests_timeout 200 CACHE INTERNAL "Limit time for tests (in seconds)")
  endif()
  if(HAVE_SICONOS_KERNEL)
    find_package(CPPUNIT REQUIRED)
#File used as main driver for cppunit tests
    set(SIMPLE_TEST_MAIN ${CMAKE_SOURCE_DIR}/kernel/tests-common/TestMain.cpp CACHE INTERNAL "")
  endif()
  if(WITH_PYTHON_WRAPPER)
    find_python_module(pytest REQUIRED)
    if(WITH_AGGRESSIVE_PYTHON_TESTS)
      set(pytest_opt "-s -v -pep8" CACHE INTERNAL "extra options for py.test")
    else()
      set(pytest_opt "-v" CACHE INTERNAL "extra options for py.test")
    endif()
  endif()
endif()

#== == Python symlinks == =
#Useful during io installation, to use links for python scripts rather
#that copying files.
#FP : in which case do we need this ? If anyone knows, please document it...
option(INSTALL_PYTHON_SYMLINKS "Install Python .py files as symlinks" OFF)

#For SiconosConfig.h
option(SICONOS_USE_MAP_FOR_HASH "Prefer std::map to std::unordered_map even if C++xy is enabled" ON)

#Check Siconos compilation with include - what - you - use
#See https: // github.com/include-what-you-use/include-what-you-use
#Set WITH_IWYU = path / to / iwyu binary file
if(WITH_IWYU)
#Clang is required for iwyu.This is a devel option, so we assume that
#you know what you are doing and that you use the same version of clang
#for both iwyu and Siconos.
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL Clang AND NOT CMAKE_CXX_COMPILER_ID STREQUAL AppleClang)
    message(FATAL_ERROR "You must compile Siconos with clang to use include-what-you-use.")
  endif()
  if(NOT CMAKE_C_COMPILER_ID STREQUAL Clang AND NOT CMAKE_C_COMPILER_ID STREQUAL AppleClang)
    message(FATAL_ERROR "You must compile Siconos with clang to use include-what-you-use.")
  endif()
  set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${WITH_IWYU};-Xiwyu;any;-Xiwyu;iwyu;-Xiwyu;args" CACHE INTERNAL "iwyu setup")
  set(CMAKE_C_INCLUDE_WHAT_YOU_USE "${WITH_IWYU};-Xiwyu;any;-Xiwyu;iwyu;-Xiwyu;args" CACHE INTERNAL "iwyu setup")
endif()
