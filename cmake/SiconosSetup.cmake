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

# Add extra logs about git references(branch name, commit number...)
# Useful for documentation and continuous integration
find_package(Git)
if(Git_FOUND)
  execute_process(COMMAND
    ${GIT_EXECUTABLE} -C ${CMAKE_SOURCE_DIR} rev-parse
    OUTPUT_VARIABLE NO_GIT
    ERROR_VARIABLE NO_GIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(NO_GIT STREQUAL "")
    set(WITH_GIT 1 CACHE INTERNAL "True if Git is available and Siconos sources are managed as a git repository.")
  else()
    set(WITH_GIT 0 CACHE INTERNAL "True if Git is available and Siconos sources are managed as a git repository.")
  endif()
else()
  set(WITH_GIT 0 CACHE INTERNAL "True if Git is available and Siconos sources are managed as a git repository.")
endif()

if(WITH_GIT)
# Get last commit id(long and short).
# Saved in SOURCE_ABBREV_GIT_SHA1 and SOURCE_GIT_SHA1
# These vars are useful for tests logs and 'write_notes' macro.
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


#Save date / time into BUILD_TIMESTAMP var
string(TIMESTAMP BUILD_TIMESTAMP)

#-- -- PYTHON SETUP -- --

# Determines the order of preference between Apple - style and unix - style package components
# --> look for python framework when all other possibilities failed.
include(FindPythonModule)
set(CMAKE_FIND_FRAMEWORK LAST)
if(WITH_PYTHON_WRAPPER)
  find_package(Python COMPONENTS Development Interpreter NumPy REQUIRED)
else()
  find_package(Python COMPONENTS Interpreter REQUIRED)
endif()
if(Python_VERSION_MAJOR VERSION_LESS 3)
  message(FATAL_ERROR "Python3 is required.")
endif()

find_python_module(packaging REQUIRED) # for siconos runtime
find_python_module(wheel REQUIRED) # for siconos runtime

if(WITH_PYTHON_WRAPPER OR WITH_DOCUMENTATION)
#-- - xml schema.Used in tests.-- -
  if(WITH_XML)
    set(SICONOS_XML_SCHEMA "${CMAKE_SOURCE_DIR}/kernel/swig/SiconosModelSchema-V3.7.xsd")
    if(NOT NO_RUNTIME_BUILD_DEP)
      find_python_module(lxml REQUIRED)
    endif()
  endif()
endif()
message(STATUS "End of Python configuration.\n")
message(STATUS "------------------------------------------------\n")

#-- - End of python conf -- -

#== == == == == = install setup == == == == == =
include(SiconosInstallSetup)
set_install_path()



#Set directory used to save cmake config files
#required to use Siconos(e.g.to call find_package(siconos))
set(SiconosConfigPackageLocation lib/cmake/siconos-${SICONOS_VERSION})

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

  # -- swig stuff --
  include(swig_setup)
  
  #== == == Create(and setup) build / install target == == ==
  add_custom_target(python-install
    COMMAND ${Python_EXECUTABLE} -m pip install -U ${CMAKE_BINARY_DIR}/wrap ${PIP_INSTALL_OPTIONS} -v 
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
  #https: // gitlab.kitware.com/cmake/cmake/issues/19714
  #set(Boost_USE_MULTITHREADED ON)
  #set(Boost_NO_BOOST_CMAKE 1)
  set(boost_min_version 1.71)
#Set the list of required boost components
  if(WITH_SERIALIZATION)
    list(APPEND boost_required_components serialization filesystem)
  endif()
  if(boost_required_components)
    set(boost_opts COMPONENTS ${boost_required_components})
  endif()

#Search boost...
  find_package(Boost ${boost_min_version} ${boost_opts} REQUIRED CONFIG)

  if(WITH_SERIALIZATION)
    set(WITH_SYSTEM_BOOST_SERIALIZATION ON CACHE INTERNAL "Siconos uses boost serialization lib.")
  endif()
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
