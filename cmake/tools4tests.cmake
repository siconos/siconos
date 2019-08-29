# =========================================================
#
# Some cmake macros to deal with tests.
# =========================================================

#========================================
# Setup for tests in a given directory
#
# 
# Usage :
# begin_tests(<SOURCE_DIR> DEPS <dep1> <dep2>)
#
# SOURCE_DIR is relative to current component path
# DEPS is optional and used to add dependencies to <COMPONENT>-test target.
# Process :
# - set CURRENT_TEST_DIR (Parent scope) to <SOURCE_DIR>
# - copy all 'data' files found in <SOURCE_DIR> into <CMAKE_CURRENT_BINARY_DIR>/<SOURCE_DIR>
# - create a library from test-utils files (if any), named <COMPONENT>-test,
#   linked (PUBLIC) with <COMPONENT> lib.
#
function(begin_tests SOURCE_DIR)
  
  set(multiValueArgs DEPS)
  cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # Fix current test dir in parent scope.
  set(CURRENT_TEST_DIR ${SOURCE_DIR} PARENT_SCOPE)
  
  if(CROSSCOMPILING_LINUX_TO_WINDOWS)
    set(EMULATOR "wine")
    set(DRIVE_LETTER "Z:")
  else()
    set(EMULATOR "")
  endif()
  
  # find and copy data files in SOURCE_DIR (*.mat, *.dat and *.xml, ...) into test directory
  file(GLOB_RECURSE _DATA_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.mat
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.dat
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.hdf5
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.npz
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.DAT
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.INI
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/*.ref)

  # Copy data files from source to binary
  foreach(_F IN LISTS _DATA_FILES)
    #  log results into TESTS_LOGFILE
    file(APPEND ${TESTS_LOGFILE} "Found ref file : ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/${_F} \n")
    file(APPEND ${TESTS_LOGFILE} "Configuring file : ${CMAKE_CURRENT_BINARY_DIR}/${SOURCE_DIR}/${_F}\n")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/${_F} ${CMAKE_CURRENT_BINARY_DIR}/${SOURCE_DIR}/${_F} COPYONLY)
  endforeach()

  # If a directory *-utils is found, its content is used
  # to create/expand the <COMPONENT>-test library.
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}-utils)
    file(GLOB_RECURSE TEST_UTILS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}-utils/*.[ch])
    if(NOT TARGET ${COMPONENT}-test) 
      # Create the target ...
      add_library(${COMPONENT}-test SHARED ${TEST_UTILS_SOURCES})
    else()
      # or just append sources if it already exists.
      target_sources(${COMPONENT}-test PRIVATE ${TEST_UTILS_SOURCES})
    endif()

    # local includes, to build <component>-test only
    target_include_directories(${COMPONENT}-test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}-utils)
    target_link_libraries(${COMPONENT}-test PUBLIC ${COMPONENT})
    # All include dirs from component are taken into account in ${COMPONENT}-lib (and so propagated to tests)
    foreach(dir IN LISTS ${COMPONENT}_DIRS)
      target_include_directories(${COMPONENT}-test PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}>)
    endforeach()

    # Link with extra deps, if any
    foreach(libtarget IN LISTS TEST_DEPS)
      target_link_libraries(${COMPONENT}-test PUBLIC ${libtarget})
    endforeach()
  else()
    # If there is no ${COMPONENT}-test but some extra DEPS.
    unset(GLOBAL_TEST_DEPS)
    set(GLOBAL_TEST_DEPS ${TEST_DEPS} PARENT_SCOPE)
  endif()
  
endfunction()

macro(test_windows_setup test_name test_sources)

  # Use the proper linker for the proper language
  # fortran -> gfortran; {c,cpp} -> link.exe
  if(BUILD_AS_CPP)
    set(EXE_FORTRAN FALSE) # Windows only
    foreach(_file IN LISTS ${test_sources)
      if(${_file} MATCHES "[.]c$")
        set_source_files_properties(${_file} PROPERTIES LANGUAGE CXX)
      endif()
      if(${_file} MATCHES "[.]f$")
        set(EXE_FORTRAN TRUE) # Windows only
      endif()
    endforeach()
  endif()
  # Windows stuff ...
  if(MSVC AND EXE_FORTRAN)
    target_link_options(${test_name} "LINKER:--as-needed")
    if(NOT V2)
      target_link_options(${_EXE} "LINKER:--subsystem,console")
      set(CMAKE_CREATE_CONSOLE_EXE "")
      # Note FP : this variable does not exist anymore in CMake.
      # What was it used for??
    endif()
    # Linker flags to be used to create executables.
    SET(CMAKE_EXE_LINKER_FLAGS "")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "")
    SET(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "")
  endif()

endmacro()

# ========================================
# Add a test into the ctest suite.
#
# Usage :
# new_test(NAME <name> SOURCES <sources list> DEPS <dependencies list> DATA <data files list)
#
# required : SOURCES
# others are optional.
#
# Process:
# - Copy data files to binary dir
# - Create executable from sources list
# - link (PRIVATE) executable with <COMPONENT> lib
# - add <COMPONENT> includes to executable (PRIVATE)
# - link (PRIVATE) executable with all libs in DEPS
# - link (PRIVATE) executable with <COMPONENT>-test (if it exists)
# - add a test (ctest) named <name>. If NAME is not set, use name of first source file (without ext).
# ========================================
function(new_test)
  set(oneValueArgs NAME HDF5)
  set(multiValueArgs SOURCES DATA DEPS)
  cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # -- set test name --
  # If not set in input args, we choose first source file name without extension.
  if(NOT TEST_NAME)
    list(GET TEST_SOURCES 0 TEST_NAME)
    get_filename_component(TEST_NAME ${TEST_NAME} NAME_WE)
  endif()

  # -- log --
  file(APPEND ${TESTS_LOGFILE} "Add test ${CURRENT_TEST_DIR}/${TEST_NAME} \n")

  # -- Check if data files are required --
  # If so, copy the corresponding file into the working directory
  if(TEST_DATA)
    foreach(_datafile IN LISTS TEST_DATA)
      file(APPEND ${TESTS_LOGFILE} "Copy data file : ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TEST_DIR}/${_datafile} \n")
      configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TEST_DIR}/${_datafile}
        ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_TEST_DIR}/${_datafile} COPYONLY
        )
    endforeach()
  endif()
 
  # ---- Build test executable and set its properties ----
  # Req : test sources
  # !!! WARNING !!!
  # Consider known limitations for Windows written
  # here https://cmake.org/cmake/help/v3.13/prop_dir/COMPILE_DEFINITIONS.html
  
  foreach(source_file IN LISTS TEST_SOURCES)
    # If source file path is not absolute, prepend current dir
    if(NOT IS_ABSOLUTE ${source_file})
      set(source_file ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TEST_DIR}/${source_file})
    endif()
    list(APPEND TEST_SOURCES_ABSPATH ${source_file})
  endforeach()
  
  if(CROSSCOMPILING_LINUX_TO_WINDOWS)
    add_executable(${TEST_NAME} WIN32 ${TEST_SOURCES_ABSPATH})
  else()
    add_executable(${TEST_NAME} ${TEST_SOURCES_ABSPATH})
  endif()
  target_compile_definitions(${TEST_NAME} PRIVATE "EMULATOR=\"${EMULATOR}\";WRAPPER=\"\"")
  if(TEST_HDF5)
    target_compile_definitions(${TEST_NAME} PRIVATE "TEST_HDF5")
  endif()
  
  
  # Set path where exe should be generated
  set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_TEST_DIR}/)

  test_windows_setup(${TEST_NAME} ${TEST_SOURCES})

  # -- link with current component and its dependencies --
  # --> add include and links
  target_link_libraries(${TEST_NAME} PRIVATE ${COMPONENT})
  foreach(dir IN LISTS ${COMPONENT}_DIRS)
    target_include_directories(${TEST_NAME} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}>)
  endforeach()

  # Link with extra deps (set for each test in new_test call)
  foreach(libtarget IN LISTS TEST_DEPS)
    target_link_libraries(${TEST_NAME} PRIVATE ${libtarget})
  endforeach()
  # ... or with GLOBAL_TEST_DEPS variable. GLOBAL_TEST_DEPS is
  # set during call to begin_test and useful only
  # when some dependencies are required by all tests
  # and if there is no <component>-test lib.
  foreach(libtarget IN LISTS GLOBAL_TEST_DEPS) 
    target_link_libraries(${TEST_NAME} PRIVATE ${libtarget})
  endforeach()

  # -- link with (optional) component-test lib.
  # At the time, only numerics used such a lib.
  if(TARGET ${COMPONENT}-test)
    target_link_libraries(${TEST_NAME} PRIVATE ${COMPONENT}-test)
  endif()
  if(WITH_CXX)
    set_target_properties(${TEST_NAME} PROPERTIES LINKER_LANGUAGE CXX)
  endif()
  
  if (LDLIBPATH)
    set_target_properties(${TEST_NAME} PROPERTIES ENVIRONMENT "${LDLIBPATH}")
  endif()

  # ---- The exe is ready, let's turn to test(s) ----
  # -- set test command --
  set(command ${TEST_NAME})
  if(CMAKE_SYSTEM_NAME MATCHES Windows)
    set(command ${command}".exe")
  endif()

  # Add the test in the pipeline
  add_test(NAME ${TEST_NAME} COMMAND ${command} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_TEST_DIR})
  set_siconos_test_properties(NAME ${TEST_NAME})
  
endfunction()

function(set_siconos_test_properties)
  set(oneValueArgs NAME)
  set(multiValueArgs PROPERTIES)
  cmake_parse_arguments(TEST "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # if the test output matches one of specified expressions below, the test will fail
  #set_tests_properties(${TEST_NAME} PROPERTIES
  #  FAIL_REGULAR_EXPRESSION "FAILURE;Exception;failed;ERROR;test unsucceeded")
  if (LDLIBPATH)
    set_tests_properties(${TEST_NAME} PROPERTIES ENVIRONMENT "${LDLIBPATH}")
  endif()
  
  set(LOCAL_USE_SANITIZER "@USE_SANITIZER@")
  if(LOCAL_USE_SANITIZER MATCHES "asan")
    set_property(TEST ${TEST_NAME} APPEND PROPERTY ENVIRONMENT ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=1:$ENV{ASAN_OPTIONS})
    set_property(TEST ${TEST_NAME} APPEND PROPERTY ENVIRONMENT LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/Build/ci-scripts/asan-supp.txt:$ENV{LSAN_OPTIONS})
  endif()
  
  # Extra tests properties (set in <component>_tests.cmake)
  if(TEST_PROPERTIES)
    set_tests_properties(${TEST_NAME} PROPERTIES "${TEST_PROPERTIES}")
  endif()
  
  # Test timer
  set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT ${tests_timeout})
  
endfunction()

# ================================================
# Build a test for a set of data and/or solvers.
#
# It usually needs :
# - a driver (.c.in file)
# - a formulation name (eg fc2d, LCP ...)
# - a 'collection' name, something to identify the set of data files
# - a list of sources files (in addition to the .c generated from the driver)
#
# Usage :
#
# Process :
# - generate <COLLECTION><SUFFIX>.c file from <DRIVER>
#   variables required (@XX@ in .c.in) : PROBLEM_FORMULATION, the formulation and PROBLEM_COLLECTION, the collection.
# - create the test named <NAME>_<COLLECTION><SUFFIX> from sources and data set.
function(new_tests_collection)
  set(oneValueArgs DRIVER NAME COLLECTION SUFFIX FORMULATION HDF5)
  set(multiValueArgs DATASET EXTRA_SOURCES)
  cmake_parse_arguments(PROBLEM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # set(TEST_SOLVER SICONOS_${PROBLEM_NAME}_${PROBLEM_SOLVER}) # Required for configure below!
  # This value is replaced in solver call in .c file.
  # Generate source file
  set(generated_driver_name ${PROBLEM_FORMULATION}_${PROBLEM_COLLECTION}${PROBLEM_SUFFIX}.c)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TEST_DIR}/${PROBLEM_DRIVER}
    ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_TEST_DIR}/${generated_driver_name})

  set(TEST_NAME_PREFIX ${PROBLEM_FORMULATION})
  set(TEST_COLLECTION ${PROBLEM_COLLECTION})
  new_test(
    SOURCES ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_TEST_DIR}/${generated_driver_name} ${PROBLEM_EXTRA_SOURCES}
    NAME ${PROBLEM_FORMULATION}_${PROBLEM_COLLECTION}${PROBLEM_SUFFIX}
    DATA_SET "${PROBLEM_DATASET}"
    HDF5 ${PROBLEM_HDF5}
    )
 
endfunction()

# Specialisation of tests_collection to lcp formulation.
function(new_lcp_tests_collection)
  set(options)
  set(oneValueArgs COLLECTION SUFFIX)
  set(multiValueArgs DATASET EXTRA_SOURCES)
  cmake_parse_arguments(PROBLEM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  new_tests_collection(
    DRIVER lcp_test_collection.c.in FORMULATION lcp COLLECTION ${PROBLEM_COLLECTION}
    EXTRA_SOURCES ${PROBLEM_EXTRA_SOURCES})
  
endfunction()

# Build plugins required for python tests
function(build_plugin)
  set(oneValueArgs FILE)
  set(multiValueArgs DEPS)
  cmake_parse_arguments(plug "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  get_filename_component(plug_name ${plug_FILE} NAME_WE)

  add_library(${plug_name} MODULE ${plug_FILE})
  target_include_directories(${plug_name} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tests/plugins/>)

  # link with extra deps.
  foreach(dep IN LISTS plug_DEPS)
    target_link_libraries(${plug_name} PRIVATE ${dep})
  endforeach()

  set_property(TARGET ${plug_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/tests)
  set_target_properties(${plug_name} PROPERTIES PREFIX "")
  if(NOT WITH_CXX)
    set_source_files_properties(${plug_FILE} PROPERTIES LANGUAGE C)
  endif(NOT WITH_CXX)
  if(WITH_MPI)
    # Note FP : temporary fix, to deal with PRIVATE deps of some components.
    # This will be reviewed later.
    target_include_directories(${plug_name} PRIVATE ${MPI_C_INCLUDE_DIRS})
  endif()
endfunction()

macro(set_ldlibpath)
  # In certain cases, ex. no rpath, or running tests with plugins,
  # libraries cannot be found at link or test time, so we add the
  # LD_LIBRARY_PATH variable.
  SET(LDLIBPATH)
  if (CMAKE_SKIP_RPATH)
    foreach(_C IN LISTS COMPONENTS)
      LIST(APPEND LDLIBPATH "${CMAKE_BINARY_DIR}/${_C}")
    ENDFOREACH()
  else()
    # Otherwise, still need the path to current component dir for tests
    # that load plugins.
    LIST(APPEND LDLIBPATH "${CMAKE_BINARY_DIR}/${COMPONENT}")
  endif()
  LIST(APPEND LDLIBPATH "${CMAKE_BINARY_DIR}/wrap/siconos/tests")
  if (NOT CMAKE_SYSTEM_NAME MATCHES WINDOWS)
    STRING(REPLACE ";" ":" LDLIBPATH "${LDLIBPATH}")
  endif()
  if (CMAKE_SYSTEM_NAME MATCHES APPLE)
    if ($ENV{DYLD_LIBRARY_PATH})
      set(LDLIBPATH "${LDLIBPATH}:$ENV{DYLD_LIBRARY_PATH}")
    endif()
    SET(LDLIBPATH "DYLD_LIBRARY_PATH=${LDLIBPATH}")
  else()
    if (CMAKE_SYSTEM_NAME MATCHES WINDOWS)
      SET(LDLIBPATH "Path=${LDLIBPATH};$ENV{Path}")
    else()
      if ($ENV{LD_LIBRARY_PATH})
        set(LDLIBPATH "${LDLIBPATH}:$ENV{LD_LIBRARY_PATH}")
      endif()
      SET(LDLIBPATH "LD_LIBRARY_PATH=${LDLIBPATH}")
    endif()
  endif()
endmacro()

# Declaration of a siconos test based on python bindings
function(add_python_test test_name test_file)
  add_test(${test_name} ${PYTHON_EXECUTABLE} ${TESTS_RUNNER} "${pytest_opt}" ${DRIVE_LETTER}${test_file})
  set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/tests)
  set_tests_properties(${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "FAILURE;Exception;[^x]failed;ERROR;Assertion")
  set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${CMAKE_BINARY_DIR}/wrap")
  if(LDLIBPATH)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "${LDLIBPATH}")
  endif()
endfunction()


# ----------------------------------------
# Prepare python tests for the current
# component
# Usage:
#   build_python_tests(path_to_tests)
#
# path_to_tests is relative to the current source dir.
# Most of the time, path_to_tests = 'tests'.
# For instance, in mechanics, tests are called in CMakeLists.txt
# in swig, current source dir is thus source_dir/mechanics/swig
# and source_dir/mechanics/swig/tests contains all the python files for tests.
#
# This routine copy the directory of tests to binary dir to allow 'py.test' run in the build.
# 
# binary dir will then look like :
# wrap/siconos
# wrap/siconos/mechanics
# wrap/siconos/mechanics/tests
#
# and running py.tests in wrap dir will end up with a run of
# all mechanics tests.
#
# Usage :
# build_python_tests(DEPS <list of targets> EXCLUDE <list of python files>)
#
#  * DEPS : list of targets that must be linked with c/c++ plugins used by python tests
#  * EXCLUDE : list of python files (path relative to current source dir) that
#    must not be run as tests.
# Both DEPS and EXCLUDE are optional.
function(build_python_tests)
  set(multiValueArgs DEPS EXCLUDE)
  cmake_parse_arguments(test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  
  # build plugins, if any
  # Note : all built libraries are saved in SICONOS_SWIG_ROOT_DIR/plugins
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/plugins)
    file(GLOB plugfiles ${CMAKE_CURRENT_SOURCE_DIR}/tests/plugins/*.cpp)
    foreach(plug ${plugfiles})
      build_plugin(FILE ${plug} DEPS ${test_DEPS})
    endforeach()
  endif()
    
  # copy test dir to binary dir (inside siconos package)
  # ---> allows py.test run in binary dir
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
    file(GLOB data4tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/tests/data
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/*)
    foreach(datafile ${data4tests})
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/data/${datafile}
	${SICONOS_SWIG_ROOT_DIR}/tests/data/${datafile} COPYONLY)
    endforeach()
  endif()
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD)
    file(GLOB data4tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD/*)
    foreach(datafile ${data4tests})
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD/${datafile}
	${SICONOS_SWIG_ROOT_DIR}/tests/CAD/${datafile} COPYONLY)
    endforeach()
  endif()
  if(CROSSCOMPILING_LINUX_TO_WINDOWS)
    set(EMULATOR "wine")
    set(DRIVE_LETTER "Z:")
  else()
    set(EMULATOR "")
  endif()

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    file(GLOB testfiles ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.py)
    foreach(excluded_test ${test_EXCLUDE})
      list(REMOVE_ITEM testfiles ${CMAKE_CURRENT_SOURCE_DIR}/${excluded_test})
    endforeach()
    foreach(file ${testfiles})
      get_filename_component(testname ${file} NAME_WE)
      get_filename_component(exename ${file} NAME)
      # Each file is copy to siconos/tests.
      # Maybe we can create a 'tests' dir for each subpackage?
      # --> Easier to deal with plugins and data if only one package
      configure_file(${file} ${SICONOS_SWIG_ROOT_DIR}/tests COPYONLY)
      set(name "python_${testname}")
      set(exename ${SICONOS_SWIG_ROOT_DIR}/tests/${exename})
      # set_ldlibpath()
      add_python_test(${name} ${exename})
    endforeach()
  endif()
endfunction()
