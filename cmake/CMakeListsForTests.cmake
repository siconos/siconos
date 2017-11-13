# -*- cmake -*-
# This is the test cmake configuration
# built from @CMAKE_SOURCE_DIR@/cmake/CMakeListsForTests.cmake
SET(SOURCE_DIR @CMAKE_CURRENT_SOURCE_DIR@/@_CURRENT_TEST_DIRECTORY@)

# Search for reference files and copy them to binary dir
FILE(GLOB TESTS_REF ${SOURCE_DIR}/*.ref)
FOREACH(_F ${TESTS_REF})
  GET_FILENAME_COMPONENT(TEST_REF ${_F} NAME)
  file(APPEND ${TESTS_LOGFILE} "Found ref file : ${_F} \n")
  file(APPEND ${TESTS_LOGFILE} "Configuring file : ${CMAKE_CURRENT_BINARY_DIR}/${TEST_REF} \n")
  CONFIGURE_FILE(${_F} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_REF}  COPYONLY)
ENDFOREACH(_F ${TESTS_REF})

IF(CMAKE_SYSTEM_NAME MATCHES Windows)
  SET(EXE_EXT ".exe")
ELSE()
  SET(EXE_EXT)
ENDIF()

if(CMAKE_SYSTEM_NAME MATCHES Windows)
  set(COMPONENT_BIN_DIR)
  get_filename_component(BASE_BIN_DIR "@CMAKE_CURRENT_BINARY_DIR@" DIRECTORY)
  foreach(_C "externals" "numerics" "kernel" "mechanics" "control" "io")
    set(COMPONENT_BIN_DIR "${COMPONENT_BIN_DIR}\;${BASE_BIN_DIR}/${_C}")
  endforeach()
  SET(COMPONENT_PATH "${COMPONENT_BIN_DIR}\;@ENV_PATH@")
endif()

# For some environment variables (LD_LIBRARY_PATH, DYLD_LIBRARY_PATH, Path)
# Needed for plugins --xhub
set(LIBFORTests @CMAKE_CURRENT_BINARY_DIR@)

FOREACH(_EXE ${_EXE_LIST_${_CURRENT_TEST_DIRECTORY}})
  
  file(APPEND ${TESTS_LOGFILE} "Adding test suite ${_CURRENT_TEST_DIRECTORY}/${_EXE} \n")

  # Set include directories for the current test :
  # all dirs from main project
  FOREACH(_D ${SICONOS_INCLUDE_DIRECTORIES})
    INCLUDE_DIRECTORIES(${_D})
  ENDFOREACH()

  # Set extra include directories for the current test :
  # 'test-utils' dir
  FOREACH(_D ${${_CURRENT_TEST_DIRECTORY}_INCLUDE_DIRECTORIES})
    INCLUDE_DIRECTORIES(${_D})
  ENDFOREACH(_D ${${_CURRENT_TEST_DIRECTORY}_INCLUDE_DIRECTORIES})

  # -- Build test executable and set its properties --- 
  # a wrapper around test
  IF(TEST_WRAP)
    ADD_EXECUTABLE(${_EXE}.ldwrap ${_EXE}.ldwrap.c)
  ENDIF(TEST_WRAP)
  
  IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
    ADD_EXECUTABLE(${_EXE} WIN32 ${${_EXE}_FSOURCES})
    SET_TARGET_PROPERTIES(${_EXE} PROPERTIES COMPILE_FLAGS "-DEMULATOR=\\\"wine\\\" -DWRAPPER=\\\"\\\"")
  ELSE()
    ADD_EXECUTABLE(${_EXE} ${${_EXE}_FSOURCES})
    IF(TEST_WRAP)
      SET_TARGET_PROPERTIES(${_EXE} PROPERTIES COMPILE_FLAGS "-DEMULATOR=\\\"\\\" -DWRAPPER=\\\".ldwrap\\\"")
    ELSE(TEST_WRAP)
      SET_TARGET_PROPERTIES(${_EXE} PROPERTIES COMPILE_FLAGS "-DEMULATOR=\\\"\\\" -DWRAPPER=\\\"\\\"")
    ENDIF(TEST_WRAP)
  ENDIF()

  #Use the proper linker for the proper language
  # fortran -> gfortran; {c,cpp} -> link.exe
  IF(BUILD_AS_CPP)
    SET(${_EXE}_FORTRAN FALSE)
    FOREACH(_TF ${${_EXE}_FSOURCES})
      IF(${_TF} MATCHES "[.]c$")
        set_source_files_properties(${_TF} PROPERTIES LANGUAGE CXX)
      ENDIF()
      IF(${_TF} MATCHES "[.]f$")
        SET(${_EXE}_FORTRAN TRUE)
      ENDIF()
    ENDFOREACH()
  ENDIF(BUILD_AS_CPP)
  
  IF(MSVC AND ${_EXE}_FORTRAN)
    SET_TARGET_PROPERTIES(${_EXE} PROPERTIES LINK_FLAGS "-Wl,--as-needed -Wl,--subsystem,console")
    SET(CMAKE_EXE_LINKER_FLAGS "")
    SET(CMAKE_CREATE_CONSOLE_EXE "")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "")
    SET(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "")
  ENDIF(MSVC AND ${_EXE}_FORTRAN)

  # -- link with current component and its dependencies --
  add_dependencies(${_EXE} @COMPONENT@)
  target_link_libraries(${_EXE} ${PRIVATE} @COMPONENT@)
  target_link_libraries(${_EXE} ${PRIVATE} ${@COMPONENT@_LINK_LIBRARIES})

  # if no rpath, tests cannot find the libraries and build fails, so
  # we add the LD_LIBRARY_PATH variable
  if (CMAKE_SKIP_RPATH)
    SET(LDLIBPATH "")
    FOREACH(_C ${COMPONENTS})
      LIST(APPEND LDLIBPATH "${CMAKE_BINARY_DIR}/${_C}")
    ENDFOREACH()
    STRING(REPLACE ";" ":" LDLIBPATH "${LDLIBPATH}")
    SET(LDLIBPATH "LD_LIBRARY_PATH=${LDLIBPATH}")
  endif()

  set(COMPONENT_TEST_LIB_ @COMPONENT_TEST_LIB@)
  if(COMPONENT_TEST_LIB_)
    add_dependencies(${_EXE} @COMPONENT_TEST_LIB@)
    target_link_libraries(${_EXE} ${PRIVATE} @COMPONENT_TEST_LIB@)
  endif()

  # Link and include for tests libraries (e.g. cppunit ...)
  FOREACH(_L ${TEST_LIBS})
    TARGET_LINK_LIBRARIES(${_EXE} ${PRIVATE} ${_L})
  ENDFOREACH()
  FOREACH(_D ${TEST_INCLUDE_DIR})
    include_directories(${_D})
  ENDFOREACH()

  IF(CPPUNIT_FOUND)
    # each test in the test suite becomes a cmake test


    IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
      SET(EMULATOR "wine")
    ELSE()
      SET(EMULATOR "")
    ENDIF()

    # --- Create a target for the current test ---
    IF(CMAKE_SYSTEM_NAME MATCHES Windows)
      ADD_CUSTOM_COMMAND(TARGET ${_EXE}
        POST_BUILD
        COMMAND env "PATH=\"${COMPONENT_PATH};\"" ${EMULATOR} ARGS ${CMAKE_CURRENT_BINARY_DIR}/${_EXE}${EXE_EXT}
        --cdash-prepare ${CMAKE_CURRENT_BINARY_DIR}/${_EXE}${EXE_EXT} > ${CMAKE_CURRENT_BINARY_DIR}/${_EXE}.cmake
        COMMENT "Generating ${_EXE}.cmake")
    ELSE()
      ADD_CUSTOM_COMMAND(TARGET ${_EXE}
        POST_BUILD
        COMMAND env 
        ARGS "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/Build/ci-scripts/asan-supp.txt:$ENV{LSAN_OPTIONS}"
        "ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=1:$ENV{ASAN_OPTIONS}"
        ${LDLIBPATH}
        ${CMAKE_CURRENT_BINARY_DIR}/${_EXE}${EXE_EXT}
        --cdash-prepare ${CMAKE_CURRENT_BINARY_DIR}/${_EXE}${EXE_EXT} > ${CMAKE_CURRENT_BINARY_DIR}/${_EXE}.cmake
        COMMENT "Generating ${_EXE}.cmake")
    ENDIF()
    
    # -- Generate a cmake macro to create a test, write into SiconosTestConfig.cmake --
    FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/SiconosTestConfig.cmake "# siconos test config file\n")
    FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/SiconosTestConfig.cmake "MACRO(ADD_CPPUNIT_TEST)\n")
    FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/SiconosTestConfig.cmake "  ADD_TEST(\${ARGV})\n")
    FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/SiconosTestConfig.cmake "  SET(_EXE \${ARGV0})\n")

    IF(CMAKE_SYSTEM_NAME MATCHES Windows)
      SET(ENV_PPTY "Path=${COMPONENT_PATH}")
    ELSEIF(NOT CMAKE_SYSTEM_NAME MATCHES APPLE)
      SET(ENV_PPTY "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${LIBFORTests}")
    ENDIF()

    IF(USE_SANITIZER MATCHES "asan")
      SET(ENV_PPTY "${ENV_PPTY};ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=1:$ENV{ASAN_OPTIONS}")
      SET(ENV_PPTY "${ENV_PPTY};LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/Build/ci-scripts/asan-supp.txt:$ENV{LSAN_OPTIONS}")
    ENDIF(USE_SANITIZER MATCHES "asan")

    IF(ENV_PPTY)
      FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/SiconosTestConfig.cmake "  SET_TESTS_PROPERTIES(\${_EXE} PROPERTIES ENVIRONMENT \"${ENV_PPTY}\")\n")
    ENDIF(ENV_PPTY)


    FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/SiconosTestConfig.cmake "ENDMACRO(ADD_CPPUNIT_TEST)\n")

    # -- testname.cmake will be included when ctest will be run.
    SET_DIRECTORY_PROPERTIES(PROPERTIES TEST_INCLUDE_FILE
      "${CMAKE_CURRENT_BINARY_DIR}/${_EXE}.cmake")
    
  ELSE()
    IF(TEST_WRAP)
      IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
        ADD_TEST(${_EXE} wine ${_EXE}.ldwrap${EXE_EXT})
      ELSE()
        ADD_TEST(${_EXE} ${_EXE}.ldwrap${EXE_EXT})
      ENDIF()
    ELSE()
      IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
        ADD_TEST(${_EXE} wine ${_EXE}${EXE_EXT})
      ELSE()
        ADD_TEST(${_EXE} ${_EXE}${EXE_EXT})
      ENDIF()
    ENDIF()

    SET_TESTS_PROPERTIES(${_EXE} PROPERTIES FAIL_REGULAR_EXPRESSION "FAILURE;Exception;failed;ERROR;test unsucceeded")

    # do not enable CXX linker if we don't build cxx code ...
    SET(WITH_CXX_ @WITH_CXX@)
    IF(WITH_CXX_)
      SET_TARGET_PROPERTIES(${_EXE} PROPERTIES LINKER_LANGUAGE CXX)
    ENDIF(WITH_CXX_)

    if(CMAKE_SYSTEM_NAME MATCHES Windows)
      set_tests_properties(${_EXE} PROPERTIES ENVIRONMENT "Path=${COMPONENT_PATH}")
    ELSEIF(NOT CMAKE_SYSTEM_NAME MATCHES APPLE)
      SET(ENV_PPTY "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${LIBFORTests}")
    endif()

    IF(USE_SANITIZER MATCHES "asan")
      set_property(TEST ${_EXE} APPEND PROPERTY ENVIRONMENT ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=1:$ENV{ASAN_OPTIONS})
      set_property(TEST ${_EXE} APPEND PROPERTY ENVIRONMENT LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/Build/ci-scripts/asan-supp.txt:$ENV{LSAN_OPTIONS})
    ENDIF(USE_SANITIZER MATCHES "asan")

    IF(${_EXE}_PROPERTIES)
      SET_TESTS_PROPERTIES(${_EXE} PROPERTIES ${${_EXE}_PROPERTIES})
    ENDIF(${_EXE}_PROPERTIES)

    if(${_EXE}_TIMEOUT)
      set_tests_properties(${_EXE} PROPERTIES TIMEOUT ${${_EXE}_TIMEOUT})
    else()
      set_tests_properties(${_EXE} PROPERTIES TIMEOUT ${tests_timeout})
    endif()

  ENDIF()

ENDFOREACH()
