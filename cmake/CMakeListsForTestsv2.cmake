# -*- cmake -*-
# This is the test cmake configuration
# built from @CMAKE_SOURCE_DIR@/cmake/CMakeListsForTests.cmake.in 
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
  get_filename_component(BASE_BIN_DIR "@CMAKE_CURRENT_BINARY_DIR@" PATH)
  foreach(_C "externals" "numerics" "kernel" "mechanics" "control" "io")
    set(COMPONENT_BIN_DIR "${COMPONENT_BIN_DIR}\;${BASE_BIN_DIR}/${_C}")
  endforeach()
ELSE()
  SET(EXE_EXT)
ENDIF()

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
  IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
    ADD_EXECUTABLE(${_EXE} WIN32 ${${_EXE}_FSOURCES})
    SET_TARGET_PROPERTIES(${_EXE} PROPERTIES COMPILE_FLAGS "-DEMULATOR=\\\"wine\\\" -DWRAPPER=\\\"\\\"")
  ELSE()
    ADD_EXECUTABLE(${_EXE} ${${_EXE}_FSOURCES})
    SET_TARGET_PROPERTIES(${_EXE} PROPERTIES COMPILE_FLAGS "-DEMULATOR=\\\"\\\" -DWRAPPER=\\\"\\\"")
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
    ENDFOREACH(_TF ${${_EXE}_FSOURCES})
  ENDIF()

  IF(MSVC AND ${_EXE}_FORTRAN)
    SET_TARGET_PROPERTIES(${_EXE} PROPERTIES LINK_FLAGS "-Wl,--as-needed")
    SET(CMAKE_EXE_LINKER_FLAGS "")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "")
    SET(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "")
  ENDIF(MSVC AND ${_EXE}_FORTRAN)


  # -- link with current component and its dependencies --
  add_dependencies(${_EXE} @COMPONENT@)
  target_link_libraries(${_EXE} ${PRIVATE} @COMPONENT@)
  target_link_libraries(${_EXE} ${PRIVATE} ${@COMPONENT@_LINK_LIBRARIES})

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

  LIST(LENGTH ${_EXE}_DATA_LIST_${_CURRENT_TEST_DIRECTORY} count)
  MATH(EXPR count "${count}-1")
  FOREACH(i RANGE ${count})
    IF(count EQUAL -1) # there is no data file ...
      SET(_TEST_NAME ${_EXE})
      IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
        ADD_TEST(${_EXE} wine ${_EXE}${EXE_EXT})
      ELSE(CROSSCOMPILING_LINUX_TO_WINDOWS)
        ADD_TEST(${_EXE} ${_EXE}${EXE_EXT})
      ENDIF(CROSSCOMPILING_LINUX_TO_WINDOWS)
    ELSE(count EQUAL -1)
      LIST(GET ${_EXE}_DATA_LIST_${_CURRENT_TEST_DIRECTORY} ${i} _DATA_FILE)
      LIST(GET ${_EXE}_NAME_LIST_${_CURRENT_TEST_DIRECTORY} ${i} _TEST_NAME)
      IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
        ADD_TEST(${_TEST_NAME} wine ${_EXE}${EXE_EXT} ./data/${_DATA_FILE})
      ELSE(CROSSCOMPILING_LINUX_TO_WINDOWS)
        ADD_TEST(${_TEST_NAME} ${_EXE}${EXE_EXT} ./data/${_DATA_FILE})
      ENDIF(CROSSCOMPILING_LINUX_TO_WINDOWS)
    ENDIF(count EQUAL -1)

    SET_TESTS_PROPERTIES(${_TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION "FAILURE;Exception;failed;ERROR;test unsucceeded")

    if(CMAKE_SYSTEM_NAME MATCHES Windows)
      set(ENV_PPTY "Path=${COMPONENT_BIN_DIR}\;@ENV_PATH@")
    endif()

    SET(LOCAL_USE_SANITIZER "@USE_SANITIZER@")

    IF(LOCAL_USE_SANITIZER MATCHES "asan")
      set(ENV_PPTY "${ENV_PPTY};ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=1:$ENV{ASAN_OPTIONS}")
      set(ENV_PPTY "${ENV_PPTY};LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/Build/ci-scripts/asan-supp.txt:$ENV{LSAN_OPTIONS}")
    ENDIF(LOCAL_USE_SANITIZER MATCHES "asan")

    IF(ENV_PPTY)
      set_tests_properties(${_TEST_NAME} PROPERTIES ENVIRONMENT "${ENV_PPTY}")
    ENDIF(ENV_PPTY)


    IF(${_TEST_NAME}_PROPERTIES)
      SET_TESTS_PROPERTIES(${_TEST_NAME} PROPERTIES ${${_TEST_NAME}_PROPERTIES})
    ENDIF(${_TEST_NAME}_PROPERTIES)
    set_tests_properties(${_TEST_NAME} PROPERTIES TIMEOUT ${tests_timeout})

  ENDFOREACH(i RANGE ${count})

ENDFOREACH(_EXE ${_EXE_LIST_${_CURRENT_TEST_DIRECTORY}})
