include(tools4tests)

if(WITH_${COMPONENT}_TESTING)

  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})

  # the main test driver
  SET(TEST_MAIN src/contactDetection/basicBroadphase/test/TestMain.cpp)

  BEGIN_TEST(src/contactDetection/basicBroadphase/test)
  NEW_TEST(testMultiBody MultiBodyTest.cpp)
  END_TEST()

  BEGIN_TEST(src/proposed)
  NEW_TEST(testContact testContact.cpp)
  END_TEST()

  IF(WITH_OCC)
    BEGIN_TEST(src/occ/test)
    NEW_TEST(testOcc OccTest.cpp)
    END_TEST()
  ENDIF()

  IF(WITH_MECHANISMS)
    MESSAGE("   ")
    MESSAGE("-------------------------------------- ************************ ")

    get_subdirectories(dirlist ${CMAKE_SOURCE_DIR}/examples/Mechanics/Mechanisms/)
    foreach(_dir ${dirlist})
      MESSAGE("${_dir}")
      set(bin_dir ${CMAKE_BINARY_DIR}/mechanics/swig/mechanics/mechanisms/test/${_dir})
      MESSAGE("mkdir ${bin_dir} ")
      file(MAKE_DIRECTORY ${bin_dir})
      file(GLOB EXAMPLES_FILES ${CMAKE_SOURCE_DIR}/examples/Mechanics/Mechanisms/${_dir}/bodydef.py
	${CMAKE_SOURCE_DIR}/examples/Mechanics/Mechanisms/${_dir}/mbtbLocalOptions.py)
      foreach(_F ${EXAMPLES_FILES})
	MESSAGE("${_F}")
	configure_file(${_F} ${bin_dir})
	file(READ ${_F} RESOURCES_FILE)
	string(REGEX REPLACE
	  "with3D=1"
	  "with3D=0 # remove X output for test#"
	  RESOURCES_FILE_MODIFIED ${RESOURCES_FILE})
	MESSAGE(${RESOURCES_FILE_MODIFIED})
	get_filename_component(_F_NAME ${_F} NAME)
	file(WRITE ${bin_dir}/${_F_NAME} ${RESOURCES_FILE_MODIFIED})
      endforeach()
      get_subdirectories(subdirlist ${CMAKE_SOURCE_DIR}/examples/Mechanics/Mechanisms/${_dir})
      
      foreach(_subdir ${subdirlist})
	if(NOT _subdir MATCHES siconos) 
	  MESSAGE("subdir ${_subdir}")
	  file(GLOB EXAMPLES_SUBDIR_FILES
	    ${CMAKE_SOURCE_DIR}/examples/Mechanics/Mechanisms/${_dir}/${_subdir}/*.*)
	  file(MAKE_DIRECTORY ${bin_dir}/${_subdir})
	  foreach(_F ${EXAMPLES_SUBDIR_FILES})
	    MESSAGE("${_F}")
	    configure_file(${_F} ${bin_dir}/${_subdir})
	  endforeach()
	endif()
      endforeach(_subdir)
      
      
      set(command_name "python ${CMAKE_BINARY_DIR}/scripts/siconos -P ${CMAKE_BINARY_DIR}/scripts/siconos-mechanisms.py .")
      MESSAGE("${command_name}")
      add_test(mechanisms_${_dir} ${CMAKE_COMMAND} -E chdir ${bin_dir} python ${CMAKE_BINARY_DIR}/scripts/siconos -P ${CMAKE_BINARY_DIR}/scripts/siconos-mechanisms.py . )
    endforeach(_dir)
 
  ENDIF()
  
endif()
