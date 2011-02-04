
SET_SOURCE_FILES_PROPERTIES(Numerics.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(Kernel.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(Numerics.i Kernel.i PROPERTIES SWIG_FLAGS "-DWall")

# do not skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

SWIG_ADD_MODULE(Numerics python Numerics.i)
SWIG_ADD_MODULE(Kernel python Kernel.i)



SWIG_LINK_LIBRARIES(Numerics ${PYTHON_LIBRARIES} 
  ${SiconosNumerics_LIBRARIES} 
  ${FORTRAN_LIBRARIES})

SWIG_LINK_LIBRARIES(Kernel ${PYTHON_LIBRARIES} 
  ${SiconosKernel_LIBRARIES} 
  ${SiconosNumerics_LIBRARIES} 
  ${FORTRAN_LIBRARIES} 
  ${Boost_LIBRARIES})

# get python dist-packages
# warning python may not know cmake_install_prefix
# installation in standard python dist package with:
# -DSTD_PYTHON_DIST_PACKAGES=TRUE 
#
# installation in specific place:
# -DPYTHON_DIST_PACKAGES=/my/specific/place

IF(STD_PYTHON_DIST_PACKAGES)
  EXECUTE_PROCESS ( COMMAND 
    python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" 
    OUTPUT_VARIABLE PYTHON_DIST_PACKAGES 
    OUTPUT_STRIP_TRAILING_WHITESPACE)
ELSE(STD_PYTHON_DIST_PACKAGES)
  IF(NOT PYTHON_DIST_PACKAGES)
    EXECUTE_PROCESS ( COMMAND 
      python -c "from distutils.sysconfig import get_python_lib; print get_python_lib(0,0,\"\")" 
      OUTPUT_VARIABLE PYTHON_DIST_PACKAGES 
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  ENDIF(NOT PYTHON_DIST_PACKAGES)
ENDIF(STD_PYTHON_DIST_PACKAGES)

# do not work :(
SET_TARGET_PROPERTIES(${SWIG_MODULE_Numerics_REAL_NAME} PROPERTIES 
  OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/swig/Siconos/numpy.i)

SET_TARGET_PROPERTIES(${SWIG_MODULE_Kernel_REAL_NAME} PROPERTIES OBJECT_DEPENDS 
  ${CMAKE_SOURCE_DIR}/src/swig/Siconos/numpy.i 
  ${CMAKE_SOURCE_DIR}/src/swig/Siconos/KernelTypes.i 
  ${CMAKE_SOURCE_DIR}/src/swig/Siconos/KernelRegistration.i)




MESSAGE(STATUS "PYTHON_DIST_PACKAGES : ${PYTHON_DIST_PACKAGES}")

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/Numerics.py 
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos)

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/Kernel.py 
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos)

CONFIGURE_FILE(__init__.py.in ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/__init__.py 
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos)

MESSAGE(STATUS "swig modules: ${SWIG_MODULE_Numerics_REAL_NAME},  ${SWIG_MODULE_Kernel_REAL_NAME}")

INSTALL(TARGETS ${SWIG_MODULE_Numerics_REAL_NAME} ${SWIG_MODULE_Kernel_REAL_NAME}
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos)

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/KernelTypes.i
	      ${CMAKE_CURRENT_SOURCE_DIR}/Kernel.i
	      ${CMAKE_CURRENT_SOURCE_DIR}/addons.hpp
  DESTINATION share/siconos-front-end/swig)

ENABLE_TESTING()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/tests/Rover1039.dat Rover1039.dat COPYONLY)


ADD_TEST(lcp ${CMAKE_BINARY_DIR}/src/swig/runtest.py ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_lcp.py)
ADD_TEST(Kernel ${CMAKE_BINARY_DIR}/src/swig/runtest.py ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_Kernel.py)
ADD_TEST(friction_contact ${CMAKE_BINARY_DIR}/src/swig/runtest.py ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_friction_contact.py)

SET_TESTS_PROPERTIES(lcp PROPERTIES FAIL_REGULAR_EXPRESSION "FAILURE;Exception;failed;ERROR;Assertion")
SET_TESTS_PROPERTIES(Kernel PROPERTIES FAIL_REGULAR_EXPRESSION "FAILURE;Exception;failed;ERROR;Assertion")
SET_TESTS_PROPERTIES(friction_contact PROPERTIES FAIL_REGULAR_EXPRESSION "FAILURE;Exception;failed;ERROR;Assertion")

SUBDIRS(tests)

