
SET_SOURCE_FILES_PROPERTIES(Numerics.i PROPERTIES CPLUSPLUS ON)
#SET_SOURCE_FILES_PROPERTIES(SiconosNumerics.i PROPERTIES SWIG_FLAGS "-includeall")

SWIG_ADD_MODULE(Numerics python Numerics.i)

SWIG_LINK_LIBRARIES(Numerics ${PYTHON_LIBRARIES} 
  ${SiconosKernel_LIBRARIES} 
  ${SiconosNumerics_LIBRARIES} 
  ${FORTRAN_LIBRARIES} 
  ${Boost_LIBRARIES})


# get python dist-packages
# warning python may not know cmake_install_prefix

execute_process ( COMMAND 
  python -c "from distutils.sysconfig import get_python_lib; print get_python_lib(0,0,\"${CMAKE_INSTALL_PREFIX}\")" 
  OUTPUT_VARIABLE PYTHON_DIST_PACKAGES 
  OUTPUT_STRIP_TRAILING_WHITESPACE)

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/Numerics.py 
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos)

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/__init__.py 
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos)

MESSAGE(STATUS "swig module: ${SWIG_MODULE_Numerics_REAL_NAME}")

INSTALL(TARGETS ${SWIG_MODULE_Numerics_REAL_NAME} 
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos)

ENABLE_TESTING()

ADD_TEST(lcp ${CMAKE_CURRENT_SOURCE_DIR}/tests/lcp.py)

CONFIGURE_FILE(__init__.py.in __init__.py)


#  PYTHONLIBS_FOUND       - have the Python libs been found
#  PYTHON_LIBRARIES       - path to the python library
#  PYTHON_INCLUDE_PATH    - path to where Python.h is found (deprecated)
#  PYTHON_INCLUDE_DIRS    - path to where Python.h is found
#  PYTHON_DEBUG_LIBRARIES - path to the debug library




MESSAGE(STATUS "PYTHON_DIST_PACKAGES : ${PYTHON_DIST_PACKAGES}")

SUBDIRS(tests)
