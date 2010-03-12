
SET_SOURCE_FILES_PROPERTIES(Numerics.i PROPERTIES CPLUSPLUS ON)
#SET_SOURCE_FILES_PROPERTIES(SiconosNumerics.i PROPERTIES SWIG_FLAGS "-includeall")

SWIG_ADD_MODULE(Numerics python Numerics.i)



SWIG_LINK_LIBRARIES(Numerics ${PYTHON_LIBRARIES} 
  ${SiconosKernel_LIBRARIES} 
  ${SiconosNumerics_LIBRARIES} 
  ${FORTRAN_LIBRARIES} 
  ${Boost_LIBRARIES})

# this fails if python libraries are not under CMAKE_PREFIX!
# so please consider distutils...
#INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/SiconosNumerics.py DESTINATION lib/${_PPN}/dist-packages)

MESSAGE(STATUS "swig module: ${SWIG_MODULE_SiconosNumerics_REAL_NAME}")
#INSTALL(TARGETS ${SWIG_MODULE_SiconosNumerics_REAL_NAME}  LIBRARY DESTINATION lib/${_PPN}/dist-packages)

ENABLE_TESTING()

ADD_TEST(lcp ${CMAKE_CURRENT_SOURCE_DIR}/tests/lcp.py)

CONFIGURE_FILE(__init__.py.in __init__.py)

SUBDIRS(tests)
