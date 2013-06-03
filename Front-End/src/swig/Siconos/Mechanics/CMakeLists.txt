MACRO(MECHANICS_ADD_MODULE _mod _dest)
  
  SET_SOURCE_FILES_PROPERTIES(${_mod}.i PROPERTIES CPLUSPLUS ON)

  SET_SOURCE_FILES_PROPERTIES(${_mod}.i 
    PROPERTIES SWIG_FLAGS "-dirprot;-I${CMAKE_SOURCE_DIR}/src/swig/Siconos;-I${CMAKE_BINARY_DIR}/src/swig/Siconos;-Wall;${SWIG_DEFS}")

  SWIG_ADD_MODULE(${_mod} python ${_mod}.i)
  
  SET(${_mod}_generated_file_fullname ${swig_generated_file_fullname})

  MESSAGE(STATUS "${_mod} generated file ${${_mod}_generated_file_fullname}")

  ADD_DEPENDENCIES(${SWIG_MODULE_${_mod}_REAL_NAME} ${SWIG_MODULE_Kernel_REAL_NAME})

  SWIG_LINK_LIBRARIES(${_mod} ${PYTHON_LIBRARIES} 
    ${${PROJECT_NAME}_LINK_LIBRARIES})
  
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${_mod}.py 
    DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos/Mechanics/${_dest})

  INSTALL(TARGETS ${SWIG_MODULE_${_mod}_REAL_NAME}
    DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos/Mechanics/${_dest})

ENDMACRO(MECHANICS_ADD_MODULE)

MECHANICS_ADD_MODULE(Bodies .)
MECHANICS_ADD_MODULE(Joints .)

INSTALL(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  DESTINATION ${PYTHON_DIST_PACKAGES}/Siconos/Mechanics)

SUBDIRS(ContactDetection)
