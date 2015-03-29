SET(SiconosNumerics_DEPS FALSE)
SET(SiconosKernel_DEPS "Numerics")
SET(SiconosMechanics_DEPS "Kernel")
SET(SiconosControl_DEPS "Kernel")
SET(SiconosIO_DEPS "IO")
SET(SiconosFCLib_DEPS "Numerics")

MACRO(SICONOS_SWIG_MODULE _module _parts)
  SICONOS_SWIG_SUBMODULE(${_module} "${_parts}" .)
ENDMACRO(SICONOS_SWIG_MODULE)

MACRO(SICONOS_SWIG_SUBMODULE _module _parts _subdir)

  SET(DEST ${PYTHON_DIST_PACKAGES}/Siconos/${_module}/${_subdir})

  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION
    ${CMAKE_BINARY_DIR}/src/swig/Siconos/${_module}/${_subdir})

  INSTALL(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
    DESTINATION ${DEST})

  FOREACH(_p ${_parts})
    # add as dependencies all the i files
    IF(DEFINED SWIG_MODULE_${mod}_EXTRA_DEPS)
      SET(SWIG_MODULE_${_p}_EXTRA_DEPS ${SWIG_MODULE_${mod}_EXTRA_DEPS})
    ENDIF()
    FILE(GLOB ${mod}_I_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.i)
    FOREACH(_f ${${mod}_I_FILES})
      LIST(APPEND SWIG_MODULE_${_p}_EXTRA_DEPS ${_f})
    ENDFOREACH()

    SET_SOURCE_FILES_PROPERTIES(${_p}.i PROPERTIES CPLUSPLUS ON)

    SET_SOURCE_FILES_PROPERTIES(${_p}.i 
      PROPERTIES SWIG_FLAGS "-dirprot;-I${CMAKE_SOURCE_DIR}/src/swig/Siconos;-I${CMAKE_BINARY_DIR}/src/swig/Siconos;-I${CMAKE_BINARY_DIR}/src/swig/Siconos/ContacDetection;-Wall;${SWIG_DEFS}")

    SWIG_ADD_MODULE(${_p} python ${_p}.i)

    # WARNING ${swig_generated_file_fullname} is overriden 
    SET(${_p}_generated_file_fullname ${swig_generated_file_fullname})

    MESSAGE(STATUS "${_p} generated file ${${_p}_generated_file_fullname}")

    IF(Siconos${_module}_DEPS)
      ADD_DEPENDENCIES(${SWIG_MODULE_${_p}_REAL_NAME} ${SWIG_MODULE_${Siconos${_module}_DEPS}_DEPNAME})
    ENDIF()

    SWIG_LINK_LIBRARIES(${_p} ${PYTHON_LIBRARIES} ${${PROJECT_NAME}_LINK_LIBRARIES})

    INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${_p}.py 
      DESTINATION ${DEST})

    INSTALL(TARGETS ${SWIG_MODULE_${_p}_REAL_NAME}
      DESTINATION ${DEST})
  ENDFOREACH()

  # we need to save this to the parent scope since we use this in other subdir and also at the toplevel
  SET(SWIG_MODULE_${_module}_DEPNAME ${SWIG_MODULE_${_module}_REAL_NAME} PARENT_SCOPE)
ENDMACRO(SICONOS_SWIG_SUBMODULE)
