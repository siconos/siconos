function(createSwigPythonModule _name _package)

  # Set .i file properties.
  set_source_files_properties(${_name}.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(${_name}.i PROPERTIES SWIG_FLAGS "-dirprot;-I${CMAKE_SOURCE_DIR}/src/swig/siconos;-I${CMAKE_BINARY_DIR}/src/swig/siconos/numerics;-I${CMAKE_SOURCE_DIR}/src/swig/siconos/numerics;-I${CMAKE_CURRENT_BINARY_DIR};-Wall"
)
  
  # Check extra deps
  set(_extraDeps ${ARGN})
  if(_extraDeps)
    foreach(_file ${_extraDeps})
      set(SWIG_MODULE_${_name}_EXTRA_DEPS ${SWIG_MODULE_${_name}_EXTRA_DEPS} ${_file}.i)
    endforeach()
  endif()
  
  print_var(SWIG_MODULE_${_name}_EXTRA_DEPS)

  # Doxygen 
  if(WITH_DOXY)
    set(SWIG_MODULE_${_name}_EXTRA_DEPS ${SWIG_MODULE_${_name}_EXTRA_DEPS} ${_package}_docstrings.i)
  endif()
  
  # Create the swig module
  swig_add_module(${_name} python ${_name}.i)
  # WARNING ${swig_generated_file_fullname} is overriden 
  # save it in a specific var if needed
  set(${_name}_generated_file_fullname ${swig_generated_file_fullname})
  message(STATUS "numerics.${_name} generated file: ${swig_generated_file_fullname}")
  print_var(SWIG_MODULE_${_name}_REAL_NAME )
  add_dependencies(${SWIG_MODULE_${_name}_REAL_NAME} common)
  if(WITH_DOXY)
    add_dependencies(${SWIG_MODULE_${_name}_REAL_NAME} ${_package}_docstrings)
  endif()

  # Links
  swig_link_libraries(
    ${_name} 
    ${PYTHON_LIBRARIES} 
    ${${PROJECT_NAME}_LINK_LIBRARIES}
    ${FORTRAN_LIBRARIES})
  

  # install
  install(PROGRAMS 
    ${CMAKE_CURRENT_BINARY_DIR}/${_name}.py 
    DESTINATION ${PYTHON_DIST_PACKAGES}/siconos/${_package}/
    )
  
  install(TARGETS ${SWIG_MODULE_${_name}_REAL_NAME} DESTINATION ${PYTHON_DIST_PACKAGES}/siconos/${_package})

endfunction()