# =======================================
# Macros and functions related to swig
# and python
# =======================================

# ----------------------------------------------------------------------
# Build a swig module from .i files
#
# Usage :
# add_siconos_swig_submodule(full_name)
#
# full_name must be a path to file.i, path relative to siconos
# python package root.
# For example, to build a swig module 'bullet' in siconos.mechanics.contact_detection,
# call this macro with full_name = mechanics/contact_detection/bullet
# and to build a swig module for kernel in siconos,
# call this macro with full_name = ./kernel
#
# ----------------------------------------------------------------------
macro(add_siconos_swig_sub_module fullname)

  get_filename_component(_name ${fullname} NAME)
  get_filename_component(_path ${fullname} PATH)
  message(" -- Build module ${_name} in directory ${_path} for parent ${COMPONENT}")
  # Add component dependencies to the current submodule deps.
  if(DEFINED SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)
    set(SWIG_MODULE_${_name}_EXTRA_DEPS ${SWIG_MODULE_${COMPONENT}_EXTRA_DEPS})
  endif()

  if(WITH_DOXY2SWIG)
    list(APPEND SWIG_MODULE_${_name}_EXTRA_DEPS ${COMPONENTS}_docstrings)
  endif()
  
  # add as dependencies all the i files
  file(GLOB ${_name}_I_FILES ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/*.i)
  foreach(_f IN LISTS ${_name}_I_FILES)
    list(APPEND SWIG_MODULE_${_name}_EXTRA_DEPS ${_f})
  endforeach()

  # set properties for current '.i' file
  set(swig_file ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/${_name}.i)

  # set output dir
  set(CMAKE_SWIG_OUTDIR "${SICONOS_SWIG_ROOT_DIR}/${_path}")

  # compile flags
  foreach(_dir IN LISTS ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)
    list(APPEND ${COMPONENT}_SWIG_DEFS "-I${_dir}")
  endforeach()
  list(REMOVE_DUPLICATES ${COMPONENT}_SWIG_DEFS)
  
  # extra per-module flags if any
  list(APPEND ${COMPONENT}_SWIG_DEFS_${_name} "${${COMPONENT}_SWIG_DEFS}")

  if(WITH_CXX AND (BUILD_AS_CPP OR NOT ${COMPONENT} MATCHES "numerics"))
    set_source_files_properties(${swig_file}
      PROPERTIES SWIG_FLAGS "${${COMPONENT}_SWIG_DEFS_${_name}}" CPLUSPLUS ON)
  else()
    # C compilation, pass SWIG_FLAGS.
    if(${COMPONENT} MATCHES "numerics")
      set_source_files_properties(${swig_file}
        PROPERTIES SWIG_FLAGS "${${COMPONENT}_SWIG_DEFS_${_name}}")
    endif()
  endif()


  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13")
    set_property(SOURCE ${swig_file} PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)
  else() #USE_TARGET_INCLUDE property is not available in old cmake versions.
    foreach(_dir IN LISTS ${COMPONENT}_DIRS)
      list(APPEND CMAKE_SWIG_FLAGS "-I${CMAKE_SOURCE_DIR}/${COMPONENT}/${_dir}")
    endforeach()
    list(REMOVE_DUPLICATES CMAKE_SWIG_FLAGS)
  endif()

  if(WITH_SERIALIZATION)
    set_source_files_properties(${swig_file}
      PROPERTIES SWIG_FLAGS "${${COMPONENT}_SWIG_DEFS_${_name}}" CPLUSPLUS ON)
  endif()

  # --- build swig module ---
  if(CMAKE_VERSION VERSION_LESS 3.8.0)
    swig_add_module(${_name} python ${swig_file})
  else()
    set(ADDITIONAL_SWIG_DEFINES ${ADDITIONAL_SWIG_DEFINES} -DBOOST_NOEXCEPT)
    swig_add_library(${_name} LANGUAGE python SOURCES ${swig_file})
  endif()

  # Link with current component
  target_link_libraries(${SWIG_MODULE_${_name}_REAL_NAME} ${COMPONENT})
  # Python and numpy
  target_link_libraries(${SWIG_MODULE_${_name}_REAL_NAME} Python3::NumPy)

  # List of siconos modules, used in __init__.py.
  # --> fix this to have 'real' modules names (e.g. control.observer ...)
  set(current_module ${SWIG_MODULE_${_name}_REAL_NAME})
  set(SICONOS_PYTHON_MODULES ""
    CACHE INTERNAL "Modules available in Siconos Python package.")
  # set(SICONOS_PYTHON_MODULES "${SICONOS_PYTHON_MODULES}, '${current_module}'"
  #   CACHE INTERNAL "Modules available in Siconos Python package.")
  
  if(${CMAKE_VERSION} VERSION_LESS "3.13")
    foreach(_dir IN LISTS ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)
      target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE ${_dir})
    endforeach()
  endif()

  if(WITH_SERIALIZATION)
    target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/io/src/serialization")
    if(NOT WITH_GENERATION)
      target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/io/src/generation")
    else()
      add_dependencies(${SWIG_MODULE_${_name}_REAL_NAME} SerializersGeneration)
      target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE "${CMAKE_BINARY_DIR}/io")
    endif()

    # SiconosFullGenerated.hpp includes files from all other components.
    # (Better way than using *_DOXYGEN_INPUTS?  ${COMPONENT}_DIR is empty here!)
    foreach(_C IN LISTS COMPONENTS)
      string(STRIP "${${_C}_DOXYGEN_INPUTS}" _dirs)
      string(REPLACE " " ";" _dirs "${_dirs}")
      foreach(_D IN LISTS _dirs)
        target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE ${_D})
      endforeach()
    endforeach()
  endif()

  # WARNING ${swig_generated_file_fullname} is overriden 
  set(${_name}_generated_file_fullname ${swig_generated_file_fullname})
  set_source_files_properties( ${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing")
  # Set path for the library generated by swig for the current module --> siconos python package path
  set_property(TARGET ${SWIG_MODULE_${_name}_REAL_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/${_path})
  message(" -- ${_name} generated (swig) file will be ${${_name}_generated_file_fullname}")

  # Set the SONAME for the SWIG module to the Siconos SONAME
  set_target_properties(${SWIG_MODULE_${name}_REAL_NAME} PROPERTIES
    NO_SONAME OFF
    VERSION "${SICONOS_SOVERSION}"
    SOVERSION "${SICONOS_SOVERSION_MAJOR}")

  IF(MSVC AND ${COMPONENT} MATCHES "kernel")
    set_source_files_properties(${${_name}_generated_file_fullname} PROPERTIES COMPILE_FLAGS "/bigobj")
  ENDIF()

  # Add a post-build step that prepends utf-8 coding indicator to .py files
  find_program(SH_COMMAND sh)
  find_program(CAT_COMMAND cat)
  find_program(MV_COMMAND mv)
  if(SH_COMMAND AND CAT_COMMAND)
    if(MV_COMMAND)
      add_custom_command(TARGET ${SWIG_MODULE_${_name}_REAL_NAME}
        POST_BUILD COMMAND ${SH_COMMAND} -c "(echo '# -*- coding: utf-8 -*-'; ${CAT_COMMAND} ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.py) > ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.tmp; ${MV_COMMAND} ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.tmp ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.py" VERBATIM)
    endif()
  endif()
  # Check dependencies and then link ...
  add_dependencies(${SWIG_MODULE_${_name}_REAL_NAME} ${COMPONENT})
  if(UNIX AND NOT APPLE)
    # do not link against the Python library on unix, it is useless
    swig_link_libraries(${_name} ${${COMPONENT}_LINK_LIBRARIES} ${COMPONENT})
  else()
    swig_link_libraries(${_name} ${Python3_LIBRARIES} ${${COMPONENT}_LINK_LIBRARIES} ${COMPONENT})
  endif()

  # set dependency of sphinx apidoc to this target
  if(WITH_DOCUMENTATION AND WITH_${COMPONENT}_DOXY2SWIG)
    include(doc_tools)
    docstrings2rst(${_path} ${_name})
  endif()
  
  # Copy __init__.py file if needed
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py.in)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py.in
      ${SICONOS_SWIG_ROOT_DIR}/${_path}/__init__.py @ONLY)
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py
      ${SICONOS_SWIG_ROOT_DIR}/${_path}/__init__.py COPYONLY)
  endif()

  # --- install python files and target ---
  # install path ...
  set(DEST "${SICONOS_PYTHON_INSTALL_DIR}/${SICONOS_PYTHON_PACKAGE}/${_path}")

  #install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${_name}.py DESTINATION ${DEST})
  install(TARGETS ${SWIG_MODULE_${_name}_REAL_NAME} LIBRARY DESTINATION ${DEST})
  
endmacro()

macro(swig_module_setup modules_list)
  if(WITH_${COMPONENT}_PYTHON_WRAPPER)
    # we should use lowercase name for python module (pep8 ...)
    message(" -- Prepare python bindings for component ${COMPONENT} ...")
    # build python bindings
    include_directories(${SICONOS_SWIG_INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)
    if(HAVE_SICONOS_IO)
      include_directories(${CMAKE_SOURCE_DIR}/io/src)
      if(HAVE_SICONOS_MECHANICS)
	include_directories(${CMAKE_SOURCE_DIR}/io/src/mechanics)
      endif()
    endif()
    
    foreach(module IN LISTS ${COMPONENT}_PYTHON_MODULES)
      add_siconos_swig_sub_module(${module})
    endforeach()
  endif()
endmacro()

include(tools4tests)
