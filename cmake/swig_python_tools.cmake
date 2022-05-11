# =======================================
# Macros and functions related to swig
# and python
# =======================================

#[=======================================================================[.rst:
Create a python module inside siconos package, from a swig (.i) file

Usage :

.. code-block:: cmake

add_swig_sub_module(
  FILE filename                # name of the swig file used to create the module,
                               # including path relative to CMAKE_CURRENT_SOURCE_DIR
		               # The package will be named after this filename.
  [DEPS]   <deps>              # list of targets on which this package depends
  [INCLUDES] <includes>        # list of paths that must be known by swig to generate the package
  [COMPILE_OPTIONS] <opts>     # compile options required by swig
  [COMPILE_DEFINITIONS] <opts> # compile definitions (-D) required by swig
    )


  
  Requirements to build module 'truc' in package siconos.osnspb:
  - osnspb/truc.i starting with %module(package="siconos.osnspb", ...
  - osnspb/__init__.py.in or osnspb/__init__.py file in ${CMAKE_CURRENT_SOURCE_DIR}/

  For example:

  add_swig_sub_module(
    FILE osnspb/truc.i
    )

  will build/create a module siconos.osnspb.truc with swig, from file ${CMAKE_CURRENT_SOURCE_DIR}/osnspb/truc.i
  and install it in ${SICONOS_PYTHON_INSTALL_DIR}/${SICONOS_PYTHON_PACKAGE}/osnspb, to allow in python:

  import siconos.osnspb.truc


  
#]=======================================================================]
function(add_swig_sub_module)

  set(oneValueArgs FILE)
  set(multiValueArgs DEPS INCLUDES COMPILE_OPTIONS COMPILE_DEFINITIONS)

  cmake_parse_arguments(target "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set (UseSWIG_TARGET_NAME_PREFERENCE STANDARD)

  # We use swig file name as the target name.
  get_filename_component(python_module_name ${target_FILE} NAME_WE)
  get_filename_component(python_module_path ${target_FILE} DIRECTORY)
    
  # --- Swig main file ---
  # Set with arg FILE of the current function. Assumed to be relative to
  # CMAKE_CURRENT_SOURCE_DIR, which must be ${COMPONENT}/swig
  set(swig_file ${CMAKE_CURRENT_SOURCE_DIR}/${target_FILE})

  if(python_module_path)
    message(" -- Build module siconos.${python_module_path}.${python_module_name} using swig file ${swig_file}.")
  else()
    message(" -- Build module siconos.${python_module_name} using swig file ${swig_file}.")
  endif()
  # If a target with that name  already exists, we prepend py to this name.
  set(pythonfile_name ${python_module_name}.py)
  #if(TARGET ${python_module_name})
  set(target_NAME py${python_module_name})
  #endif()


  # --- Set properties for the current swig file ---

  # Properties on sources

  set_property(SOURCE ${swig_file} PROPERTY CPLUSPLUS ON) # c++

  #set_property(SOURCE ${swig_file} PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON) # propagates includes from target to swig.

  
  # if(WITH_SERIALIZATION)
  #   set_source_files_properties(${swig_file}
  #     PROPERTIES SWIG_FLAGS "${${COMPONENT}_SWIG_DEFS_${_name}}" CPLUSPLUS ON)
  # endif()

  # --- build swig module ---
  #set(ADDITIONAL_SWIG_DEFINES ${ADDITIONAL_SWIG_DEFINES} -DBOOST_NOEXCEPT)

  # ---- Build the library ---  
  swig_add_library(${target_NAME}
    TYPE MODULE
    LANGUAGE python
    OUTPUT_DIR "${SICONOS_SWIG_ROOT_DIR}/${python_module_path}" # where to write the language specific files
    OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR}   # where the generated source file will be placed 
    SOURCES ${swig_file})
  
  set_property(TARGET ${target_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/${python_module_path})
  #  set_property(TARGET ${target_NAME} PROPERTY SWIG_COMPILE_DEFINITIONS )

  # Forward includes of the target for the swig command
  set_property(TARGET ${target_NAME} PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON)
  if(target_INCLUDES)
    set_property(TARGET ${target_NAME} PROPERTY SWIG_INCLUDE_DIRECTORIES ${target_INCLUDES})
  endif()

  if(target_COMPILE_OPTIONS)
    set_property(TARGET ${target_NAME} PROPERTY SWIG_COMPILE_OPTIONS ${target_COMPILE_OPTIONS})
  endif()
  if(target_COMPILE_DEFINITIONS)
    set_property(TARGET ${target_NAME} PROPERTY SWIG_COMPILE_DEFINITIONS ${target_COMPILE_DEFINITIONS})
  endif()
  
  # Link with current component
  foreach(dep IN LISTS target_DEPS)
    target_link_libraries(${target_NAME} PRIVATE ${dep})
  endforeach()
  # Python and numpy
  target_link_libraries(${target_NAME} PRIVATE Python3::NumPy)
  target_include_directories(${target_NAME} PRIVATE ${SICONOS_SWIG_SOURCE_DIR})
  target_include_directories(${target_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${target_NAME} PRIVATE ${target_INCLUDES})

  # List of siconos modules, used in __init__.py.
  
  # if(WITH_SERIALIZATION)
  #   target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/io/src/serialization")
  #   target_link_libraries(${SWIG_MODULE_${_name}_REAL_NAME} Boost::serialization)
  #   if(NOT WITH_GENERATION)
  #     target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/io/src/generation")
  #   else()
  #     add_dependencies(${SWIG_MODULE_${_name}_REAL_NAME} SerializersGeneration)
  #     target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE "${CMAKE_BINARY_DIR}/io")
  #   endif()

  #   # SiconosFullGenerated.hpp includes files from all other components.
  #   # (Better way than using *_DOXYGEN_INPUTS?  ${COMPONENT}_DIR is empty here!)
  #   foreach(_C IN LISTS COMPONENTS)
  #     string(STRIP "${${_C}_DOXYGEN_INPUTS}" _dirs)
  #     string(REPLACE " " ";" _dirs "${_dirs}")
  #     foreach(_D IN LISTS _dirs)
  #       target_include_directories(${SWIG_MODULE_${_name}_REAL_NAME} PRIVATE ${_D})
  #     endforeach()
  #   endforeach()
  # endif()

  # WARNING ${swig_generated_file_fullname} is overriden 
  #set(${_name}_generated_file_fullname ${swig_generated_file_fullname})
  #set_source_files_properties( ${swig_generated_file_fullname}
  #    PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing")
  # Set path for the library generated by swig for the current module --> siconos python package path
  #set_property(TARGET ${SWIG_MODULE_${_name}_REAL_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/${_path})
  #message(" -- ${_name} generated (swig) file will be ${${_name}_generated_file_fullname}")


  # Failed on Macos. Something is different for version on Apple. See https://cmake.org/cmake/help/v3.16/prop_tgt/SOVERSION.html?highlight=mach.
  if(NOT APPLE)
    set_property(TARGET ${target_NAME} PROPERTY NO_SONAME OFF)
    set_property(TARGET ${target_NAME} PROPERTY VERSION "${SICONOS_SOVERSION}")
    set_property(TARGET ${target_NAME} PROPERTY SOVERSION "${SICONOS_SOVERSION_MAJOR}")
  endif()

  # IF(MSVC AND ${COMPONENT} MATCHES "kernel")
  #   set_source_files_properties(${${_name}_generated_file_fullname} PROPERTIES COMPILE_FLAGS "/bigobj")
  # ENDIF()

  # Add a post-build step that prepends utf-8 coding indicator to .py files
  # find_program(SH_COMMAND sh)
  # find_program(CAT_COMMAND cat)
  # find_program(MV_COMMAND mv)
  # if(SH_COMMAND AND CAT_COMMAND)
  #   if(MV_COMMAND)
  #     add_custom_command(TARGET ${SWIG_MODULE_${_name}_REAL_NAME}
  #       POST_BUILD COMMAND ${SH_COMMAND} -c "(echo '# -*- coding: utf-8 -*-'; ${CAT_COMMAND} ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.py) > ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.tmp; ${MV_COMMAND} ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.tmp ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.py" VERBATIM)
  #   endif()
  # endif()
  # # Check dependencies and then link ...
  # add_dependencies(${SWIG_MODULE_${_name}_REAL_NAME} ${COMPONENT})
  # if(UNIX AND NOT APPLE)
  #   # do not link against the Python library on unix, it is useless
  #   swig_link_libraries(${_name} ${${COMPONENT}_LINK_LIBRARIES} ${COMPONENT})
  # else()
  #   swig_link_libraries(${_name} ${Python3_LIBRARIES} ${${COMPONENT}_LINK_LIBRARIES} ${COMPONENT})
  # endif()

  if(WITH_DOCUMENTATION AND WITH_DOXY2SWIG)
    include(doc_tools)
    # Create the target to generate rst/sphinx from docstrings (swig outputs).
    docstrings2rst(PATH ${python_module_path} NAME ${python_module_name})
    add_dependencies(${python_module_name}_autodoc ${target_NAME})
  endif()
  
  # Copy __init__.py file if needed
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${python_module_path}/__init__.py.in)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${python_module_path}/__init__.py.in
      ${SICONOS_SWIG_ROOT_DIR}/${python_module_path}/__init__.py @ONLY)
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${python_module_path}/__init__.py)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${python_module_path}/__init__.py
      ${SICONOS_SWIG_ROOT_DIR}/${python_module_path}/__init__.py COPYONLY)
  endif()

  # --- install python files and target ---
  # install path ...
  set(DEST "${SICONOS_PYTHON_INSTALL_DIR}/${SICONOS_PYTHON_PACKAGE}/${python_module_path}")

  # install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${_name}.py DESTINATION ${DEST})
  install(TARGETS ${target_NAME} LIBRARY DESTINATION ${DEST})
  
endfunction()




include(tools4tests)
