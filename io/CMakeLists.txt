#================================================================
# cmake utilities to build IO component
#================================================================
#
# The objective is to call component_setup to create the target <COMPONENT>.
# Before, it's necessary to set:
# 
# - COMPONENT component name
# - <COMPONENT>_DIRS: the list of paths (relative to CMAKE_CURRENT_SOURCE_DIR) that
#   contain source files
# - <COMPONENT>_EXCLUDE_SRCS: the list of files, in <COMPONENT>_DIRS, that must be excluded
#   from build.
# - <COMPONENT>_INTERFACE_INCLUDE_DIRECTORIES: a list of directories
#   to populate the interface of the target for include directories at build time

# Component name (i.e. target name)
set(COMPONENT io)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# ------ source directories for current component ------
# What is needed by component to compile ?
# List here all directories that contain sources files
# for current component.
# Path must be relative to component path (i.e. to CMAKE_CURRENT_SOURCE_DIR)
set(${COMPONENT}_DIRS)

# By default this component has no sources, everything is optional :
# - serialization/generation part : sources in serialization and generation
# - mechanicsIO : sources in mechanics 

if(HAVE_SICONOS_MECHANICS)
  list(APPEND ${COMPONENT}_DIRS src/mechanics) 
endif()

# --- Serialization setup ---

if(WITH_SERIALIZATION)
  list(APPEND ${COMPONENT}_DIRS src/serialization)
  
  if(NOT WITH_GENERATION)
    # in that case, we use headers in src/generation.
    # Warning : they must regularly updated (build siconos with generation ON!)
    # to take into account all changes in siconos interface.
    # --- Some headers are generated in current binary dir ---
    list(APPEND ${COMPONENT}_DIRS src/generation)
  else()
    set(GENERATED_HEADER ${CMAKE_CURRENT_BINARY_DIR}/SiconosFullGenerated.hpp)
    set(GENERATED_HEADER_KERNEL ${CMAKE_CURRENT_BINARY_DIR}/SiconosFullKernelGenerated.hpp)

    # list components 'targeted' by generation
    set(io_targets ${COMPONENTS})
    list(REMOVE_ITEM io_targets numerics io externals)
    string(REPLACE ";" "," io_targets_list "${io_targets}")
    set(GENERATED_TARGETS --targets=${io_targets_list})

    if(HAVE_SICONOS_KERNEL)
      list(APPEND GENERATED_DEPENDENCIES ${CMAKE_SOURCE_DIR}/kernel/src/SiconosKernel.hpp)
    endif()
    if(HAVE_SICONOS_MECHANICS)
      list(APPEND GENERATED_DEPENDENCIES ${CMAKE_SOURCE_DIR}/mechanics/src/collision/native/SpaceFilter.hpp)
    endif()
    if(HAVE_SICONOS_CONTROL)
      list(APPEND GENERATED_DEPENDENCIES ${CMAKE_SOURCE_DIR}/control/src/SiconosControl.hpp)
    endif()
    list(APPEND GENERATED_INCLUDES -I${CMAKE_BINARY_DIR}) # For SiconosConfig.h
    if(SuiteSparse_CXSparse_INCLUDE_DIR)
      list(APPEND GENERATED_INCLUDES -I${SuiteSparse_CXSparse_INCLUDE_DIR})
    endif()
    
    if(MPI_FOUND)
      foreach(_DIR IN LISTS MPI_CXX_INCLUDE_PATH)
        list(APPEND GENERATED_INCLUDES -I${_DIR})
      endforeach()
    endif()

    set(GENERATED_COMMAND COMMAND ${CMAKE_SOURCE_DIR}/io/tools/build_from_doxygen.py
      ${GENERATED_TARGETS} ${GENERATED_INCLUDES} --output=${GENERATED_HEADER}
      --source="${CMAKE_CURRENT_SOURCE_DIR}/../" --build="${CMAKE_BINARY_DIR}")

		## CMake >= 3.12
    # set(io_targets_doxy2xml "${io_targets}")
    # list(TRANSFORM io_targets_doxy2xml APPEND "-doxy2xml")
    set(io_targets_doxy2xml)
		foreach(_T IN LISTS io_targets)
			list(APPEND io_targets_doxy2xml "${_T}-doxy2xml")
		endforeach()

    add_custom_command(
      OUTPUT ${GENERATED_HEADER}
      DEPENDS ${io_targets_doxy2xml} ${GENERATED_DEPENDENCIES}
      ${CMAKE_SOURCE_DIR}/io/tools/build_from_doxygen.py
      ${CMAKE_SOURCE_DIR}/io/tools/builder_common.py
      ${GENERATED_COMMAND})

    SET(GENERATED_COMMAND_KERNEL COMMAND ${CMAKE_SOURCE_DIR}/io/tools/build_from_doxygen.py
      --targets=kernel ${GENERATED_INCLUDES}
      --output=${GENERATED_HEADER_KERNEL}
      --source="${CMAKE_CURRENT_SOURCE_DIR}/../" --build="${CMAKE_BINARY_DIR}")

    ADD_CUSTOM_COMMAND(
      OUTPUT ${GENERATED_HEADER_KERNEL}
      DEPENDS kernel-doxy2xml ${GENERATED_DEPENDENCIES}
      ${GENERATED_COMMAND_KERNEL})
  endif()
endif()

# -- Documentation --
# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
# set(${COMPONENT}_EXCLUDE_DOXY)

# ------ include interface ------
# What is needed at build time
# by other targets to compile with current component.
# 
# It means that a call to
#  target_link_libraries(truc PRIVATE io)
# will imply -I<dirs> with dirs listed in
# ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES.
#
# If all dirs are required, just set
# set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES ${${COMPONENT}_DIRS})
# set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES ${${COMPONENT}_DIRS})

# ---- Final setup for the library ----
if(io_DIRS) # the component may be empty if mechanics and serialization are off
  include(ComponentSetup)
  create_siconos_component(${COMPONENT})
  
  # --- Extra setup for the component ---
  # Right now boost::serialization does some nasty things... so we have to
  # disable some undefined behavior checks ...
  if(USE_SANITIZER MATCHES "asan")
    target_compile_options(io PRIVATE "-fno-sanitize=vptr")
    target_compile_options(io PRIVATE "-fno-sanitize=alignment")
    target_compile_options(io PRIVATE "-fno-sanitize=null")
  endif()
  
  if(WITH_GENERATION)
    # add dependencies for header generation
    add_custom_target(SerializersGeneration ALL
      DEPENDS ${GENERATED_HEADER} ${GENERATED_HEADER_KERNEL})
    add_dependencies(_kernel SerializersGeneration)
    add_dependencies(${COMPONENT} SerializersGeneration)
  endif()
  
  # Links with other Siconos components
  
  if (CMAKE_SKIP_RPATH)
    # if no RPATH, then linking does not work for tests without specifying externals
    target_link_libraries(io PRIVATE externals)
  endif()
  
  target_link_libraries(io PRIVATE numerics)
  target_link_libraries(io PUBLIC kernel)
  
  if(HAVE_SICONOS_CONTROL)
    target_link_libraries(io PRIVATE control)
  endif()
  
  if(HAVE_SICONOS_MECHANICS)
    target_link_libraries(io PUBLIC mechanics)
  endif()
  
  if(HAVE_SICONOS_MECHANISMS)
    # FP : do we really need this link?
    target_link_libraries(io PUBLIC mechanisms)
  endif()

  # Links with non-Siconos libraries
  if(SICONOS_HAS_OCE)
    if(HAVE_SICONOS_MECHANICS OR HAVE_SICONOS_MECHANISMS)
      target_link_libraries(io PRIVATE ${OCE_LIBRARIES})
    endif()
  endif()
  
  # --- Search component dependencies ---
  if(WITH_SERIALIZATION)
    # boost must have been searched for in SiconosSetup ...
    target_link_libraries(io PUBLIC Boost::boost Boost::serialization Boost::filesystem)
    # TEST_SERIALIZATION_VECTOR_BUG() : not needed anymore, since we required boost >= 1.61
    # These test macros have been moved to siconos-junk.

    target_include_directories(io PRIVATE serialization)
    if(NOT WITH_GENERATION)
      target_include_directories(io PRIVATE generation)
    else()
      target_include_directories(io PRIVATE ${CMAKE_BINARY_DIR}/io)
    endif()
  endif()
  
  if(WITH_VTK)
    # https://cmake.org/cmake/help/latest/module/FindVTK.html
    find_package(VTK )
    target_link_libraries(io PRIVATE vtkIO;vtkCommon;vtkGraphics;vtkRendering)
  endif()

  # --- python bindings ---
  if(WITH_${COMPONENT}_PYTHON_WRAPPER)
    add_subdirectory(swig)
  endif()
  
  # ---- Installation ----
  # Call siconos_component_install_setup(<COMPONENT>)
  # to prepare installation of the current target.
  #
  # Before, it's necessary to set:
  # 
  # - <COMPONENT>_INSTALL_INTERFACE_INCLUDE_DIRECTORIES with all directories
  #    that contain headers files that must be installed.
  # 
  # set(${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES)
  if(WITH_SERIALIZATION)
    list(APPEND ${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES src/serialization)
    if(NOT WITH_GENERATION)
      list(APPEND ${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES src/generation)
    endif()
  endif()
  
  if(HAVE_SICONOS_MECHANICS)
    list(APPEND ${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES src/mechanics)
  endif()
  siconos_component_install_setup(${COMPONENT})
  
  # --- tests ---
  include(${COMPONENT}_tests.cmake)
endif()




