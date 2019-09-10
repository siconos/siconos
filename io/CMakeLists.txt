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

set(COMPONENT io)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  src/.)

if(HAVE_SICONOS_MECHANICS)
  list(APPEND ${COMPONENT}_DIRS src/mechanics) 
endif()

if(WITH_VTK)
  list(APPEND ${COMPONENT}_DIRS VTK)
endif()

# --- Serialization setup ---

if(WITH_SERIALIZATION)
  if(WITH_GENERATION)
    set(GENERATED_HEADER ${CMAKE_CURRENT_BINARY_DIR}/SiconosFullGenerated.hpp)
    set(GENERATED_HEADER_KERNEL ${CMAKE_CURRENT_BINARY_DIR}/SiconosFullKernelGenerated.hpp)
    # FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullGenerated.hpp)
    # FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullKernelGenerated.hpp)

    set(io_targets ${COMPONENTS})

    list(REMOVE_ITEM io_targets numerics io externals)
    print_var(io_targets)

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
    
    if(MPI_FOUND)
      foreach(_DIR IN LISTS MPI_CXX_INCLUDE_PATH)
        list(APPEND GENERATED_INCLUDES -I${_DIR})
      endforeach()
    endif()

    set(GENERATED_COMMAND COMMAND ${CMAKE_SOURCE_DIR}/io/tools/build_from_doxygen.py
      ${GENERATED_TARGETS} ${GENERATED_INCLUDES} --output=${GENERATED_HEADER}
      --source="${CMAKE_CURRENT_SOURCE_DIR}/../" --build="${CMAKE_BINARY_DIR}")

    add_custom_command(
      OUTPUT ${GENERATED_HEADER}
      DEPENDS kernel-doxy2xml ${GENERATED_DEPENDENCIES}
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
set(${COMPONENT}_EXCLUDE_DOXY src/utils/SiconosTools
  # if included, it breaks swig process for mechanics
  )

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
set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES ${${COMPONENT}_DIRS})

# wrap - Set extra dependencies
if(WITH_${COMPONENT}_PYTHON_WRAPPER)
  set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)
  # List of modules to be built in the python bindings
  set(${COMPONENT}_PYTHON_MODULES "${COMPONENT}")
endif()#


# ---- Final setup for the library ----
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

list(APPEND ${COMPONENT}_LINK_LIBRARIES kernel numerics)
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
  target_link_libraries(io PUBLIC mechanisms)
endif()
# # remove debug & optimized inside LINK_LIBRARIES...
# list(REMOVE_ITEM ${COMPONENT}_LINK_LIBRARIES "debug" "optimized")

# Links with non-Siconos libraries

if(SICONOS_HAS_OCE)
  target_link_libraries(io PRIVATE ${OCE_LIBRARIES})
endif()
# --- Search component dependencies ---
if(WITH_SERIALIZATION)
  include(serialization_vector_test)
  find_package(Boost 1.61 COMPONENTS serialization REQUIRED)
  find_package(Boost 1.61 COMPONENTS filesystem REQUIRED)
  target_link_libraries(io PUBLIC Boost::boost)
  set(WITH_SYSTEM_BOOST_SERIALIZATION ON)
  TEST_SERIALIZATION_VECTOR_BUG()
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

if(WITH_PYTHON_WRAPPER)
  add_dependencies(${COMPONENT} ${COMPONENT}_docstrings)
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
set(${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES
  )

siconos_component_install_setup(${COMPONENT})

# --- specific install ----
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullGenerated.hpp
  DESTINATION include/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullNumerics.hpp
    DESTINATION include/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/Register.hpp
    DESTINATION include/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosRestart.hpp
    DESTINATION include/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFull.hpp
    DESTINATION include/${PROJECT_NAME})
if(HAVE_SICONOS_MECHANICS)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/mechanics/MechanicsIO.hpp
    DESTINATION include/${PROJECT_NAME})
endif()

# --- tests ---
include(${COMPONENT}_tests.cmake)





