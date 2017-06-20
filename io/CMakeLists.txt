#================================================================
# cmake utilities to build IO component
#================================================================

set(COMPONENT io)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})

# --- Search component dependencies ---
if(WITH_SERIALIZATION)
  find_package(Boost 1.47 COMPONENTS serialization REQUIRED)
  find_package(Boost 1.47 COMPONENTS filesystem REQUIRED)
  # cmake 2.8.12.2 (ubuntu 14.04) need this
  set(${COMPONENT}_LINK_LIBRARIES  "${${COMPONENT}_LINK_LIBRARIES};${Boost_LIBRARIES}")
endif()
# Right now boost::serialization does some nasty things... so we have to
# disable some undefined behavior checks ...
IF(USE_SANITIZER MATCHES "asan")
  ADD_CXX_OPTIONS("-fno-sanitize=vptr")
  ADD_CXX_OPTIONS("-fno-sanitize=alignment")
  ADD_CXX_OPTIONS("-fno-sanitize=null")
ENDIF(USE_SANITIZER MATCHES "asan")
# Note FP : we should move this to CXXCompilerSetup?
# Answer xhub : no this is local, because boost::serialization is hackish

# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  src/.)

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS
  )
if(HAVE_SICONOS_MECHANICS)
  list(APPEND ${COMPONENT}_DIRS src/mechanics)
endif()

if(WITH_VTK)
  list(APPEND ${COMPONENT}_DIRS VTK)
  set(${COMPONENT}_LINK_LIBRARIES "${${COMPONENT}_LINK_LIBRARIES};vtkIO;vtkCommon;vtkGraphics;vtkRendering")
endif()

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR)

# --- Serialization setup ---
set(GENERATED_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullGenerated.hpp)
set(GENERATED_HEADER_KERNEL ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullKernelGenerated.hpp)


if(WITH_SERIALIZATION)
  if(WITH_GENERATION)
    FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullGenerated.hpp)
    FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/src/SiconosFullKernelGenerated.hpp)

    set(io_targets ${COMPONENTS})

    # REMOVE_item not on cmake 2.8 ...
    list(REMOVE_ITEM io_targets numerics io externals)
    print_var(io_targets)

    string(REPLACE ";" "," io_targets_list "${io_targets}")
    SET(GENERATED_TARGETS --targets=${io_targets_list})

    if(HAVE_SICONOS_KERNEL)
      set(GENERATED_DEPENDENCIES ../kernel/src/SiconosKernel.hpp)
    endif()
    if(HAVE_SICONOS_MECHANICS)
      set(GENERATED_DEPENDENCIES ${GENERATED_DEPENDENCIES}
        ../mechanics/src/collision/native/SpaceFilter.hpp)
    endif()
    if(HAVE_SICONOS_CONTROL)
      set(GENERATED_DEPENDENCIES ${GENERATED_DEPENDENCIES} ../control/src/SiconosControl.hpp)
    endif()

    # append local headers directories to GENERATED_INCLUDES)
    foreach(_DIR ${${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES})
      list(APPEND GENERATED_INCLUDES -I${_DIR})
    endforeach()

    if(MPI_FOUND)
      FOREACH(_DIR ${MPI_CXX_INCLUDE_PATH})
        SET(GENERATED_INCLUDES ${GENERATED_INCLUDES} -I${_DIR})
      ENDFOREACH()
    ENDIF()

    SET(GENERATED_COMMAND COMMAND ${CMAKE_SOURCE_DIR}/io/tools/build_from_doxygen.py
      ${GENERATED_TARGETS} ${GENERATED_INCLUDES} --output=${GENERATED_HEADER}
      --source="${CMAKE_CURRENT_SOURCE_DIR}/../" --build="${CMAKE_BINARY_DIR}")

    ADD_CUSTOM_COMMAND(
      OUTPUT ${GENERATED_HEADER}
      DEPENDS doxygen ${GENERATED_DEPENDENCIES}
      ${GENERATED_COMMAND})

    SET(GENERATED_COMMAND_KERNEL COMMAND ${CMAKE_SOURCE_DIR}/io/tools/build_from_doxygen.py
      --targets=kernel ${GENERATED_INCLUDES}
      --output=${GENERATED_HEADER_KERNEL}
      --source="${CMAKE_CURRENT_SOURCE_DIR}/../" --build="${CMAKE_BINARY_DIR}")

    ADD_CUSTOM_COMMAND(
      OUTPUT ${GENERATED_HEADER_KERNEL}
      DEPENDS doxygen ${GENERATED_DEPENDENCIES}
      ${GENERATED_COMMAND_KERNEL})
  endif()
endif()

# wrap - Set extra dependencies
if(WITH_${COMPONENT}_PYTHON_WRAPPER)
  set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)
  # List of modules to be built in the python bindings
  set(${COMPONENT}_PYTHON_MODULES "${COMPONENT}")
endif()#


# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE CXX)
set(${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
list(APPEND ${COMPONENT}_LINK_LIBRARIES kernel)
list(APPEND ${COMPONENT}_LINK_LIBRARIES numerics)

if(HAVE_SICONOS_CONTROL)
  list(APPEND ${COMPONENT}_LINK_LIBRARIES control)
endif()
if(HAVE_SICONOS_MECHANICS)
  list(APPEND ${COMPONENT}_LINK_LIBRARIES mechanics)
endif()

# remove debug & optimized inside LINK_LIBRARIES...
foreach(_L ${${COMPONENT}_LINK_LIBRARIES})
  if(NOT ${_L} STREQUAL "debug")
    if(NOT ${_L} STREQUAL "optimized")
      list(APPEND _LL ${_L})
    endif()
  endif()
endforeach()
set(${COMPONENT}_LINK_LIBRARIES ${_LL})

# add dependencies for header generation
add_custom_target(SerializersGeneration ALL
  DEPENDS ${GENERATED_HEADER} ${GENERATED_HEADER_KERNEL})

include(LibraryProjectSetup)

library_project_setup()

add_dependencies(${COMPONENT} SerializersGeneration)

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
