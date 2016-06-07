#===============================================================================
# cmake utility to compile and install siconos.
#
#
# Siconos development team - September 2015.
#
#===============================================================================

# ============= Global 'standard' cmake Settings =============

# Set minimum version for cmake
cmake_minimum_required(VERSION 3.0.2)

# Set cmake modules directory (i.e. the one which contains all user-defined FindXXX.cmake files among other things)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Force out-of-source build
include(OutOfSourcesBuild)

# Some useful macros
include(SiconosTools)

# Check/set uninitialized cache variables (--> for docker config)
get_uninitialized_vars(CMAKE_ARGS_FROM_CACHE)

# cmake project name
set(PROJECT_NAME siconos)
set(PROJECT_NAME_CAP Siconos)

if(USER_OPTIONS_FILE)
  get_filename_component(input_file ${USER_OPTIONS_FILE} NAME)
  get_filename_component(path_to_file ${USER_OPTIONS_FILE} PATH)
  find_file(filefound NAMES ${input_file}
    PATHS ${path_to_file} ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    )
  if(NOT ${filefound}-NOTFOUND)
    configure_file(${path_to_file}/${input_file} ${CMAKE_SOURCE_DIR}/${input_file})
  else()
    message(FATAL_ERROR "Cannot find file ${USER_OPTIONS_FILE}")
  endif()
  message("\n !!!!! Load user-defined options set from file ${filefound} !!!!! \n")
  include(${input_file})
  # copy user file to binary dir, just for 'visual' check.
  configure_file(${input_file} ${input_file} COPYONLY)
else()
  include(default_options)
endif()

# Configurations used for continuous integration
if(CI_CONFIG)
  string(REPLACE "," ";" CI_CONFIG_LIST ${CI_CONFIG})
  FOREACH(_CI ${CI_CONFIG_LIST})
    include(CI/config/${_CI}.cmake)
  ENDFOREACH(_CI CI_CONFIG_LIST)
endif()

# ---------- Docker or not docker ... ------------------------
if(WITH_DOCKER)
  find_program(HAVE_DOCKER docker)
  if(HAVE_DOCKER)
    set(DOCKER_MKSENV_COMMAND ${CMAKE_SOURCE_DIR}/CI/machinery/mksenv.py)
    set(DOCKER_MKSENV_INPUT ${CMAKE_SOURCE_DIR}/CI/config/siconos.yml)
    set(DOCKER_REPOSITORY siconos)
    set(DOCKER_CMAKE_FLAGS ${CMAKE_ARGS_FROM_CACHE};-DWITH_DOCKER=0;-V)
    set(DOCKER_SHARED_DIRECTORIES ${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/examples/Control/PID;${CMAKE_SOURCE_DIR}/kernel/src/model/test;${DOCKER_SHARED_DIRECTORIES})
    set(DOCKER_CTEST_DRIVER ${CMAKE_SOURCE_DIR}/CI/machinery/CTestDriver.cmake)
    include(Docker)
    
    add_docker_targets(
      DOCKER_IMAGE ${DOCKER_IMAGE}
      DOCKER_TEMPLATE ${DOCKER_TEMPLATE}
      DOCKER_IMAGE_DIR ${DOCKER_IMAGE_DIR}
      DOCKER_REPOSITORY ${DOCKER_REPOSITORY}
      DOCKER_CMAKE_FLAGS ${DOCKER_CMAKE_FLAGS}
      DOCKER_SHARED_DIRECTORIES ${DOCKER_SHARED_DIRECTORIES}
      DOCKER_CTEST_DRIVER ${DOCKER_CTEST_DRIVER})
      
    project(${PROJECT_NAME} NONE)
  else()
    message(FATAL_ERROR "Cannot find docker.")
  endif()
else()
  # Set target names from components list
  assert(COMPONENTS)
  set(AVAILABLE_COMPONENTS "externals numerics kernel control mechanics io")
  foreach(component_name ${COMPONENTS})
    string(TOUPPER ${component_name} upname)
    list(FIND ${AVAILABLE_COMPONENTS} ${component_name} comp_exists)
    if(NOT comp_exists)
      message(FATAL_ERROR "Unknown siconos component : ${component_name}")
    endif()
    set(HAVE_SICONOS_${upname} TRUE)
  endforeach()
  
  # ============= The project =============
  # Set project name and project languages
  # => this automatically defines:
  #   - ${PROJECT_NAME}_BINARY_DIR : where you have run cmake, i.e. the place for compilation
  #   - ${PROJECT_NAME}_SOURCE_DIR : where sources (.f and .h and this CMakeLists.txt) are located
  # Note that because of OutOfSourceBuild, binary_dir and source_dir must be different.
  if(WITH_CXX)
    # Set CXX compilation flags
    project(${PROJECT_NAME} CXX C Fortran)
    include(CXXCompilerSetup)
  else() # we may need to compile only numerics and only in C.
    project(${PROJECT_NAME} C Fortran)
  endif()
  # Set C compilation flags
  include(CCompilerSetup)
  # Fortran compiler configuration
  include(FCompilerSetup)
  
  # Common setup
  include(SiconosSetup)
  
  # ========= Search for 'common' dependencies =========
  include(siconos_find_common_dep)
  
  # ========= Python bindings =========
  if(WITH_PYTHON_WRAPPER)
    add_subdirectory(wrap)
  endif()
  
  # ========= Tests setup =========
  if(WITH_TESTING)
    include(CTest)
    include(SiconosCTest)
    # log file for tests setup
    file(WRITE ${CMAKE_BINARY_DIR}/tests.log "--- Siconos tests setup --- \n")
  endif()
  
  # ========= Documentation setup =========
  # This includes doxygen_warnings,
  # doxygen doc from source files, sphinx doc
  # and docstrings from doxygen, depending on the
  # enabled options.
  add_subdirectory(Docs)
  
  # ========= Setup each required component =========
  foreach(COMPONENT ${COMPONENTS})
    message("\n--------------------------------------------------")
    add_subdirectory(${COMPONENT})
  endforeach()
  
  finalize_doxygen()
  
  include(FeatureSummary)
  feature_summary(WHAT ALL)
  
  # =========== Generate siconos scripts ===========
  if(EXISTS ${CMAKE_SOURCE_DIR}/scripts/CMakeLists.txt.in)
    configure_file(scripts/CMakeLists.txt.in scripts/CMakeLists-temp.txt @ONLY)
    configure_file(${CMAKE_BINARY_DIR}/scripts/CMakeLists-temp.txt ${CMAKE_BINARY_DIR}/scripts/CMakeLists.txt @ONLY)
    install(FILES ${CMAKE_BINARY_DIR}/scripts/CMakeLists.txt DESTINATION share/${PROJECT_NAME})
  endif()
  if(EXISTS ${CMAKE_SOURCE_DIR}/scripts/siconos.py.in)
    message("-- Generate siconos script ...")
    configure_file(scripts/siconos.py.in scripts/siconos @ONLY)
    install(PROGRAMS ${CMAKE_BINARY_DIR}/scripts/siconos DESTINATION bin)
  endif()
  
  # ============= Save config for other cmake projects ============
  include(CMakePackageConfigHelpers)
  
  # Generate ${PROJECT_NAME}Config.cmake
  configure_package_config_file(Config.cmake.in ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/cmake)
  install(FILES ${CMAKE_BINARY_DIR}/SiconosConfig.h DESTINATION include/${PROJECT_NAME})

  # Generate configVersion.cmake file.
  write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}
    COMPATIBILITY ExactVersion
    )
  
  # install both of them
  install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake ${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION share/${PROJECT_NAME}/cmake)
  
  # ============= Summary =============
  message(STATUS "\n============================================ Summary ============================================")
  message(STATUS "${PROJECT_NAME} version ${SICONOS_VERSION} is now ready for compilation and installation.")
  message(STATUS "To proceed run 'make' and 'make install' and optionaly 'make test'.")
  message(STATUS "C++ Compiler : ${CMAKE_CXX_COMPILER}")
  message(STATUS "C Compiler : ${CMAKE_C_COMPILER}")
  message(STATUS "Fortran Compiler : ${CMAKE_Fortran_COMPILER}")
  message(STATUS "Compilation mode is : ${CMAKE_BUILD_TYPE}")
  message(STATUS "Code Sources are in : ${CMAKE_SOURCE_DIR}")
  message(STATUS "Blas is ${WITH_BLAS}. ")
  message("    Blas libraries : ${BLAS_LIBRARIES}.")
  message("    Blas headers : ${BLAS_HEADER} in ${BLAS_INCLUDE_DIRS}.")
  message(STATUS "Lapack is ${WITH_LAPACK}.")
  message("    Lapack libraries : ${LAPACK_LIBRARIES}.")
  message("    Lapack headers : ${LAPACK_HEADER} in ${LAPACK_INCLUDE_DIRS}.")
  message(STATUS "Python interpreter is ${PYTHON_EXECUTABLE}")
  message(STATUS "Python libraries are ${PYTHON_LIBRARIES}")
  message(STATUS "All linking libraries are : ${SICONOS_LINK_LIBRARIES}.")
  message(STATUS "The project will be installed in ${CMAKE_INSTALL_PREFIX}")
  if(WITH_PYTHON_WRAPPER)
    message(STATUS "${SICONOS_PYTHON_PACKAGE} python package will be installed in ${SICONOS_PYTHON_INSTALL_DIR}")
  endif()
  message(STATUS "To get more information about dependencies, config or else, ")
  message(STATUS "check CMakeCache.txt file or re-run cmake with -DPRINT_ENV=ON.")
  message(STATUS "=================================================================================================\n")
  
  # Log build info into ${CMAKE_BINARY_DIR}/Testing/Notes/Build
  write_notes()
endif()
