# Set minimum version for cmake
cmake_minimum_required(VERSION 3.0.2)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Force out-of-source build
include(OutOfSourcesBuild)

# Some usefull macros
include(SiconosTools)

# Check/set uninitialized cache variables (--> for docker config)
get_uninitialized_vars(CMAKE_ARGS_FROM_CACHE)

include(WindowsSiconosSetup)

if(WITH_DOCKER)
  find_program(HAVE_DOCKER docker)
  if(HAVE_DOCKER)
    set(DOCKER_MKSENV_COMMAND ${CMAKE_SOURCE_DIR}/CI/machinery/mksenv.py)
    set(DOCKER_MKSENV_INPUT ${CMAKE_SOURCE_DIR}/CI/config/siconos.yml)
    set(DOCKER_REPOSITORY siconos)
    set(DOCKER_CMAKE_FLAGS ${CMAKE_ARGS_FROM_CACHE};-DWITH_DOCKER=0;-V)
    set(DOCKER_SHARED_DIRECTORIES ${CMAKE_SOURCE_DIR}/..;${DOCKER_SHARED_DIRECTORIES})
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
  project(siconos_examples CXX)
  include(CXXCompilerSetup)
  include(CCompilerSetup)
  
  set(EXAMPLES_DIRECTORIES
    Numerics
    Biology
    Control
    Electronics
    Mechanics
    Robotics/HuMAns_pa10
    Spheres/Lagrangian
    Spheres/NewtonEuler
    )
  
  set(NO_TEST_FILES QGL tools)

  # FIND_PACKAGE(Siconos)
  # STRING(FIND "${installed_targets}" "mechanics" mechanics_installed )
  # IF (${mechanics_installed} LESS 0)
  #   list(APPEND NO_TEST_FILES
  #     JointsTests)
  #   list(APPEND NO_TEST_FILES
  #     NewtonEuler)
  #  endif()
  
  
  # sharing a double between an executable and a library is hackisch
  IF(MSVC)
    SET(NO_TEST_FILES ${NO_TEST_FILES} RockingBlock Pendulum)
  ENDIF(MSVC)
  
  # --- Search siconos package path ---
  # --> based on path to script siconos
  find_file(SICONOS_SCRIPT siconos
    PATHS ENV PATH
    NO_CMAKE_PATH
    NO_DEFAULT_PATH)
  
  if(SICONOS_SCRIPT_NOT_FOUND)
    message(FATAL_ERROR "Can not find siconos script. Please set your PATH environment variable properly.")
  endif()
  
  message("============================================================================================")
  message("Found siconos script : ${SICONOS_SCRIPT}")
  get_filename_component(SICONOS_INSTALL_DIR ${SICONOS_SCRIPT} PATH)
  get_filename_component(SICONOS_INSTALL_DIR ${SICONOS_INSTALL_DIR} PATH)
  message("The examples will be run with siconos installed in : ${SICONOS_INSTALL_DIR}")
  message("============================================================================================")
  
  list(APPEND CMAKE_MODULE_PATH ${SICONOS_INSTALL_DIR}/share/siconos/cmake)
  
  # include siconos config
  include(SiconosVersion)
  find_package(siconos ${SICONOS_VERSION} REQUIRED)
  
  include(CTest)

  # python specific part
  if(HAS_PYTHON_WRAPPER)
    set(PYTHON_EXECUTABLE ${SICONOS_PYTHON_EXECUTABLE})
    include(FindPackageHandleStandardArgs)
    include(FindPythonModule)
    find_python_module(matplotlib)
  endif(HAS_PYTHON_WRAPPER)

  # Exclude some tests if LpSolve has not been found
  if(NOT HAS_ONE_LP_SOLVER)
    set(NO_TEST_FILES ${NO_TEST_FILES} SMCElectroPneumaticItw Twisting)
  endif()

  IF(HAS_BULLET)
    LIST(APPEND EXAMPLES_DIRECTORIES Mechanics/ContactDetection)

    if(HAS_PYTHON_WRAPPER)
      find_python_module(h5py)
      if(NOT h5py_FOUND)
        LIST(APPEND NO_TEST_FILES "BulletIO")
      endif(NOT h5py_FOUND)
    else(HAS_PYTHON_WRAPPER)
      LIST(APPEND NO_TEST_FILES "BulletIO")
    endif(HAS_PYTHON_WRAPPER)
ENDIF(HAS_BULLET)

  # Scan each dir in examples_directories and collect source files
  foreach(_D ${EXAMPLES_DIRECTORIES})
    message("Scan example directory ${_D}")
    set(dirlist)
    get_subdirectories(dirlist ${CMAKE_CURRENT_SOURCE_DIR}/${_D})
    foreach(_dir ${dirlist})
      set(src_dir ${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_dir})
      file(GLOB EXAMPLES_P ${src_dir}/*.cpp ${src_dir}/*.py)
      file(GLOB EXAMPLES_DATA ${src_dir}/*.ref)
      if(EXAMPLES_P)
        foreach(_P ${EXAMPLES_P})
	  # Full path to current file ...
	  get_filename_component(_RP ${_P} ABSOLUTE)
	  if(NOT _RP MATCHES Plugin) # To exclude plugins
	    set(TEST_ME TRUE)
	    # Exclude 'no test' files
	    foreach(_NT ${NO_TEST_FILES})
	      if(_RP MATCHES "${_NT}")
	        set(TEST_ME FALSE)
	      endif()
	    endforeach()
	    if(TEST_ME)
	      get_filename_component(EXAMPLE_DIR ${_RP} PATH)
	      get_filename_component(EXAMPLE_NAME ${_RP} NAME_WE)
	      get_filename_component(ext ${_P} EXT)
	      message("Add ${_D}/${_dir}/${EXAMPLE_NAME}")
	      set(bin_dir ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_dir}/${EXAMPLE_NAME})
	      file(MAKE_DIRECTORY ${bin_dir})
	      foreach(datafile ${EXAMPLES_DATA})
	        configure_file(${datafile} ${bin_dir})
	      endforeach()
	      if(ext MATCHES ".py")
 	        set(target_name ${_dir}_${EXAMPLE_NAME}_py)
 	      else()
	        set(target_name ${_dir}_${EXAMPLE_NAME})
	        add_custom_target(${target_name}_build ALL ${CMAKE_COMMAND} -E chdir ${bin_dir} ${SICONOS_SCRIPT} --noexec --silent ${_P}
		  WORKING_DIRECTORY ${bin_dir})
	      endif()
 	      add_test(${target_name} ${CMAKE_COMMAND} -E chdir ${bin_dir} ${SICONOS_SCRIPT} ${_P})
 	    endif()
	  endif()
        endforeach()
      endif()
    endforeach()
  endforeach()
endif()    

# The following should work, but isn't ... so I used NO_TEST_FILES -- xhub
#SET_TESTS_PROPERTIES("SliderCrankD1MinusLinear" PROPERTIES WILL_FAIL TRUE)
