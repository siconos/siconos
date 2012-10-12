# cf http://stackoverflow.com/questions/435708/any-way-in-cmake-to-require-gcc-version-4
# with cmake version >= 2.8.8 we can use CMAKE_C_COMPILER_VERSION and
# CMAKE_CXX_COMPILER_VERSION
IF(CMAKE_COMPILER_IS_GNUCXX)
  EXEC_PROGRAM(
    ${CMAKE_CXX_COMPILER}
    ARGS                    --version
    OUTPUT_VARIABLE _compiler_output)
  STRING(REGEX REPLACE ".*([0-9]\\.[0-9]\\.[0-9]).*" "\\1"
    gcc_compiler_version ${_compiler_output})
  MESSAGE(STATUS "C++ compiler version: ${gcc_compiler_version} [${CMAKE_CXX_COMPILER}]")
  
  IF(gcc_compiler_version VERSION_LESS "4.5")
    MESSAGE(FATAL_ERROR "gcc greater than 4.5 needed")
  ENDIF()
ENDIF()

FIND_PACKAGE(SWIG 2.0.7 REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonInterp)

FIND_PACKAGE(PythonLibs)

MESSAGE("PYTHON EXECUTABLE : ${PYTHON_EXECUTABLE}")
MESSAGE("Python_ADDITIONAL_VERSIONS : ${Python_ADDITIONAL_VERSIONS}")
MESSAGE("Python_FRAMEWORKS : ${Python_FRAMEWORKS}")

MESSAGE("PYTHON INCLUDE DIR : ${PYTHON_INCLUDE_DIR}")
MESSAGE("PYTHON LIBRARIES : ${PYTHON_LIBRARIES}")

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

COMPILE_WITH(Python_Numpy REQUIRED)

# Numerics + Kernel 
COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)

COMPILE_WITH(SiconosIO)

IF(SiconosIO_FOUND)
  SET(HAVE_SICONOS_IO TRUE)
ENDIF(SiconosIO_FOUND)

COMPILE_WITH(Boost REQUIRED)
FIND_PACKAGE(Boost COMPONENTS graph)

IF(NOT Boost_GRAPH_FOUND)
  MESSAGE(FATAL_ERROR "boost graph library is needed.")
ENDIF(NOT Boost_GRAPH_FOUND)

COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(GMP REQUIRED)
COMPILE_WITH(ATLAS)

COMPILE_WITH(LAPACK REQUIRED COMPLETE)
COMPILE_WITH(BLAS REQUIRED)

# contact detection
IF(WITH_BULLET)
  COMPILE_WITH(Bullet REQUIRED)
  IF(BULLET_FOUND)
    SET(HAVE_BULLET TRUE)
    MESSAGE( STATUS " Bullet include dirs : ${BULLET_INCLUDE_DIRS}" )
    SET(SiconosKernel_DIRS ${SiconosKernel_DIRS} multiBodyTools/contactDetection/bullet)
    IF(BULLET_USE_DOUBLE_PRECISION)
      APPEND_FLAGS(CMAKE_CXX_FLAGS "-DBT_USE_DOUBLE_PRECISION")
    ENDIF(BULLET_USE_DOUBLE_PRECISION)
  ENDIF(BULLET_FOUND)
ENDIF(WITH_BULLET)

#SET(CMAKE_SWIG_FLAGS "")

IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
  SET(EMULATOR "wine ")
ELSE()
  SET(EMULATOR)
ENDIF()

#set correct lib install path on linux
IF(CMAKE_SYSTEM_NAME MATCHES Linux)
  IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(LIB_SUFFIX "")
  ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(LIB_SUFFIX "64")
  ENDIF()
ENDIF()

# For some environment variables (LD_LIBRARY_PATH, DYLD_LIBRARY_PATH, Path)
#GET_FILENAME_COMPONENT(SICONOS_NUMERICS_PATH SiconosNumerics_FOUND PATH)
SET(SICONOS_NUMERICS_PATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

CONFIGURE_FILE(setup.py.in setup.py)

CONFIGURE_FILE(runtest.py.in ${CMAKE_BINARY_DIR}/src/swig/runtest.py)

find_python_module(pytest REQUIRED)

SUBDIRS(Siconos)

ADD_CUSTOM_TARGET(distutils-install
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py build
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py install
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Installing with distutils")

# config.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/FrontEndConfig.h DESTINATION include/Siconos/FrontEnd)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/Siconos/KernelRegistration.i DESTINATION include/Siconos/FrontEnd)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/Siconos/Common.i DESTINATION include/Siconos/FrontEnd)
INSTALL(FILES ${CMAKE_BINARY_DIR}/src/swig/Siconos/Kernel-docstrings.i DESTINATION include/Siconos/FrontEnd)


# To have a full print of variables, use -DPRINT_ENV=ON
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "${PROJECT_NAME} version ${VERSION} is now ready for compilation and installation.")
message(STATUS "To proceed run 'make' and 'make install' and optionaly 'make test'.")
message(STATUS " C++ Compiler : ${CMAKE_CXX_COMPILER}")
message(STATUS " C Compiler : ${CMAKE_C_COMPILER}")
message(STATUS " Compilation mode is : ${CMAKE_BUILD_TYPE}")
message(STATUS " Code Sources are in : ${CMAKE_SOURCE_DIR}")
message(STATUS " The project will be installed in ${CMAKE_INSTALL_PREFIX}")
message(STATUS " The Numerics library used is : ${SiconosNumerics_FOUND}")
message(STATUS " The Kernel library used is : ${SiconosKernel_FOUND}")
message(STATUS " The Python libraries are : ${PYTHON_LIBRARIES}")
message(STATUS " The Python includes are : ${PYTHON_INCLUDE_DIR}")
message(STATUS " To get more information about dependencies, config or else, ")
message(STATUS "check CMakeCache.txt file or re-run cmake with -DPRINT_ENV=ON.")
message(STATUS "=================================================================================================\n")
