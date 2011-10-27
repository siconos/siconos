CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)

INCLUDE(SiconosProject)

SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Front-End provides Python bindings to the Siconos library.")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/Front-End provides Python bindings.")

SET(DEBIAN_PACKAGE_DEPENDS 
  "g++, python, libboost-dev")

SICONOS_PROJECT(FrontEnd 3 4 0)

# To have a full print of variables, use -DPRINT_ENV=ON
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "${PROJECT_NAME} version ${VERSION} is now ready for compulation and installation.")
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
