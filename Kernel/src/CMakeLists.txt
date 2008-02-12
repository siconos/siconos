#
# Here comes Numerics
#
PROJECT(SiconosKernel C CXX)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosKernel_VERSION "${VERSION}")

#
# Published headers
#
# if <PROJECT>_HDRS unset -> all headers under <PROJECT>_DIRS

#
# Project directories 
#
SET(SiconosKernel_DIRS
  .
  controlTools
  model
  modelingTools
  plugin
  simulationTools
  utils
  utils/SiconosAlgebra
  utils/SiconosException
  utils/SiconosMemory
  utils/SiconosSharedLibrary
  utils/SiconosTools
  xml)


#
# Unstable sources
#
SET(SiconosKernel_Unstable_SRCS)

#
# Installation layout
#
SET(SiconosKernel_INCLUDE_DIR include/Siconos/Kernel)

#
# compilation setup 
#
SET(SiconosKernel_LINK ${SiconosNumerics_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LIBXML2_LIBRARIES})

LIBRARY_PROJECT_SETUP()

#ADD_DEFINITIONS(-std=c99) : warning with f95
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIRS} ${SiconosNumerics_INCLUDE_DIRS})
ADD_DEFINITIONS(${LIBXML2_DEFINITIONS})

# test directories
#SUBDIRS(netlib/odepack/test)
#SUBDIRS(NonSmoothSolvers/test)
