# ########## Project setup ##########
PROJECT(NonSmoothSolvers)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.4)

# ######### General setup ##########

SET(SOLVDIRS
  NSSTools
  Relay
  LCP
  MLCP
  QP
  NCP
  InterfaceToPathFerris
  FrictionContact2D
  FrictionContact3D
  pfc3D)


INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${SOLVDIRS}) 
  

ADD_DEFINITIONS(-std=c99)

SUBDIRS(${SOLVDIRS})

# ########## NonSmoothSolvers library ##########
# Sources:
SET(NonSmoothSolvers_library_SRCS
    NonSmoothNewton.c
)

# Headers:
SET(NonSmoothSolvers_library_HDRS
    NonSmoothDrivers.h
    NonSmoothNewton.h
)


# actual target:
ADD_LIBRARY(NonSmoothSolvers ${NonSmoothSolvers_library_HDRS} ${NonSmoothSolvers_library_SRCS})
SET_TARGET_PROPERTIES(NonSmoothSolvers PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS NonSmoothSolvers DESTINATION lib)
