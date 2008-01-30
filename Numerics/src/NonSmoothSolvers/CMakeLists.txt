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
ADD_DEPENDENCIES(NonSmoothSolvers ${SOLVDIRS})
