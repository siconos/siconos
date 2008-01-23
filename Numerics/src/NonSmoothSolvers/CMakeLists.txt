# ########## Project setup ##########
PROJECT(NonSmoothSolvers)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

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

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(NonSmoothSolvers ${NonSmoothSolvers_library}_SCRS ${NonSmoothSolvers_library_HDRS})
SET_TARGET_PROPERTIES(NonSmoothSolvers PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS NonSmoothSolvers DESTINATION lib)
