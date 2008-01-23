# ########## Project setup ##########
PROJECT(QP)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## QP library ##########
# Sources:
SET(QP_library_SRCS
)

# Headers:
SET(QP_library_HDRS
    QP_Solvers.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(QP ${QP_library}_SCRS ${QP_library_HDRS})
SET_TARGET_PROPERTIES(QP PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS QP DESTINATION lib)
