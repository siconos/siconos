# ########## Project setup ##########
PROJECT(odepack)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## odepack library ##########
# Sources:
SET(odepack_library_SRCS
)

# Headers:
SET(odepack_library_HDRS
    odepack.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(odepack ${odepack_library}_SCRS ${odepack_library_HDRS})
SET_TARGET_PROPERTIES(odepack PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS odepack DESTINATION lib)
