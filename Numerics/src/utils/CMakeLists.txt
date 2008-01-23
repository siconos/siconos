# ########## Project setup ##########
PROJECT(utils)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## utils library ##########
# Sources:
SET(utils_library_SRCS
    Numerics_Options.c
)

# Headers:
SET(utils_library_HDRS
    Numerics_Options.h
    LA.h
    blaslapack.h
    f2c.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(utils ${utils_library}_SCRS ${utils_library_HDRS})
SET_TARGET_PROPERTIES(utils PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS utils DESTINATION lib)
