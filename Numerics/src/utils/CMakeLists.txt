# ########## Project setup ##########
PROJECT(utils)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

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

# actual target:
ADD_LIBRARY(utils  ${utils_library_HDRS} ${utils_library_SRCS})
SET_TARGET_PROPERTIES(utils PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS utils DESTINATION lib)
