# ########## Project setup ##########
PROJECT(netlib)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## netlib library ##########
# Sources:
SET(netlib_library_SRCS
)

# Headers:
SET(netlib_library_HDRS
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(netlib ${netlib_library}_SCRS ${netlib_library_HDRS})
SET_TARGET_PROPERTIES(netlib PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS netlib DESTINATION lib)
