# ########## Project setup ##########
PROJECT(NCP)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## NCP library ##########
# Sources:
SET(NCP_library_SRCS
    NCP_compute_error.c
    NCP_FischerBurmeister.c
    NCP_block_compute_error.c
)

# Headers:
SET(NCP_library_HDRS
    NCP.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(NCP ${NCP_library}_SCRS ${NCP_library_HDRS})
SET_TARGET_PROPERTIES(NCP PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS NCP DESTINATION lib)
