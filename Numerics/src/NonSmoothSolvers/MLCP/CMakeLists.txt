# ########## Project setup ##########
PROJECT(MLCP)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## MLCP library ##########
# Sources:
SET(MLCP_library_SRCS
    mlcp_driver.c
    mlcp_pgs.c
    mlcp_filter_result.c
    mlcp_compute_error.c
    mlcp_path.c
    mlcp_rpsor.c
    mlcp_rpgs.c
    mlcp_psor.c
)

# Headers:
SET(MLCP_library_HDRS
    MLCP_Solvers.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(MLCP ${MLCP_library}_SCRS ${MLCP_library_HDRS})
SET_TARGET_PROPERTIES(MLCP PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS MLCP DESTINATION lib)
