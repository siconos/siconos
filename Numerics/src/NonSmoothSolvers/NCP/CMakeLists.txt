# ########## Project setup ##########
PROJECT(NCP)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)


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

# actual target:
ADD_LIBRARY(NCP ${NCP_library_HDRS} ${NCP_library_SRCS} )
SET_TARGET_PROPERTIES(NCP PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS NCP DESTINATION lib)
