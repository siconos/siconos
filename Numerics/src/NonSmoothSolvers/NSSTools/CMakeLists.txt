# ########## Project setup ##########
PROJECT(NSSTools)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## NSSTools library ##########
# Sources:
SET(NSSTools_library_SRCS
    abs_part.c
    max_part.c
    maxf.c
    minf.c
    diffns.c
    min_part.c
    freeSpBlMat.c
    clapack_dtrtrs.c
    Solver_Options.c
    pos_part.c
    NumericsMatrix.c
    SparseBlockMatrix.c
)

# Headers:
SET(NSSTools_library_HDRS
    Solver_Options.h
    NumericsMatrix.h
    SparseBlockMatrix.h
    NSSTools.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(NSSTools ${NSSTools_library}_SCRS ${NSSTools_library_HDRS})
SET_TARGET_PROPERTIES(NSSTools PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS NSSTools DESTINATION lib)
