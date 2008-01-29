# ########## Project setup ##########
PROJECT(Numerics)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.4)

# ######### General setup ##########
INCLUDE_DIRECTORIES(utils)

SUBDIRS(utils NonSmoothSolvers netlib)

# Headers:
SET(src_library_HDRS
    Solver_Options.h
    SiconosNumerics.h
)

ADD_LIBRARY(Numerics ${src_library_HDRS})
SET_TARGET_PROPERTIES(Numerics PROPERTIES VERSION 0.0.0 LINKER_LANGUAGE C)

TARGET_LINK_LIBRARIES(Numerics NonSmoothSolvers/NonSmoothSolvers.a)
