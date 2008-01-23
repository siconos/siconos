# ########## Project setup ##########
PROJECT(FrictionContact3D)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## FrictionContact3D library ##########
# Sources:
SET(FrictionContact3D_library_SRCS
    FrictionContact3D_AlartCurnier.c
    FrictionContact3D_Newton.c
    FrictionContact3D_driver.c
    FrictionContact3D2NCP_Glocker.c
    FrictionContact3D_projection.c
    FrictionContact3D_solver_SBS.c
    FrictionContact3D_nsgs_SBS.c
    FrictionContact3D_nsgs.c
)

# Headers:
SET(FrictionContact3D_library_HDRS
    FrictionContact3D_AlartCurnier.h
    FrictionContact3D_Solvers.h
    FrictionContact3D_Newton.h
    FrictionContact3D_projection.h
    FrictionContact3D_Problem.h
    FrictionContact3D2NCP_Glocker.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(FrictionContact3D ${FrictionContact3D_library}_SCRS ${FrictionContact3D_library_HDRS})
SET_TARGET_PROPERTIES(FrictionContact3D PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS FrictionContact3D DESTINATION lib)
