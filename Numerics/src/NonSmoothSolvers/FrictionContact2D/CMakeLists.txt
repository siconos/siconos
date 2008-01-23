# ########## Project setup ##########
PROJECT(FrictionContact2D)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## FrictionContact2D library ##########
# Sources:
SET(FrictionContact2D_library_SRCS
    pfc_2D_cpg.c
    dfc_2D2cond_2D.c
    dfc_2D2lcp.c
    dfc_2D_driver.c
    lcp2dfc_2D.c
    pfc_2D_driver.c
    projf.c
    pfc_2D_nlgs.c
    pfc_2D_latin.c
    cond_2D2dfc_2D.c
    pfc_2D_projc.c
    pfc_2D_projf.c
    dfc_2D_latin.c
    projc.c
)

# Headers:
SET(FrictionContact2D_library_HDRS
    FrictionContact2D_Solvers.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(FrictionContact2D ${FrictionContact2D_library}_SCRS ${FrictionContact2D_library_HDRS})
SET_TARGET_PROPERTIES(FrictionContact2D PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS FrictionContact2D DESTINATION lib)
