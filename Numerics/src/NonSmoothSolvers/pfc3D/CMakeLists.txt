# ########## Project setup ##########
PROJECT(pfc3D)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## pfc3D library ##########
# Sources:
SET(pfc3D_library_SRCS
    pfc_3D_nlgs.c
    pfc_3D_nsgs.c
    pfc_3D_Alart_Curnier.c
    pfc_3D_projection.c
    filter_result_pfc.c
    pfc_3D_solver.c
    pfc_3D_Fischer_Burmeister.c
    pfc_3D_cpg.c
    pfc_3D_nsgs_new.c
    pfc_3D_solver_block.c
    pfc_3D_projc.c
    pfc_3D_nlgsnewton.c
    pfc_3D_nsgs_block.c
    pfc_3D_newton.c
    filter_result_pfc_block2.c
    filter_result_pfc_block.c
    pfc_3D_projf.c
)

# Headers:
SET(pfc3D_library_HDRS
    pfc_3D_Solvers.h
    pfc_3D_Fischer_Burmeister.h
    pfc_3D_Alart_Curnier.h
)


# actual target:
ADD_LIBRARY(pfc3D ${pfc3D_library_HDRS} ${pfc3D_library_SRCS})
SET_TARGET_PROPERTIES(pfc3D PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS pfc3D DESTINATION lib)
