# ########## Project setup ##########
PROJECT(LCP)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## LCP library ##########
# Sources:
SET(LCP_library_SRCS
    lcp_newton_FB.c
    filter_result_LCP_block.c
    lcp_nsqp.c
    lcp_cpg.c
    lcp_newton_min.c
    lcp_latin.c
    filter_result_LCP.c
    lcp_path.c
    lcp_pgs.c
    lcp_driver.c
    lcp_qp.c
    lcp_latin_w.c
    lcp_lexicolemke.c
    lcp_compute_error.c
    lcp_psor.c
    lcp_rpgs.c
    lcp_solver_block.c
)

# Headers:
SET(LCP_library_HDRS
    LCP_Solvers.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(LCP ${LCP_library}_SCRS ${LCP_library_HDRS})
SET_TARGET_PROPERTIES(LCP PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS LCP DESTINATION lib)
