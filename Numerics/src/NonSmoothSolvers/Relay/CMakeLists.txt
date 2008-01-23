# ########## Project setup ##########
PROJECT(Relay)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## Relay library ##########
# Sources:
SET(Relay_library_SRCS
    pr_nlgs.c
    pr_driver.c
    dr_nlgs.c
    pr_latin.c
    pr_gsnl.c
    dr_latin.c
    dr_driver.c
)

# Headers:
SET(Relay_library_HDRS
    Relay_Solvers.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(Relay ${Relay_library}_SCRS ${Relay_library_HDRS})
SET_TARGET_PROPERTIES(Relay PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS Relay DESTINATION lib)
