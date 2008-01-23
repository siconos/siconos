# ########## Project setup ##########
PROJECT(src)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## src library ##########
# Sources:
SET(src_library_SRCS
    dummy.cpp
)

# Headers:
SET(src_library_HDRS
    Solver_Options.h
    SiconosNumerics.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(src ${src_library}_SCRS ${src_library_HDRS})
SET_TARGET_PROPERTIES(src PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS src DESTINATION lib)
