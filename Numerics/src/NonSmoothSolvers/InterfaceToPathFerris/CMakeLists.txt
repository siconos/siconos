# ########## Project setup ##########
PROJECT(InterfaceToPathFerris)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# ########## InterfaceToPathFerris library ##########
# Sources:
SET(InterfaceToPathFerris_library_SRCS
    SimpleLCP_Path.c
)

# Headers:
SET(InterfaceToPathFerris_library_HDRS
    SimpleLCP.h
)

CONFIGURE_FILE(Makefile.in ${CMAKE_BINARY_DIR}/Makefile)


# actual target:
ADD_LIBRARY(InterfaceToPathFerris ${InterfaceToPathFerris_library}_SCRS ${InterfaceToPathFerris_library_HDRS})
SET_TARGET_PROPERTIES(InterfaceToPathFerris PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS InterfaceToPathFerris DESTINATION lib)
