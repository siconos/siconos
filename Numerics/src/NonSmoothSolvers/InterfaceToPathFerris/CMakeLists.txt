# ########## Project setup ##########
PROJECT(InterfaceToPathFerris)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

INCLUDE_DIRECTORIES(include)
# ########## InterfaceToPathFerris library ##########
# Sources:
SET(InterfaceToPathFerris_library_SRCS
    SimpleLCP_Path.c
)

# Headers:
SET(InterfaceToPathFerris_library_HDRS
    SimpleLCP.h
)

# actual target:
ADD_LIBRARY(InterfaceToPathFerris ${InterfaceToPathFerris_library_HDRS} ${InterfaceToPathFerris_library_SRCS})
SET_TARGET_PROPERTIES(InterfaceToPathFerris PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS InterfaceToPathFerris DESTINATION lib)
