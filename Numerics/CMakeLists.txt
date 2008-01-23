# ########## Project setup ##########
PROJECT(Siconos-Numerics)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# Sub directories 
SUBDIRS(src)
INCLUDE_DIRECTORIES(src)

# ########## src library ##########
# Sources:
SET(src_library_SRCS
    dummy.cpp
)

