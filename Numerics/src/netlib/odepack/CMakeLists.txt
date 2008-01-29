# ########## Project setup ##########
PROJECT(odepack)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# ########## odepack library ##########
# Sources:
SET(odepack_library_SRCS
  opkda1.f
  opkda2.f
  opkdmain.f
)

# Headers:
SET(odepack_library_HDRS
    odepack.h
)

#F95 Bug
CONFIGURE_FILE(odepack.h DPCG.mod.proxy COPYONLY) 
CONFIGURE_FILE(odepack.h DPCGS.mod.proxy COPYONLY) 
CONFIGURE_FILE(odepack.h DSOLBT.mod.proxy COPYONLY) 
CONFIGURE_FILE(odepack.h DUSOL.mod.proxy COPYONLY) 
CONFIGURE_FILE(odepack.h if.mod.proxy COPYONLY) 
CONFIGURE_FILE(odepack.h in.mod.proxy COPYONLY) 
CONFIGURE_FILE(odepack.h the.mod.proxy COPYONLY) 

# actual target:
ADD_LIBRARY(odepack ${odepack_library_HDRS} ${odepack_library_SRCS})
SET_TARGET_PROPERTIES(odepack PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS odepack DESTINATION lib)
