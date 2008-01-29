# ########## Project setup ##########
PROJECT(QP)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# ########## QP library ##########
# Sources:
SET(QP_library_SRCS
  ql0001.f
)

# Headers:
SET(QP_library_HDRS
    QP_Solvers.h
)


# actual target:
ADD_LIBRARY(QP ${QP_library_HDRS} ${QP_library_SRCS})
SET_TARGET_PROPERTIES(QP PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS QP DESTINATION lib)
