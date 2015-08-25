INCLUDE(FindSiconosGeneric)
FIND_SICONOS_COMPONENT(Kernel SiconosKernel.hpp)

FOREACH(_DIR ${SiconosKernel_LOCATION_DIRS})
  FIND_PATH(SiconosKernel_EXE_DIR bin/siconos
    HINTS ${_DIR}
    ENV PATH)
ENDFOREACH()


IF(SiconosKernel_EXE_DIR)
  SET(SiconosKernel_EXE_DIR "${SiconosKernel_EXE_DIR}/bin")
  MESSAGE(STATUS "siconos executable found in ${SiconosKernel_EXE_DIR}")
ELSE(SiconosKernel_EXE_DIR)
  IF(SiconosKernel_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required siconos executable not found!")
  ENDIF(SiconosKernel_FIND_REQUIRED)
ENDIF(SiconosKernel_EXE_DIR)

FOREACH(_DIR ${SiconosKernel_LOCATION_DIRS})
  FIND_PATH(SiconosKernel_CMAKE_DIR share/siconos-kernel/cmake
    HINTS ${_DIR}
    ENV PATH)
ENDFOREACH()

IF(SiconosKernel_CMAKE_DIR)
  SET(SiconosKernel_CMAKE_DIR "${SiconosKernel_CMAKE_DIR}/share/siconos-kernel/cmake")
ENDIF(SiconosKernel_CMAKE_DIR)


