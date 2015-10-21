IF(MSVC)
  SET(CMAKE_CXX_CREATE_SHARED_MODULE
"<CMAKE_LINKER> ${CMAKE_CL_NOLOGO} <OBJECTS> ${CMAKE_START_TEMP_FILE} /out:<TARGET> /pdb:<TARGET_PDB> /dll /version:<TARGET_VERSION_MAJOR>.<TARGET_VERSION_MINOR> <LINK_FLAGS> <LINK_LIBRARIES> ${CMAKE_END_TEMP_FILE}")
  # we need CMAKE_NM for the EXPORTS hack
  # On other Module this variable is set through FindBinutils since we ask for
  # a Fortran compiler --xhub
  find_program(CMAKE_NM NAMES ${_CMAKE_TOOLCHAIN_PREFIX}nm HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
  include(Platform/Windows-GNU) # for proper prefixes and suffixes
ENDIF(MSVC)

IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
  SET(EXTRA_EXT ".a")
ELSE()
  SET(EXTRA_EXT)
ENDIF()

