# Atlas and boost headers
FILE(GLOB ATLAS_HDRS AtlasLocal/atlas/*.h)
FILE(GLOB BOOST_HDRS boost/numeric/bindings/*/*.hpp)

INSTALL(FILES ${ATLAS_HDRS};${BOOST_HDRS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Siconos/Kernel)
