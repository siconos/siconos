SICONOS_SWIG_MODULE(Kernel "Kernel;")

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../KernelTypes.i
              ${CMAKE_CURRENT_SOURCE_DIR}/Kernel.i
              ${CMAKE_CURRENT_SOURCE_DIR}/addons.hpp
  DESTINATION share/siconos-front-end/swig)


