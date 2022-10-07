set(kernel_root_dir ${CMAKE_SOURCE_DIR}/kernel/src)

pybind11_add_module(algebra 
  ${CMAKE_CURRENT_SOURCE_DIR}/algebra/algebra.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/algebra/vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/algebra/matrix.cpp
)
target_link_libraries(algebra PUBLIC kernel numerics)

pybind11_add_module(modeling 
  ${CMAKE_CURRENT_SOURCE_DIR}/modeling/modeling.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/modeling/dynamical_systems.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/modeling/nonsmoothlaws.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/modeling/relations.cpp
)

target_link_libraries(modeling PUBLIC kernel numerics)

pybind11_add_module(integrators 
  ${CMAKE_CURRENT_SOURCE_DIR}/integrators/integrators.cpp
)

target_link_libraries(integrators PUBLIC kernel numerics)


install(TARGETS
  algebra
  modeling
  integrators
  COMPONENT python
  LIBRARY DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/"
  )
