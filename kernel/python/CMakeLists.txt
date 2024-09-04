set(kernel_root_dir ${CMAKE_SOURCE_DIR}/kernel/src)

pybind11_add_module(modeling 
  ${CMAKE_CURRENT_SOURCE_DIR}/modeling/modeling.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/modeling/dynamical_systems_wrapper.cpp
)

target_link_libraries(modeling PUBLIC kernel numerics)

pybind11_add_module(integrators 
  ${CMAKE_CURRENT_SOURCE_DIR}/integrators/integrators.cpp
)

target_link_libraries(integrators PUBLIC kernel numerics)

pybind11_add_module(simulation
  ${CMAKE_CURRENT_SOURCE_DIR}/simulation/simulation.cpp
)

target_link_libraries(simulation PUBLIC kernel numerics)

pybind11_add_module(nonsmooth_formulations
  ${CMAKE_CURRENT_SOURCE_DIR}/nonsmooth_formulations/nonsmooth_formulations.cpp
)

target_link_libraries(nonsmooth_formulations PUBLIC kernel numerics)

install(TARGETS
  modeling
  integrators
  simulation
  nonsmooth_formulations
  COMPONENT python
  LIBRARY DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/"
  )
