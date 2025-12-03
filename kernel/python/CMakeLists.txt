set(kernel_root_dir ${CMAKE_SOURCE_DIR}/kernel/src)


create_pybind11_module(
  MODULE_NAME pb11_template
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/pb11_template.cpp
  DEPENDENCIES kernel numerics
)
create_pybind11_module(
  MODULE_NAME geometry 
  SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry_pybind11.cpp
  DEPENDENCIES kernel numerics
)
target_include_directories(geometry PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

create_pybind11_module(
  MODULE_NAME modeling 
  SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen2python_pybind11.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/modeling/modeling_pybind11.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/modeling/dynamical_systems_pybind11.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/modeling/relations_pybind11.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/modeling/nonsmoothlaws_pybind11.cpp
  DEPENDENCIES kernel numerics
)
target_include_directories(modeling PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

create_pybind11_module(MODULE_NAME  integrators 
SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/integrators/integrators_pybind11.cpp
DEPENDENCIES kernel numerics
)

create_pybind11_module(MODULE_NAME simulation
SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/simulation/simulation_pybind11.cpp
DEPENDENCIES kernel numerics
)

create_pybind11_module(MODULE_NAME nonsmooth_formulations
SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/nonsmooth_formulations/nonsmooth_formulations_pybind11.cpp
DEPENDENCIES kernel numerics
)

create_pybind11_module(MODULE_NAME input
SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/io/io_pybind11.cpp
DEPENDENCIES kernel numerics
)

add_custom_target(py${COMPONENT}
    DEPENDS pb11_template
    geometry
    modeling
    integrators
    simulation
    nonsmooth_formulations
    input
  ) 

 

install(TARGETS
  pb11_template
  geometry
  modeling
  integrators
  simulation
  nonsmooth_formulations
  input
  COMPONENT python
  LIBRARY DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/"
  )

# --- Tests ---
if(WITH_TESTING)
  build_python_tests()
endif()
