#================================================================
# cmake utilities to build python-bindings for mechanics component
#================================================================

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES
  "mechanics/joints.i"
  "mechanics/collision/base.i"
  "mechanics/collision/native.i"
  "mechanics/collision/bodies.i")

find_python_module(h5py)
find_python_module(vtk REQUIRED) # Required for native (if not present, python_test_native_collision will fail)
  
configure_file(mechanics/collision/tools.py ${SICONOS_SWIG_ROOT_DIR}/mechanics/collision/tools.py)
configure_file(mechanics/collision/convexhull.py  ${SICONOS_SWIG_ROOT_DIR}/mechanics/collision/convexhull.py )
configure_file(mechanics/collision/__init__.py  ${SICONOS_SWIG_ROOT_DIR}/mechanics/collision/__init__.py )
include(swig_python_tools)

foreach(fname IN LISTS ${COMPONENT}_PYTHON_MODULES)
  add_swig_sub_module(
    FILE ${fname}
    DEPS mechanics
    INCLUDES ${CMAKE_SOURCE_DIR}/kernel/swig
    #COMPILE_OPTIONS ${swig_mechanics_options}
    )
endforeach()

if(SICONOS_HAS_BULLET)
  add_swig_sub_module(
    FILE mechanics/collision/bullet.i
    DEPS "mechanics;BULLET::BULLET"
    INCLUDES ${CMAKE_SOURCE_DIR}/kernel/swig
    COMPILE_OPTIONS "-w312,322,325,350,351,362,383,389,394,395,401,402,403,503,512,520"
    COMPILE_DEFINITIONS ${BULLET_DEFINITIONS}
    )
endif()

if(WITH_OCE)
  add_swig_sub_module(
    FILE mechanics/occ.i
    DEPS mechanics
    INCLUDES ${CMAKE_SOURCE_DIR}/kernel/swig
    )
endif()

if(WITH_FREECAD)
  # -- FreeCAD --
  # For python bindings.
  find_package(FeedCAD COMPONENTS Part REQUIRED)
  #  FP : which swig module needs freecad? Is this part outdated? 
  # target_link_libraries(?? PRIVATE ${FreeCAD_LIBRARIES})
endif()


# --- Tests ---
if(WITH_TESTING)
  if(NOT WITH_OCE)
    list(APPEND python_excluded_tests tests/test_occ_inertia.py)
  else()
    find_python_module(OCC REQUIRED)
  endif()
  build_python_tests(
    # DEPS ${COMPONENT} # some plugins in tests need to be linked with kernel
    EXCLUDE ${python_excluded_tests})
endif()

if(NOT NO_RUNTIME_BUILD_DEP)
  find_python_module(pyhull)
endif()
 
