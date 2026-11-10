include(tools4tests)

if(WITH_TESTING)
  add_custom_target(vkernel-tests echo "Start vkernel tests")

  begin_tests(src/siconos/simul/tests DEPS
              "numerics;kernel;mechanics;CPPUNIT::CPPUNIT")
  new_test(SOURCES RtOsiTest.cpp ${SIMPLE_TEST_MAIN})
  # end_tests()

  begin_tests(src/siconos/storage/tests DEPS "vkernel;CPPUNIT::CPPUNIT")
  new_test(SOURCES StorageTest.cpp ${SIMPLE_TEST_MAIN})
  # the schedules (storage/schedule.hpp) and the Moreau-Jean time step as text
  new_test(NAME vkernel_schedule SOURCES ScheduleTest.cpp ${SIMPLE_TEST_MAIN} DEPS
           "numerics;kernel;mechanics")
  new_test(SOURCES t1.cpp ${SIMPLE_TEST_MAIN} DEPS "numerics;kernel;mechanics")
  # end_tests()

  if(WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()

  # Tests that produce a trajectory also compare it to a committed reference
  # file (tests/ref_check.hpp) at a relative Frobenius error of 1e-8.
  # bouncing_ball's reference is the tutorial's BouncingBallTS.ref (external
  # truth); the others lock the current behaviour against regression.
  begin_tests(src/siconos/tests DEPS
              "vkernel;numerics;CPPUNIT::CPPUNIT;Eigen3::Eigen;pthread")
  new_test(
    NAME
    vkernel_bouncing_ball
    SOURCES
    bouncing_ball.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_ball.ref)
  new_test(
    NAME
    vkernel_bouncing_ball_registry
    SOURCES
    bouncing_ball_registry.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_ball_registry.ref)
  new_test(
    NAME
    vkernel_bouncing_ball_split
    SOURCES
    bouncing_ball_split.cpp
    bouncing_ball_split_sim.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_ball_split.ref)
  new_test(
    NAME
    vkernel_unbounded_bouncing_ball
    SOURCES
    unbounded_bouncing_ball.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    unbounded_bouncing_ball.ref)
  new_test(
    NAME
    vkernel_bouncing_balls
    SOURCES
    bouncing_balls.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_balls.ref)
  new_test(
    NAME
    vkernel_bouncing_disk
    SOURCES
    bouncing_disk.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_disk.ref
    DEPS
    "${COMPACTNSEARCH}")
  new_test(
    NAME
    vkernel_bouncing_disk_line
    SOURCES
    bouncing_disk_line.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_disk_line.ref
    DEPS
    "${COMPACTNSEARCH}")
  new_test(NAME vkernel_fremond_disk SOURCES fremond_disk.cpp ${SIMPLE_TEST_MAIN}
           DEPS "${COMPACTNSEARCH}")
  new_test(NAME vkernel_fremond_rigid SOURCES fremond_rigid.cpp ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_registry_remove SOURCES registry_remove.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_nslaws_multi_group SOURCES nslaws_multi_group.cpp
           ${SIMPLE_TEST_MAIN} DEPS "${COMPACTNSEARCH}")
  new_test(NAME vkernel_polydisperse_disks SOURCES polydisperse_disks.cpp
           ${SIMPLE_TEST_MAIN} DEPS "${COMPACTNSEARCH}")

  begin_tests(
    src/siconos/collision/tests DEPS
    "vkernel;numerics;CPPUNIT::CPPUNIT;Eigen3::Eigen;pthread;${COMPACTNSEARCH}")
  new_test(
    NAME
    vkernel_bouncing_disks
    SOURCES
    bouncing_disks.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_disks.ref)
  new_test(
    NAME
    vkernel_fixed_disks
    SOURCES
    fixed_disks.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    fixed_disks.ref)
  new_test(NAME vkernel_cxd_batch_smoke SOURCES cxd_batch_smoke.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytope_cxd SOURCES polytope_cxd.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_bouncing_polytopes SOURCES bouncing_polytopes.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_churn_polytopes SOURCES churn_polytopes.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytope_fsegment SOURCES polytope_fsegment.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytopemesh SOURCES polytopemesh.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytopefdisk SOURCES polytopefdisk.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytopedisk SOURCES polytopedisk.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_rigid_shapes SOURCES rigid_shapes.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_mixed_bodies SOURCES mixed_bodies.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytope_on_fem SOURCES polytope_on_fem.cpp
           ${SIMPLE_TEST_MAIN})
  # the same with the Moreau-Jean scheme on the global problem
  new_test(NAME vkernel_polytope_on_fem_gosi SOURCES polytope_on_fem.cpp
           ${SIMPLE_TEST_MAIN} COMPILE_DEFINITIONS VKERNEL_GOSI)
  new_test(NAME vkernel_io_contact_points SOURCES io_contact_points.cpp
           ${SIMPLE_TEST_MAIN})
  # the analytic cases of the NSCD scheme: sliding block, incline, column
  # (one test per case)
  set(nscd_dem_cases testDemColumn testDemStability testDemSlidingSquare
                     testDemInclineSticks)
  list(TRANSFORM nscd_dem_cases PREPEND NscdCasesTest::)
  set(nscd_cases
      testSlidingSquare
      testInclineSticks
      testInclineSlides
      testColumn
      testAnitescuTasoraFrictionless
      testAnitescuTasoraSeparation
      testAnitescuTasoraWarmStart)
  list(TRANSFORM nscd_cases PREPEND NscdCasesTest::)
  new_test(NAME vkernel_nscd_cases SOURCES nscd_cases.cpp ${SIMPLE_TEST_MAIN}
           CASES ${nscd_cases} ${nscd_dem_cases})
  # the DEM cases on the dem integrator, with the storage of its own scheme
  new_test(NAME vkernel_nscd_cases_dem SOURCES nscd_cases.cpp ${SIMPLE_TEST_MAIN}
           COMPILE_DEFINITIONS NSCD_DEM_ONLY CASES ${nscd_dem_cases})

  # the same deposit run twice with a multithreaded neighbor search must give
  # bitwise identical states
  new_test(NAME vkernel_reproducible_disks SOURCES reproducible_disks.cpp
           ${SIMPLE_TEST_MAIN})
  set_property(TEST vkernel_reproducible_disks APPEND PROPERTY ENVIRONMENT
                                                           OMP_NUM_THREADS=4)

  # topology::sort_interactions (storage::permute) on both storage backends: the
  # registry, and the compile time database (_static_db). The neighbor search
  # order, hence the interactions, must be reproducible (the neighbor lists are
  # sorted, whatever the number of threads), and the assertions of the library
  # headers must stay enabled in a Release build.
  new_test(NAME vkernel_sort_interactions_registry SOURCES
           sort_interactions.cpp ${SIMPLE_TEST_MAIN})
  new_test(
    NAME
    vkernel_sort_interactions_static_db
    SOURCES
    sort_interactions.cpp
    ${SIMPLE_TEST_MAIN}
    COMPILE_DEFINITIONS
    SORT_STATIC_DB)

  # the tests were written against the vkernel/soa flags, which include the
  # C++23 standard they need for std::print
  set(vkernel_ref_tests
      vkernel_bouncing_ball
      vkernel_bouncing_ball_registry
      vkernel_bouncing_ball_split
      vkernel_unbounded_bouncing_ball
      vkernel_bouncing_balls
      vkernel_bouncing_disk
      vkernel_bouncing_disk_line
      vkernel_fremond_disk
      vkernel_fremond_rigid
      vkernel_bouncing_disks
      vkernel_fixed_disks
      vkernel_registry_remove
      vkernel_nslaws_multi_group
      vkernel_polydisperse_disks
      vkernel_cxd_batch_smoke
      vkernel_polytope_cxd
      vkernel_bouncing_polytopes
      vkernel_churn_polytopes
      vkernel_polytope_fsegment
      vkernel_polytopemesh
      vkernel_polytopefdisk
      vkernel_polytopedisk
      vkernel_rigid_shapes
      vkernel_mixed_bodies
      vkernel_reproducible_disks
      vkernel_polytope_on_fem
      vkernel_polytope_on_fem_gosi
      vkernel_io_contact_points
      vkernel_nscd_cases
      vkernel_nscd_cases_dem
      vkernel_sort_interactions_registry
      vkernel_sort_interactions_static_db)
  foreach(t IN LISTS vkernel_ref_tests)
    set_target_properties(${t} PROPERTIES COMPILE_FLAGS "${flags}")
  endforeach()
  set_target_properties(
    vkernel_sort_interactions_registry vkernel_sort_interactions_static_db
    PROPERTIES COMPILE_FLAGS "${flags} -UNDEBUG")

  set_target_properties(RtOsiTest PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")
  set_target_properties(t1 PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")

  # The python scripts of src/siconos/tests, on the build tree: a nonos
  # package next to the build (the module and its python files), the siconos
  # package of the build, its libraries before the install ones (_nonos is
  # built with the install rpath). Each script checks its results and exits
  # with a nonzero status on a failure.
  if(WITH_PYB11_WRAPPER)
    set(nonos_test_path ${CMAKE_CURRENT_BINARY_DIR}/python)
    foreach(f __init__.py bridge.py simoptions.py)
      configure_file(src/siconos/py/nonos/${f} ${nonos_test_path}/nonos/${f} COPYONLY)
    endforeach()
    add_custom_command(
      TARGET _nonos
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_nonos> ${nonos_test_path}/nonos/)

    # new_python_test(NAME <test> SCRIPT <file in src/siconos/tests> [PYTEST]
    #                 [ARGS ...]): the script, or its tests with PYTEST
    function(new_python_test)
      cmake_parse_arguments(pt "PYTEST" "NAME;SCRIPT" "ARGS" ${ARGN})
      set(rundir ${CMAKE_CURRENT_BINARY_DIR}/python-tests/${pt_NAME})
      file(MAKE_DIRECTORY ${rundir})
      set(script ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/tests/${pt_SCRIPT})
      if(pt_PYTEST)
        set(command ${Python_EXECUTABLE} -m pytest -p no:cacheprovider ${script} ${pt_ARGS})
      else()
        set(command ${Python_EXECUTABLE} ${script} ${pt_ARGS})
      endif()
      add_test(NAME ${pt_NAME} COMMAND ${command} WORKING_DIRECTORY ${rundir})
      set(libs
          "$<TARGET_FILE_DIR:numerics>:$<TARGET_FILE_DIR:kernel>:$<TARGET_FILE_DIR:mechanics>:$<TARGET_FILE_DIR:externals>:$<TARGET_FILE_DIR:externals_fortran>"
      )
      # the neighbor search of a GPU build, a library of its own: else the
      # installed one
      if(TARGET cuNSearch)
        string(APPEND libs ":$<TARGET_FILE_DIR:cuNSearch>")
      endif()
      set_tests_properties(
        ${pt_NAME}
        PROPERTIES ENVIRONMENT
                   "PYTHONPATH=${nonos_test_path}:${SICONOS_PB11_BINARY_DIR};LD_LIBRARY_PATH=${libs}:$ENV{LD_LIBRARY_PATH};MPLBACKEND=Agg;MPLCONFIGDIR=${rundir}/mpl;SICONOS_BUILD_DIR=${CMAKE_BINARY_DIR};SICONOS_CHECK_REF=1;PYTHONDONTWRITEBYTECODE=1"
      )
    endfunction()

    new_python_test(NAME vkernel_python_fremond_oblique_disk SCRIPT fremond_oblique_disk.py)
    new_python_test(NAME vkernel_python_sliding_square_schemes SCRIPT sliding_square_schemes.py)
    new_python_test(NAME vkernel_python_add_objects SCRIPT test_add_objects.py PYTEST)
    new_python_test(NAME vkernel_python_nonos SCRIPT test_nonos.py PYTEST)
    new_python_test(NAME vkernel_python_remove_static_segment SCRIPT test_remove_static_segment.py
                    PYTEST)
    # fem.py and fem_tensor_analytical.py
    new_python_test(NAME vkernel_python_fem_tensors SCRIPT fem_tensors_test.py PYTEST)
    set(mesh_n5 ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/tests/mesh_N5.msh)
    new_python_test(NAME vkernel_python_two_fem_bodies SCRIPT two_fem_bodies.py ARGS ${mesh_n5})
    new_python_test(NAME vkernel_python_two_fem_bodies_gosi SCRIPT two_fem_bodies.py
                    ARGS --gosi ${mesh_n5})
    new_python_test(NAME vkernel_python_two_fem_bodies_dem SCRIPT two_fem_bodies.py
                    ARGS --dem ${mesh_n5})
    # the scripts checked against a reference trajectory (ref_trajectory.py)
    new_python_test(NAME vkernel_python_one_disk SCRIPT one_disk.py ARGS vnative)
    new_python_test(NAME vkernel_python_fixed_disk SCRIPT fixed-disk.py ARGS vnative 2)
    new_python_test(NAME vkernel_python_bouncing_disk SCRIPT bouncing_disk.py ARGS vnative 2)
    new_python_test(NAME vkernel_python_2d_flow SCRIPT 2d_flow_vnative.py ARGS 5 3 1)
    new_python_test(NAME vkernel_python_convexhull2d SCRIPT convexhull2d_vnative.py ARGS 1)
    new_python_test(NAME vkernel_python_disks SCRIPT disks.py ARGS vnative 2)
    new_python_test(NAME vkernel_python_2d_sphere_recirculation SCRIPT 2d_sphere_recirculation.py
                    ARGS vnative 0.5 0.5 20 3 0.5)
    new_python_test(NAME vkernel_python_bench_fem_disks SCRIPT bench_fem_disks.py
                    ARGS ${mesh_n5} 50 bench_fem_disks.hdf5 2.5 2.05 0.1 mj)
    new_python_test(NAME vkernel_python_bench_fem_disks_gosi SCRIPT bench_fem_disks.py
                    ARGS ${mesh_n5} 50 bench_fem_disks.hdf5 2.5 2.05 0.1 gosi)
    # the mesh of the FEM tests, generated again
    new_python_test(NAME vkernel_python_square_mesh SCRIPT square_mesh.py ARGS 5 mesh_N5.msh
                    ${mesh_n5})
  endif()
endif()
