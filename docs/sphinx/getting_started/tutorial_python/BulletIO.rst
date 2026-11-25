.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _BulletIO_example:

Collections of rigid bodies with Bullet-based contact detection (Siconos/Mechanics)
-----------------------------------------------------------------------------------------

Authors: Maurice Bremond, Vincent Acary 2013--2016

The source code of this set of examples is `available here <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/-/tree/main/examples/mechanics/GeometricPrimitives?ref_type=heads>`_.


Description of the physical problems: collections of rigid bodies with contact and Coulomb friction
```````````````````````````````````````````````````````````````````````````````````````````````````

In this set of examples (cubes.py, n_cubes.py, ...), we model collections of rigid bodies
associated with shapes (primitives (sphere, cube, capsule, etc.), convex hulls, meshes) that interact
through contact and friction.

The examples make extensive use of the mechanics_run Python code, which creates and manages the geometrical
data and the results of the simulation in an HDF5 file.

The results can be viewed with the siconos_vview and siconos_vexport Python scripts.

Brief description of the examples:

* cube.py is the simplest example, with a few cubes falling onto a fixed rigid plane. The scene is built
  with calls to mechanics_run methods such as

  * addPrimitiveShape, to create a primitive shape


  .. code-block:: python

      with Hdf5(mode='r+') as io:
	io.addPrimitiveShape('Ground', 'Box', (100, 100, .5))

  * addConvexShape, to create a convex shape (see also addPrimitiveShape, addNewtonImpactFrictionNSL)

  .. code-block:: python

     with Hdf5(mode='r+') as io:
	io.addConvexShape('Cube_1', [
	(-1.0, 1.0, -1.0),
	(-1.0, -1.0, -1.0),
	(-1.0, -1.0, 1.0),
	(-1.0, 1.0, 1.0),
	(1.0, 1.0, 1.0),
	(1.0, 1.0, -1.0),
	(1.0, -1.0, -1.0),
	(1.0, -1.0, 1.0)])

  * addObject, to associate a Newton-Euler dynamical system with a shape

  .. code-block:: python

      with Hdf5(mode='r+') as io:
        io.addObject('cube', [Contactor('Cube')], translation=[0, 0, 2],
	velocity=[10, 0, 0, 1, 1, 1],
	mass=1)

  The computation is launched with the run() method, here with its default arguments:

  .. code-block:: python

      with Hdf5(mode='r+') as io:
        io.run(with_timer=False,
            time_stepping=None,
            space_filter=None,
            body_class=None,
            shape_class=None,
            face_class=None,
            edge_class=None,
            gravity_scale=1,
            t0=0,
            T=10,
            h=0.0005,
            multipoints_iterations=True,
            theta=0.50001,
            Newton_max_iter=20,
            set_external_forces=None,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=100000,
            tolerance=1e-8,
            numerics_verbose=False,
            output_frequency=None)

* n_cubes.py: an extension of cubes.py that builds a rectangular pile of cubes.

* cube_scene.py, cube_simulation.py, n_cubes_scene.py, n_cubes_simulation.py: the same examples,
  but the construction of the scene and the simulation are in separate files. This way, one scene file can be produced
  and then used for several simulations. Furthermore, if the result of a simulation is stored in the HDF5 file, the next simulation restarts from the last
  state of the system. Several simulations can thus be chained in time.

* MultipleContactorsAndNSlaws.py: two shapes are attached to the same mechanical system. This makes it possible to build clusters and to get
  complex contact shapes.

* bar.py, bar_contact.py: simple simulations where the performance of the Lie integrator matters for the long-term, stable
  simulation of rotating bodies. They also test large inertia ratios in the solver.

* convexhull.py: shows how to build a collection of convex hull shapes.

* PairWise_test_bullet.py: tests contacts between pairs of shapes of various kinds (primitives, convex hulls). Meshes still have to be added.

* chute.py: a more complex example where a collection of polyhedra falls into a hopper.
