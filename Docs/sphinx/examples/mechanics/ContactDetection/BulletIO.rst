.. _BulletIO_example:

Collections of rigid bodies with Bullet based contact detection (Siconos/Mechanics)
-----------------------------------------------------------------------------------------

Author: Maurice Bremond, Vincent Acary 2013--2016

You may refer to the source code of this set of  examples, `found here <https://github.com/siconos/siconos/tree/master/examples/Mechanics/MechanicsIO/>`_.


Description of the physical problems : rigid bodies collection with contact and Coulomb friction
````````````````````````````````````````````````````````````````````````````````````````````````

In this set of examples (cubes.py, n_cubes.py, ...), we model collections of rigid bodies
associated with shapes (primitive (sphere, cube capsule, etc), convex hull, mesh) that interact
with contact and friction (see :ref:`fig-n_cubes`).

The example use extensively the mechanics_run python code that enable to create and manage the geometrical
data and the results of the simulation in hdf5 file.

The results cab be viewed with the siconos_vview and the siconos_vexport python scripts

A brief description of examples is as follows :

* cube.py is  the simplest example with few cubes that fall down on a rigid fixed plane. The scene is built
  using calls to mechanics_run methods such as
  
  * addPrimitiveShape for creating a primitive shape ::
    
      with Hdf5(mode='r+') as io:
	io.addPrimitiveShape('Ground', 'Box', (100, 100, .5))
    
  * addConvexShape for creating a convex shape addPrimitiveShape, addNewtonImpactFrictionNSL) ::
      
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
      
  * addObject for associating Newton Euler Dynamical System to a shape::

      with Hdf5(mode='r+') as io:
        io.addObject('cube', [Contactor('Cube')], translation=[0, 0, 2],
	velocity=[10, 0, 0, 1, 1, 1],
	mass=1)

  The computation is launched using the method run() with default arguments::

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
	    
* n_cubes.py . This is an extension of cubes.py where it is possible to build a reactangular pile of cubes
  
* cube_scene.py, cube_simulation.py n_cubes_scene.py n_cubes_simulation.py . These examples are a different treatment of the same example
  but the construction of the scene and the simulation are separated into differents. In that way, it is possible to produce a scene file
  and then to perfoms several simulations. Furthermore, if the result of the simulation is stored in the hdf5 file, we restart from the last
  state of the system to continue the simulation. it is alo possible to concatenate in time several simulations in that way.

* MultipleContactorsAndNSlaws.py . This  is an example where two shapes are fixed to the same mechanical system. It allows to build clusters and to get
  some complex contact shapes.

* bar.py, bar_contact.py  This simple example provides one with a simple simulation where the performance of Lie integrator are important for long term
  simulation of rotating bodie in a stable way. It allows also to test large inertia ratio in the solver.

* convexhull.py . This test shos how to build a collection of convexhull shapes

* PairWise_test_bullet.py . This test file performs the test of contacting pairs of various nature (primitive, convex hull). Meshes have to be added.

* chute.py. This example is a little more complex one where the a collection of polyhedra falls down into an hopper.
  
