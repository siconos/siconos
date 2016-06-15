import os,sys

import numpy
import math
import pickle

import random
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull

from siconos.io.mechanics_io import Hdf5
#sys.path.append('../..')
#from mechanics_io import Hdf5
import siconos.numerics as Numerics

unscaled_polyhedron_size=0.5
unscaled_density=2500

scale=1.0/unscaled_polyhedron_size*1000.0

polyhedron_size = unscaled_polyhedron_size*scale
density =  unscaled_density/(scale**3)

box_height = 3.683*scale
box_width  = 3.430*scale

#create some bodies

# Creation of the hdf5 file for input/output
with Hdf5() as io:

  angle= 0.0 ; math.pi/4.0
  translation = [0.0, 0.0 ,0.0]
  orientation = [math.cos(angle/2.0), 0.0, math.sin(angle/2.0), 0.0]
  orientation = [1.0, 0.0, 0.0, 0.0]
  #orientation = [0.0, 1.0, 0.0, 0.0] #--> cos(\theta/2.0) = 0, sin(\theta/2.0) =1 ==> \theta = pi

  ######### left_up
  v1 = numpy.array([0, 0 , 1.0*box_height])
  v2 = numpy.array([box_width,box_width,0.0])
  v3 = numpy.array([box_width,-box_width,0.0])
  v4 = numpy.array([-box_width,-box_width,0.0])
  v5 = numpy.array([-box_width,box_width,0.0])

  left_up_vertices=numpy.array([v1,v2,v3,v4, v5])
  print left_up_vertices
  io.addConvexShape('Left_up',left_up_vertices )
  io.addObject('left_up', [Contactor('Left_up')],
               translation=translation)

  ######### right_up
  io.addConvexShape('Right_up',left_up_vertices )
  translation = [0.0, 2.0*box_width ,0.0]
  io.addObject('right_up', [Contactor('Right_up')],
               translation=translation)

  n_polyhedron=1
  n_row=1
  n_col=1
  x_shift=3.0

  angle =math.pi/2.0
  orientation_polyhedron = [math.cos(angle/2.0), math.sin(angle/2.0), 0.0, 0.0]
  orientation_polyhedron = [math.cos(angle/2.0), 0.0, math.sin(angle/2.0), 0.0]
  #orientation_polyhedron = [1.0, 0.0, 0.0, 0.0]
  #orientation_polyhedron = [0.0, 1.0, 0.0, 0.0] #--> cos(\theta/2.0) = 0, sin(\theta/2.0) =1 ==> \theta = pi


  for i in range(n_row):
    for j in range(n_col):
      for n in range(n_polyhedron):
        polyhedron_size_rand = polyhedron_size*(1.0 + 0.5*( random.random()-1.0))
        polyhedron_vertices=[ (-polyhedron_size_rand, polyhedron_size_rand, 0.0),
                        (-polyhedron_size_rand, -polyhedron_size_rand, 0.0),
                        (polyhedron_size_rand, -polyhedron_size_rand, 0.0),
                        (polyhedron_size_rand, polyhedron_size_rand, 0.0),
                        (0.0,0.0, 10.0*polyhedron_size_rand)]
        print('polyhedron_vertices', polyhedron_vertices)
        ch = ConvexHull(polyhedron_vertices)
        cm = ch.centroid()
        print('cm', cm)

        # correction of vertices such that o is the centroid
        polyhedron_vertices = numpy.array(polyhedron_vertices)[:]-cm[:]
        print('corrected polyhedron_vertices', polyhedron_vertices)
        ch = ConvexHull(polyhedron_vertices)
        cm = ch.centroid()
        print('cm', cm)

        # computation of inertia and volume
        inertia,volume=ch.inertia(cm)
        print ('inertia,volume',inertia,volume)
        mass = volume*density
        inertia = inertia*density
        print ('inertia,mass',inertia,mass)

        # Definition of a polyhedron as a convex shape
        io.addConvexShape('PolyhedronCS'+str(n)+'_'+str(i)+'_'+str(j), polyhedron_vertices)

        trans=[box_width/5.0+i*x_shift*polyhedron_size, x_shift*(j+2)*polyhedron_size, box_height+polyhedron_size*x_shift*n]
        #trans= [0,0,0]
        io.addObject('polyhedron'+str(n)+'_'+str(i)+'_'+str(j), [Contactor('PolyhedronCS'+str(n)+'_'+str(i)+'_'+str(j))],
                     translation=trans, orientation = orientation_polyhedron,
                     velocity=[0, 0, 0, 0, 0, 0],
                     mass=mass,inertia=inertia)

  # Definition of a non smooth law. As no group ids are specified it
  # is between contactors of group id 0.
  io.addNewtonImpactFrictionNSL('contact', mu=0.3)

step=10000
hstep=1e-4

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+', collision_margin=0.01) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
  io.run(with_timer=False,
         time_stepping=None,
         space_filter=None,
         body_class=None,
         shape_class=None,
         face_class=None,
         edge_class=None,
         gravity_scale=1.0/scale,
         t0=0,
         T=step*hstep,
         h=hstep,
         multipoints_iterations=True,
         theta=0.50001,
         Newton_max_iter=10,
         set_external_forces=None,
         solver=Numerics.SICONOS_FRICTION_3D_NSGS,
         itermax=1000,
         tolerance=1e-4,
         numerics_verbose=False,
         violation_verbose=True,
         output_frequency=100)
