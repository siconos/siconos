import os,sys
import numpy
import math
import pickle
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

from siconos.mechanics.collision.convexhull import ConvexHull

from siconos.mechanics.collision.bullet import SiconosBulletOptions
options = SiconosBulletOptions()
options.worldScale = 1.0
options.contactBreakingThreshold = 0.01

plan_thickness = 0.05

density = 2679.1838



def normal_plane(p1,p2,p3):

  x1=p1[0]
  y1=p1[1]
  z1=p1[2]
  x2=p2[0]
  y2=p2[1]
  z2=p2[2]
  x3=p3[0]
  y3=p3[1]
  z3=p3[2]

  vector1 = [x2 - x1, y2 - y1, z2 - z1]
  vector2 = [x3 - x1, y3 - y1, z3 - z1]

  cross_product = [vector1[1] * vector2[2] - vector1[2] * vector2[1], -1 * (vector1[0] * vector2[2] - vector1[2] * vector2[0]), vector1[0] * vector2[1] - vector1[1] * vector2[0]]

  a = cross_product[0]
  b = cross_product[1]
  c = cross_product[2]
  d = - (cross_product[0] * x1 + cross_product[1] * y1 + cross_product[2] * z1)

  return numpy.array([a,b,c])/numpy.linalg.norm([a,b,c])

#create plans

# Creation of the hdf5 file
with Hdf5(use_compression=True) as io:

  ######### amont
  v0 = numpy.array([5.00,1.6131,1.0751])
  v1 = numpy.array([-2.50, 1.6131, 1.0751])
  v2 = numpy.array([-2.50,0.9354, 0.3535])
  v3 = numpy.array([ 5.00, 0.9354, 0.3535])

  amont_normal = normal_plane(v1,v2,v3)
  print('amont_normal=', amont_normal)

  v0_extruded = v0 + numpy.dot(plan_thickness,amont_normal)
  v1_extruded = v1 + numpy.dot(plan_thickness,amont_normal)
  v2_extruded = v2 + numpy.dot(plan_thickness,amont_normal)
  v3_extruded = v3 + numpy.dot(plan_thickness,amont_normal)

  amont_vertices=numpy.array([v0,v1,v2,v3, v0_extruded,v1_extruded,v2_extruded,v3_extruded])
  print amont_vertices

  io.addConvexShape('amont',amont_vertices )
  io.addObject('amont', [Contactor('amont')],
               translation=[1.50, -1.45, -1.5331])

  ######### aval
  v4 = numpy.array([-2.50,0,0])
  v5 = numpy.array([ 5.00,0,0])

  aval_normal = normal_plane(v2,v4,v3)
  print('aval_normal=', aval_normal)

  v4_extruded = v4 + numpy.dot(plan_thickness,aval_normal)
  v5_extruded = v5 + numpy.dot(plan_thickness,aval_normal)

  aval_vertices=numpy.array([v2,v3,v4,v5,v2_extruded,v3_extruded,v4_extruded,v5_extruded])
  print aval_vertices

  io.addConvexShape('aval',aval_vertices )
  io.addObject('aval', [Contactor('aval')],
               translation=[1.50, -1.45, -1.5331])

  ######### sol
  v6 = numpy.array([5.00,-5.00,0])
  v7 = numpy.array([-2.50,-5.00,0])

  sol_normal = normal_plane(v4,v6,v5)
  print('sol_normal=', sol_normal)

  v6_extruded = v6 - [plan_thickness, 0.0 ,0.] + numpy.dot(plan_thickness,sol_normal)
  v7_extruded = v7 + [plan_thickness, 0.0 ,0.] + numpy.dot(plan_thickness,sol_normal)

  sol_vertices=numpy.array([v4-[plan_thickness, 0.0 ,0.],v5+[plan_thickness, 0.0 ,0.],
                                  v6-[plan_thickness, 0.0 ,0.],v7+[plan_thickness, 0.0 ,0.],
                                  v4_extruded-[plan_thickness, 0.0 ,0.],v5_extruded+[plan_thickness, 0.0 ,0.],
                                  v6_extruded,v7_extruded])
  print sol_vertices

  io.addConvexShape('sol',sol_vertices )
  io.addObject('sol', [Contactor('sol')],
               translation=[1.50, -1.45, -1.5331])


  n_cube=1
  n_row=50
  n_col=1
  cube_size =0.0144
  x_shift=0.030
  x_translate = 0.1
  sphere_count =0
  spheres = []
  radius = 0.005
  for i in range(n_row):
    for j in range(n_col):
      for n in range(n_cube):
        # Definition of a cube
        vertices = [ (-cube_size, cube_size, -cube_size),
                     (-cube_size, -cube_size, -cube_size),
                     (-cube_size, -cube_size, cube_size),
                     (-cube_size, cube_size, cube_size),
                     (cube_size, cube_size, cube_size),
                     (cube_size, cube_size, -cube_size),
                     (cube_size, -cube_size, -cube_size),
                     (cube_size, -cube_size, cube_size)]
        io.addConvexShape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j), vertices)

        for v in vertices:
          sphere_count += 1
          spheres.append('Sphere%03d'%sphere_count)
          io.addPrimitiveShape(spheres[-1], 'Sphere', [radius])

        # computation of inertia and volume
        ch = ConvexHull(vertices)
        inertia,volume=ch.inertia(ch.centroid())
        #print inertia, volume
        #raw_input()
        #angle_init = math.pi/4.0
        #angle_init = 0.001
        angle_init = 0.0
        #contactor = [Shape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j))]
        contactor = []
        for sph,loc in zip(spheres, vertices):
          contactor.append(Contactor(sph, relative_translation=loc))

        io.addObject('cube'+str(n)+'_'+str(i)+'_'+str(j),
                     contactor,
                     translation=[i*(x_translate+x_shift*cube_size), x_shift*j*(x_translate+cube_size), (x_translate+cube_size*x_shift)*n],
                     velocity=[0, 0, 0, 0, 0, 0],
                     orientation=[math.cos(angle_init/2.0),0,math.sin(angle_init/2.0),0],
                     mass=volume*density,
                     inertia=inertia*density)


  # Definition of a non smooth law
  io.addNewtonImpactFrictionNSL('contact', e=0.01, mu=0.9)


step=10000
hstep=0.0005

with Hdf5(mode='r+') as io:


  io.run(with_timer=False,
         time_stepping=None,
         space_filter=None,
         body_class=None,
         shape_class=None,
         face_class=None,
         edge_class=None,
         t0=0,
         T=step*hstep,
         h=hstep,
         multipoints_iterations=True,
         theta=0.50001,
         Newton_max_iter=1,
         set_external_forces=None,
         solver=Numerics.SICONOS_FRICTION_3D_NSGS,
         itermax=100,
         tolerance=1e-8,
         numerics_verbose=False,
         output_frequency=10)
