import os,sys

import numpy
import math
import pickle

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

box_height = 3.683
box_length = 6.900
box_width  = 3.430

plan_thickness = 0.1

body_collection={}

body_collection['plan_id']= {}
id_plan=0


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




#create some bodies

# Creation of the hdf5 file for input/output
with Hdf5(use_compression=True) as io:

  ######### left_up
  id_plan=id_plan+1
  body_collection['plan_id']["left_up"]=id_plan
  v1 = numpy.array([0, 0 , box_height-1.200])
  v2 = numpy.array([4.370-4.370*1.200/box_height,0.0, 1.200])
  v3 = numpy.array([ 6.900, 0,1.200])

  left_up_normal = normal_plane(v1,v2,v3)
  print('left_up_normal=', left_up_normal)
  v1_extruded = v1 + numpy.dot(plan_thickness,left_up_normal)
  v2_extruded = v2 + numpy.dot(plan_thickness,left_up_normal)
  v3_extruded = v3 + numpy.dot(plan_thickness,left_up_normal)

  left_up_vertices=numpy.array([v1,v2,v3,v1_extruded,v2_extruded,v3_extruded])
  print left_up_vertices

  v1 = numpy.array([0, 0 , box_height])
  v2 = numpy.array([4.370-4.370*1.200/box_height,0.0, 1.200])
  v3 = numpy.array([ box_length, 0,1.200])
  v1_extruded = v1 + numpy.dot(plan_thickness,left_up_normal)
  v2_extruded = v2 + numpy.dot(plan_thickness,left_up_normal)
  v3_extruded = v3 + numpy.dot(plan_thickness,left_up_normal)

  left_up_vertices=numpy.array([v1,v2,v3,v1_extruded,v2_extruded,v3_extruded])
  print left_up_vertices

  io.addConvexShape('Left_up',left_up_vertices )
  io.addObject('left_up', [Contactor('Left_up')],
               translation=[0, 0, 0])


  ######### left_middle
  id_plan=id_plan+1
  body_collection['plan_id']["left_middle"]=id_plan

  v4 = numpy.array([4.370,1.280, 0.0])
  v5 = numpy.array([ 6.900-1.770, 1.280,0.0])

  left_middle_normal = normal_plane(v2,v4,v3)
  print('left_middle_normal=', left_middle_normal)

  v4_extruded = v4 + numpy.dot(plan_thickness,left_middle_normal)
  v5_extruded = v5 + numpy.dot(plan_thickness,left_middle_normal)

  left_middle_vertices=numpy.array([v2,v3,v4,v5,v2_extruded,v3_extruded,v4_extruded,v5_extruded])
  print left_middle_vertices

  io.addConvexShape('Left_middle',left_middle_vertices )
  io.addObject('left_middle', [Contactor('Left_middle')],
               translation=[0, 0, 0])

  ######### left_down
  id_plan=id_plan+1
  body_collection['plan_id']["left_down"]=id_plan

  v6 = numpy.array([4.370,box_width, -.6])
  v7 = numpy.array([6.900-1.770, box_width,-.6])

  left_down_normal = normal_plane(v4,v6,v5)
  print('left_down_normal=', left_down_normal)

  v6_extruded = v6 - [plan_thickness, 0.0 ,0.] + numpy.dot(plan_thickness,left_down_normal)
  v7_extruded = v7 + [plan_thickness, 0.0 ,0.] + numpy.dot(plan_thickness,left_down_normal)

  left_down_vertices=numpy.array([v4-[plan_thickness, 0.0 ,0.],v5+[plan_thickness, 0.0 ,0.],
                                  v6-[plan_thickness, 0.0 ,0.],v7+[plan_thickness, 0.0 ,0.],
                                  v4_extruded-[plan_thickness, 0.0 ,0.],v5_extruded+[plan_thickness, 0.0 ,0.],
                                  v6_extruded,v7_extruded])
  print left_down_vertices

  io.addConvexShape('Left_down',left_down_vertices )
  io.addObject('left_udown', [Contactor('Left_down')],
               translation=[0, 0, 0])

  ######### right_up
  id_plan=id_plan+1
  body_collection['plan_id']["right_up"]=id_plan

  v8 = numpy.array([0, box_width , box_height])
  v9 = numpy.array([ box_length, box_width,1.200])

  v10 = numpy.array([ 6.900-1.770, box_width,0.0])
  v11 =  numpy.array([4.370,box_width, 0.0])


  right_up_normal = normal_plane(v8,v9,v10)
  print('right_up_normal=', right_up_normal)

  v8_extruded = v8 + numpy.dot(plan_thickness,right_up_normal)
  v9_extruded = v9 + numpy.dot(plan_thickness,right_up_normal)
  v10_extruded = v10 + numpy.dot(plan_thickness,right_up_normal)
  v11_extruded = v11 + numpy.dot(plan_thickness,right_up_normal)

  right_up_vertices=numpy.array([v8,v9,v10,v11,v8_extruded,v9_extruded,v10_extruded,v11_extruded])
  print right_up_vertices

  io.addConvexShape('Right_up',right_up_vertices )
  io.addObject('right_up', [Contactor('Right_up')],
               translation=[0, 0, 0])

  ######### rear_up
  id_plan=id_plan+1
  body_collection['plan_id']["rear_up"]=id_plan


  rear_up_normal = normal_plane(v1,v8,v4)
  print('rear_up_normal=', rear_up_normal)

  v1_extruded = v1 + numpy.dot(plan_thickness,rear_up_normal)
  v2_extruded = v2 + numpy.dot(plan_thickness,rear_up_normal)
  v8_extruded = v8 + numpy.dot(plan_thickness,rear_up_normal)
  v4_extruded = v4 + numpy.dot(plan_thickness,rear_up_normal)
  v11_extruded = v11 + numpy.dot(plan_thickness,rear_up_normal)

  rear_up_vertices=numpy.array([v1-[0.0,plan_thickness,0.0],v2-[0.0,plan_thickness,0.0],
                                v8+[0.0,plan_thickness,0.0],v4-[0.0,plan_thickness,0.0],v11+[0.0,plan_thickness,0.0],
                                v1_extruded-[0.0,plan_thickness,0.0],v2_extruded-[0.0,plan_thickness,0.0],
                                v8_extruded+[0.0,plan_thickness,0.0],v4_extruded-[0.0,plan_thickness,0.0],
                                v11_extruded+[0.0,plan_thickness,0.0]])
  print rear_up_vertices

  io.addConvexShape('Rear_up',rear_up_vertices )
  io.addObject('rear_up', [Contactor('Rear_up')],
               translation=[0, 0, 0])


  ######### rear_down
  id_plan=id_plan+1
  body_collection['plan_id']["rear_down"]=id_plan
  #v12 = numpy.array([ 6.900-1.770, box_width,-.6])
  #v13 = numpy.array([4.370,box_width, -.6])


  rear_down_normal = normal_plane(v4,v11,v6)
  print('rear_down_normal=', rear_down_normal)

  v4_extruded = v4 + numpy.dot(plan_thickness,rear_down_normal)
  v11_extruded = v11 + numpy.dot(plan_thickness,rear_down_normal)
  v6_extruded = v6 + numpy.dot(plan_thickness,rear_down_normal)

  rear_down_vertices=numpy.array([v4,v11,v6,v4_extruded,v11_extruded,v6_extruded])
  print rear_down_vertices


  io.addConvexShape('Rear_down',rear_down_vertices )
  io.addObject('rear_down', [Contactor('Rear_down')],
               translation=[0, 0, 0])

  ######### front_up
  id_plan=id_plan+1
  body_collection['plan_id']["front_up"]=id_plan


  front_up_normal = normal_plane(v3,v5,v9)
  print('front_up_normal=', front_up_normal)

  v3_extruded = v3 + numpy.dot(plan_thickness,front_up_normal)
  v5_extruded = v5 + numpy.dot(plan_thickness,front_up_normal)
  v9_extruded = v9 + numpy.dot(plan_thickness,front_up_normal)
  v10_extruded = v10 + numpy.dot(plan_thickness,front_up_normal)

  front_up_vertices=numpy.array([v3-[0.0,plan_thickness,0.0],v5-[0.0,plan_thickness,0.0],
                                 v9+[0.0,plan_thickness,0.0],v10+[0.0,plan_thickness,0.0],
                                 v3_extruded-[0.0,plan_thickness,0.0],v5_extruded-[0.0,plan_thickness,0.0],
                                 v9_extruded+[0.0,plan_thickness,0.0],v10_extruded+[0.0,plan_thickness,0.0]])
  print front_up_vertices


  io.addConvexShape('Front_up',front_up_vertices )
  io.addObject('front_up', [Contactor('Front_up')],
               translation=[0, 0, 0])

  ######### front_down
  id_plan=id_plan+1
  body_collection['plan_id']["front_down"]=id_plan


  front_down_normal = normal_plane(v5,v7,v10)
  print('front_down_normal=', front_down_normal)

  v7_extruded = v7 + numpy.dot(plan_thickness,front_down_normal)
  v5_extruded = v5 + numpy.dot(plan_thickness,front_down_normal)
  v10_extruded = v10 + numpy.dot(plan_thickness,front_down_normal)

  front_down_vertices=numpy.array([v5,v7,v10,v5_extruded,v7_extruded,v10_extruded])
  print front_down_vertices

  io.addConvexShape('Front_down',front_down_vertices )
  io.addObject('front_down', [Contactor('Front_down')],
               translation=[0, 0, 0])


  n_cube=5
  n_row=1
  n_col=1
  cube_size =0.25
  x_shift=3.0
  for i in range(n_row):
    for j in range(n_col):
      for n in range(n_cube):
        # Definition of a cube as a convex shape
        io.addConvexShape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j), [ (-cube_size, cube_size, -cube_size),
                                                                   (-cube_size, -cube_size, -cube_size),
                                                                   (-cube_size, -cube_size, cube_size),
                                                                   (-cube_size, cube_size, cube_size),
                                                                   (cube_size, cube_size, cube_size),
                                                                   (cube_size, cube_size, -cube_size),
                                                                   (cube_size, -cube_size, -cube_size),
                                                                   (cube_size, -cube_size, cube_size)])


        io.addObject('cube'+str(n)+'_'+str(i)+'_'+str(j), [Contactor('CubeCS'+str(n)+'_'+str(i)+'_'+str(j))],
                     translation=[2+i*x_shift*cube_size, 2+x_shift*j*cube_size, 3+cube_size*x_shift*n],
                     velocity=[0, 0, 0, 0, 0, 0],
                     mass=1)


  # Definition of a non smooth law. As no group ids are specified it
  # is between contactors of group id 0.
  io.addNewtonImpactFrictionNSL('contact', mu=0.3)

  print body_collection
  f = open('body_collection.dict', 'w')
  pickle.dump(body_collection,f)
  f.close()


step=1000
hstep=0.005

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

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
         gravity_scale=0.1,
         t0=0,
         T=step*hstep,
         h=hstep,
         multipoints_iterations=True,
         theta=0.50001,
         Newton_max_iter=1,
         set_external_forces=None,
         solver=Numerics.SICONOS_FRICTION_3D_NSGS,
         itermax=100,
         tolerance=1e-4,
         numerics_verbose=False,
         output_frequency=10)
