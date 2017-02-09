import os,sys

import numpy
import math
import pickle

import random

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
#sys.path.append('../..')
#from mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel
# WARNING : in 3D by default z-axis is upward
# this is very important to direct PLANx objects

dim = 3

unscaled_bar_length=1.5
aspect_ratio=100.0
unscaled_bar_height=unscaled_bar_length/aspect_ratio
unscaled_bar_width=unscaled_bar_length/aspect_ratio

unscaled_volume = unscaled_bar_length*unscaled_bar_height*unscaled_bar_width
unscaled_density=1000

unscaled_mass=unscaled_volume*unscaled_density

print('unscaled_mass',unscaled_mass)


scale=1.0/unscaled_bar_length*1.0

density =  unscaled_density/(scale**3)

bar_height =  unscaled_bar_height*scale
bar_length =  unscaled_bar_length*scale
bar_width  =  unscaled_bar_width*scale

body_collection={}
body_collection['plan_id']= {}
id_plan=0
# scale =1
# mass :3.375000e-01
# Inertia :
# 3.600000e-04, 0.000000e+00, 0.000000e+00,
# 0.000000e+00, 1.195312e-01, 0.000000e+00,
# 0.000000e+00, 0.000000e+00, 1.195312e-01,

#create some bodies
# Creation of the hdf5 file for input/output
with Hdf5() as io:
  volume = bar_height * bar_length * bar_width
  mass = volume*density
  print('mass', mass)
  print('scale', scale)
  # raw_input()
  # Definition of a cube as a convex shape
  io.addConvexShape('Bar', [ (-bar_length,  bar_width, -bar_height),
                             (-bar_length, -bar_width, -bar_height),
                             (-bar_length, -bar_width,  bar_height),
                             (-bar_length,  bar_width,  bar_height),
                             ( bar_length , bar_width,  bar_height),
                             ( bar_length,  bar_width, -bar_height),
                             ( bar_length, -bar_width, -bar_height),
                             ( bar_length ,-bar_width,  bar_height)])

  angle= math.pi/8.0
  trans=[0,0,4.0*scale]
  ori = [math.cos(angle/2.0),0.0,math.sin(angle/2.0),0]
  axis = numpy.zeros(3)
  angle_test = Kernel.axisAngleFromQuaternion(trans+ori, axis)
  print angle_test,axis
  print('ori initial', ori)
  io.addObject('bar', [Contactor('Bar')],
               translation=trans,
               orientation = ori,
               velocity=[0, 0, 0, 0, 0.0, 0],
               mass=mass)

  # Definition of the ground shape
  io.addPrimitiveShape('Ground', 'Box', (5*scale, 5*scale, 0.1*scale))
  angleground=  -math.pi/4.0
  axis = [1.0, 0.0, 0.0]
  origround = [math.cos(angleground/2.0),
               axis[0] * math.sin(angleground/2.0),
               axis[1] * math.sin(angleground/2.0),
               axis[2] * math.sin(angleground/2.0)]
  io.addObject('ground', [Contactor('Ground')],
               translation=[0, 0, 0.0],
               orientation = origround)


  # Definition of a non smooth law. As no group ids are specified it
  # is between contactors of group id 0.
  io.addNewtonImpactFrictionNSL('contact', mu=0.3)

  print body_collection

  f = open('body_collection.dict', 'w')
  pickle.dump(body_collection,f)
  f.close()


step=2000
hstep=0.001

gravity_scale=1.0/scale
import scipy.constants as constants

def apply_forces(body):
  g = constants.g / gravity_scale
  weight = [0, 0, - body.scalarMass() * g]
  body.setFExtPtr(weight)



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
         gravity_scale=gravity_scale,
         t0=0,
         T=step*hstep,
         h=hstep,
         multipoints_iterations=True,
         theta=1.0,
         Newton_max_iter=10,
         set_external_forces=apply_forces,
         solver=Numerics.SICONOS_FRICTION_3D_NSGS,
         itermax=1000,
         tolerance=1e-12,
         numerics_verbose=False,
         violation_verbose=True,
         output_frequency=10)
