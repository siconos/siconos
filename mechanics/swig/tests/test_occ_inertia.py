#!/usr/bin/env python
#

from siconos.mechanics.collision.tools import Contactor, Shape, Volume, Material
from siconos.io.mechanics_run import MechanicsHdf5Runner
from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere
from OCC.gp import gp_Pnt, gp_Ax2, gp_Dir
#import mechanics_run
#mechanics_run.set_backend('occ')
import siconos.io.mechanics_run

siconos.io.mechanics_run.set_backend('occ')
import numpy as np
import math

# create sphere
radius = 1.0
# The sphere center
X1 = 0.0
Y1 = 0.0
Z1 = 0.0
# create OCC.gp.gp_Pnt-Point from vector
point = gp_Pnt( X1, Y1, Z1 )
sphere_r1 = BRepPrimAPI_MakeSphere( point, radius )
sphere_r1_shape = sphere_r1.Shape()

radius = 0.01
sphere_r01 = BRepPrimAPI_MakeSphere( point, radius )
sphere_r01_shape = sphere_r01.Shape()

radius = 0.001
sphere_r001 = BRepPrimAPI_MakeSphere( point, radius )
sphere_r001_shape = sphere_r001.Shape()






from OCC.BRep import BRep_Builder
from OCC.TopoDS import TopoDS_Compound

builder = BRep_Builder()
comp = TopoDS_Compound()
builder.MakeCompound(comp)

from OCC.BRep import BRep_Builder
from OCC.TopoDS import TopoDS_Compound

builder.Add(comp, sphere_r1_shape)

radius=1.0
point = gp_Pnt( 4., 0., 0. )
sphere_r1_t = BRepPrimAPI_MakeSphere( point, radius )
sphere_r1_t_shape = sphere_r1_t.Shape()

builder.Add(comp, sphere_r1_t_shape)


# this cylinder is defined for the visualisation of the axis edge
# it may be used for contact also.
axis = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, -.05, 0),gp_Dir(0, 1
                                            , 0)), .001, .1).Shape()

density=7750.0
steel = Material(density=density)
water = Material(density=1000)

steel_cm = Material(density=density*1e-6)

tol  =1e-10



def test_sphere_1m():
    with MechanicsHdf5Runner() as io:
        io.add_occ_shape('sphere_r1_shp', sphere_r1_shape)

        ###############
        # Sphere 1m
        ###############

        obj= io.add_object('sphere_bdy',
                      [Volume(shape_name='sphere_r1_shp',
                              instance_name='sphere_r1_shp_bdy',
                              parameters=steel)],
                       translation=[0., 0., 0.],
                       velocity=[0., 0., 0., 0., 0., 0.])

        radius =1.0
        volume = 4/3. * math.pi * radius**3
        mass= volume * steel.density

        inertia = 2/5*mass*radius**2
        print('volume :', volume)
        print('mass :', mass)
        print('inertia :', inertia)

        # print(obj.attrs['id'])
        c_mass=obj.attrs['mass']
        c_inertia = obj.attrs['inertia'][0][0]
        assert(abs(mass-c_mass) < tol)
        assert(abs(inertia-c_inertia) < tol)


def test_two_sphere_1m():
    with MechanicsHdf5Runner() as io:

        io.add_occ_shape('two_spheres_r1_shp', comp)

        #################
        # Two spheres 1m
        #################

        obj= io.add_object('two_spheres_comp',
                      [Volume(shape_name='two_spheres_r1_shp',
                              instance_name='comp',
                              parameters=steel)],
                           translation=[0., 0., 0.],
                           velocity=[0., 0., 0., 0., 0., 0.])

        # obj= io.add_object('Two_Sphere_com',
        #               [Volume(shape_name='two_sphere_R=1_shp',
        #                       instance_name='balance_wheel_shp_bdy',
        #                       parameters=steel)],
        #                translation=[0., 0., 0.],
        #                velocity=[0., 0., 0., 0., 0., 0.])
        radius =1.0
        volume = 2* 4/3. * math.pi * radius**3
        mass= volume * steel.density
        print('volume :', volume)
        print('mass :', mass)
        d = 2
        inertia_zy = 2*(2/5*mass/2.0*radius**2 + mass/2.0*d**2)
        print('inertia (z and y) :', inertia_zy)
        inertia_x = (2/5*mass*radius**2)
        print('inertia (x) :', inertia_x)
        c_mass=obj.attrs['mass']
        c_inertia_zy = obj.attrs['inertia'][1][1]
        c_inertia_x = obj.attrs['inertia'][0][0]

        assert(abs(mass-c_mass) < tol)
        assert(abs(inertia_x-c_inertia_x) < tol)
        assert(abs(inertia_zy-c_inertia_zy) < tol)


def test_two_sphere_1m_steel_water():
    with MechanicsHdf5Runner() as io:
        io.add_occ_shape('sphere_r1_shp', sphere_r1_shape)

        #################
        # Two spheres 1m steel and water
        #################

        obj = io.add_object('two_spheres_bdy',
                      [Volume(shape_name='sphere_r1_shp',
                              instance_name='sphere_r1_shp_bdy',
                              parameters=steel,
                              relative_translation=[0., -2., 0.]),
                       Volume(shape_name='sphere_r1_shp',
                              instance_name='sphere_r1_shp_bdy',
                              parameters=water,
                              relative_translation=[0, +2., 0.])],
                       translation=[0., 0., 0.],
                       velocity=[0., 0., 0., 0., 0., 0.])

        radius =1.0
        volume =  4/3. * math.pi * radius**3
        mass1= volume * steel.density
        mass2= volume * water.density
        mass= mass1 + mass2

        print('volume :', volume)
        print('mass :', mass)
        inertia_y = 2/5*mass*radius**2
        print('inertia (y) :', inertia_y)

        d1=2-1.54285714286
        d2=2+1.54285714286
        inertia_zx = mass1*d1**2 + 2/5*mass1*radius**2 +  mass2*d2**2 + 2/5*mass2*radius**2
        print('inertia (z and x) :', inertia_zx)
        print(obj.attrs['id'])
        c_mass=obj.attrs['mass']
        c_inertia_zx = obj.attrs['inertia'][0][0]
        c_inertia_y = obj.attrs['inertia'][1][1]

        assert(abs(mass-c_mass) < 1e-08)
        assert(abs(inertia_y-c_inertia_y) < tol)
        assert(abs(inertia_zx-c_inertia_zx) < tol)


def test_sphere_01m():
    with MechanicsHdf5Runner() as io:

        io.add_occ_shape('sphere_r01_shp', sphere_r01_shape)

        ###############
        # Sphere 0.01 m
        ###############

        print('############### sphere of Radius 0.01 m')
        obj= io.add_object('Sphere01_com',
                           [Volume(shape_name='sphere_r01_shp',
                                   instance_name='sphere_r01_shp_bdy',
                                   parameters=steel)],
                           translation=[0., 0., 0.],
                       velocity=[0., 0., 0., 0., 0., 0.])

        radius =0.01
        volume = 4/3. * math.pi * radius**3
        mass= volume * steel.density

        inertia = 2/5*mass*radius**2
        print('volume :', volume)
        print('mass :', mass)
        print('inertia :', inertia)
        c_mass=obj.attrs['mass']
        c_inertia = obj.attrs['inertia'][0][0]
        assert(abs(mass-c_mass) < tol)
        assert(abs(inertia-c_inertia) < tol)

def test_sphere_001m():
    with MechanicsHdf5Runner() as io:

        io.add_occ_shape('sphere_r001_shp', sphere_r001_shape)

        ###############
        # Sphere 0.001 m
        ###############

        # to compute the center of mass and perform the right translation
        print('############### sphere of Radius 0.01 m')
        obj= io.add_object('Sphere001_com',
                           [Volume(shape_name='sphere_r001_shp',
                                   instance_name='sphere_r001_shp_bdy',
                                   parameters=steel)],
                           translation=[0., 0., 0.],
                       velocity=[0., 0., 0., 0., 0., 0.])

        radius =0.001
        volume = 4/3. * math.pi * radius**3
        mass= volume * steel.density

        inertia = 2/5*mass*radius**2
        print('volume :', volume)
        print('mass :', mass)
        print('inertia :', inertia)
        c_mass=obj.attrs['mass']
        c_inertia = obj.attrs['inertia'][0][0]
        assert(abs(mass-c_mass) < tol)
        assert(abs(inertia-c_inertia) < tol)


def test_sphere_1cm():
    with MechanicsHdf5Runner() as io:
        io.add_occ_shape('sphere_r1_shp', sphere_r1_shape)
        ###############
        # Sphere 1 cm
        ###############

        # to compute the center of mass and perform the right translation
        print('############### sphere of Radius 1cm')
        obj = io.add_object('sphere_cm_bdy',
                      [Volume(shape_name='sphere_r1_shp',
                              instance_name='sphere_r1_shp_bdy',
                              parameters=steel_cm)],
                       translation=[0., 0., 0.],
                       velocity=[0., 0., 0., 0., 0., 0.])

        radius =1.0
        volume = 4/3. * math.pi * radius**3
        mass= volume * steel_cm.density

        inertia = 2/5*mass*radius**2 # * 1e-4 # Warning the inertia unit is kg.m^2, so it is not sufficient to scale to density.

        print('volume :', volume)
        print('mass :', mass)
        print('inertia :', inertia)
        c_mass=obj.attrs['mass']
        c_inertia = obj.attrs['inertia'][0][0]
        assert(abs(mass-c_mass) < tol)
        assert(abs(inertia-c_inertia) < tol)

if __name__ == '__main__':
    test_sphere_1m()
    test_two_sphere_1m()
    test_two_sphere_1m_steel_water()
    test_sphere_01m()
    test_sphere_001m()
    test_sphere_1cm()
