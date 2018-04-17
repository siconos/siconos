#!/usr/bin/env python

__all__ = ['create_chute']

import os, sys

import numpy
import math
import pickle

from siconos.mechanics.collision.tools import Contactor
import siconos.numerics as Numerics

# WARNING : in 3D by default z-axis is upward
# this is very important to direct PLANx objects

dim = 3

box_height = 3.683
box_length = 6.900
box_width  = 3.430

plane_thickness = 0.1

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
def create_chute(io, box_height = box_height,
                 box_length = box_length,
                 box_width = box_width,
                 plane_thickness = plane_thickness,
                 scale = 1.0, trans = [0,0,0]):
    print('Creation of the hopper')
    box_height *= scale
    box_length *= scale
    box_width *= scale
    plane_thickness *= scale

    ######### left_up
    v1 = numpy.array([0, 0, box_height])
    v2 = numpy.array([4.370*scale-4.370*1.200*scale*scale/box_height, 0.0, 1.200*scale])
    v3 = numpy.array([box_length, 0, 1.200*scale])
    left_up_normal = normal_plane(v1,v2,v3)
    v1_extruded = v1 + numpy.dot(plane_thickness,left_up_normal)
    v2_extruded = v2 + numpy.dot(plane_thickness,left_up_normal)
    v3_extruded = v3 + numpy.dot(plane_thickness,left_up_normal)

    left_up_vertices=numpy.array([v1,v2,v3,v1_extruded,v2_extruded,v3_extruded])
    # print left_up_vertices

    io.add_convex_shape('Left_up',left_up_vertices )
    io.add_object('left_up', [Contactor('Left_up')],
                 translation = trans)


    ######### left_middle
    v4 = numpy.array([4.370*scale, 1.280*scale, 0.0])
    v5 = numpy.array([(6.900-1.770)*scale, 1.280*scale, 0.0])

    left_middle_normal = normal_plane(v2,v4,v3)
    # print('left_middle_normal=', left_middle_normal)

    v4_extruded = v4 + numpy.dot(plane_thickness, left_middle_normal)
    v5_extruded = v5 + numpy.dot(plane_thickness, left_middle_normal)

    left_middle_vertices=numpy.array([v2,v3,v4,v5,v2_extruded,v3_extruded,v4_extruded,v5_extruded])
    # print left_middle_vertices

    io.add_convex_shape('Left_middle',left_middle_vertices )
    io.add_object('left_middle', [Contactor('Left_middle')],
                 translation = trans)

    ######### left_down
    v6 = numpy.array([4.370*scale, box_width, -.6*scale])
    v7 = numpy.array([(6.900-1.770)*scale, box_width, -.6*scale])

    left_down_normal = normal_plane(v4,v6,v5)
    # print('left_down_normal=', left_down_normal)

    v6_extruded = v6 - [plane_thickness, 0.0, 0.] + numpy.dot(plane_thickness,
                                                              left_down_normal)
    v7_extruded = v7 + [plane_thickness, 0.0, 0.] + numpy.dot(plane_thickness,
                                                              left_down_normal)

    left_down_vertices = numpy.array(
        [v4-[plane_thickness, 0.0, 0.],
         v5+[plane_thickness, 0.0, 0.],
         v6-[plane_thickness, 0.0, 0.],
         v7+[plane_thickness, 0.0, 0.],
         v4_extruded-[plane_thickness, 0.0, 0.],
         v5_extruded+[plane_thickness, 0.0, 0.],
         v6_extruded,v7_extruded])
    # print left_down_vertices

    io.add_convex_shape('Left_down',left_down_vertices )
    io.add_object('left_udown', [Contactor('Left_down')],
                 translation = trans)

    ######### right_up
    v8 = numpy.array([0, box_width, box_height])
    v9 = numpy.array([box_length, box_width, 1.200*scale])

    v10 = numpy.array([(6.900-1.770)*scale, box_width, 0.0])
    v11 =  numpy.array([4.370*scale, box_width, 0.0])

    right_up_normal = normal_plane(v8,v9,v10)
    # print('right_up_normal=', right_up_normal)

    v8_extruded = v8 + numpy.dot(plane_thickness,right_up_normal)
    v9_extruded = v9 + numpy.dot(plane_thickness,right_up_normal)
    v10_extruded = v10 + numpy.dot(plane_thickness,right_up_normal)
    v11_extruded = v11 + numpy.dot(plane_thickness,right_up_normal)

    right_up_vertices = numpy.array(
        [v8,v9,v10,v11,v8_extruded,v9_extruded,v10_extruded,v11_extruded])
    # print right_up_vertices

    io.add_convex_shape('Right_up',right_up_vertices )
    io.add_object('right_up', [Contactor('Right_up')],
                 translation = trans)

    ######### rear_up
    rear_up_normal = normal_plane(v1,v8,v4)
    # print('rear_up_normal=', rear_up_normal)

    v1_extruded = v1 + numpy.dot(plane_thickness,rear_up_normal)
    v2_extruded = v2 + numpy.dot(plane_thickness,rear_up_normal)
    v8_extruded = v8 + numpy.dot(plane_thickness,rear_up_normal)
    v4_extruded = v4 + numpy.dot(plane_thickness,rear_up_normal)
    v11_extruded = v11 + numpy.dot(plane_thickness,rear_up_normal)

    rear_up_vertices = numpy.array(
        [v1-[0.0,plane_thickness,0.0],
         v2-[0.0,plane_thickness,0.0],
         v8+[0.0,plane_thickness,0.0],
         v4-[0.0,plane_thickness,0.0],
         v11+[0.0,plane_thickness,0.0],
         v1_extruded-[0.0,plane_thickness,0.0],
         v2_extruded-[0.0,plane_thickness,0.0],
         v8_extruded+[0.0,plane_thickness,0.0],
         v4_extruded-[0.0,plane_thickness,0.0],
         v11_extruded+[0.0,plane_thickness,0.0]])
    # print rear_up_vertices

    io.add_convex_shape('Rear_up',rear_up_vertices )
    io.add_object('rear_up', [Contactor('Rear_up')],
                 translation = trans)


    ######### rear_down
    #v12 = numpy.array([(6.900-1.770)*scale, box_width,-.6*scale])
    #v13 = numpy.array([4.370*scale, box_width, -.6*scale])


    rear_down_normal = normal_plane(v4,v11,v6)
    # print('rear_down_normal=', rear_down_normal)

    v4_extruded = v4 + numpy.dot(plane_thickness, rear_down_normal)
    v11_extruded = v11 + numpy.dot(plane_thickness, rear_down_normal)
    v6_extruded = v6 + numpy.dot(plane_thickness, rear_down_normal)

    rear_down_vertices=numpy.array([v4,v11,v6,v4_extruded,v11_extruded,v6_extruded])
    # print rear_down_vertices

    io.add_convex_shape('Rear_down',rear_down_vertices )
    io.add_object('rear_down', [Contactor('Rear_down')],
                 translation = trans)

    ######### front_up
    front_up_normal = normal_plane(v3,v5,v9)
    # print('front_up_normal=', front_up_normal)

    v3_extruded = v3 + numpy.dot(plane_thickness,front_up_normal)
    v5_extruded = v5 + numpy.dot(plane_thickness,front_up_normal)
    v9_extruded = v9 + numpy.dot(plane_thickness,front_up_normal)
    v10_extruded = v10 + numpy.dot(plane_thickness,front_up_normal)

    front_up_vertices = numpy.array(
        [v3-[0.0,plane_thickness,0.0],v5-[0.0,plane_thickness,0.0],
         v9+[0.0,plane_thickness,0.0],v10+[0.0,plane_thickness,0.0],
         v3_extruded-[0.0,plane_thickness,0.0],v5_extruded-[0.0,plane_thickness,0.0],
         v9_extruded+[0.0,plane_thickness,0.0],v10_extruded+[0.0,plane_thickness,0.0]])
    # print front_up_vertices

    io.add_convex_shape('Front_up',front_up_vertices )
    io.add_object('front_up', [Contactor('Front_up')],
                 translation = trans)

    ######### front_down
    front_down_normal = normal_plane(v5,v7,v10)
    # print('front_down_normal=', front_down_normal)

    v7_extruded = v7 + numpy.dot(plane_thickness,front_down_normal)
    v5_extruded = v5 + numpy.dot(plane_thickness,front_down_normal)
    v10_extruded = v10 + numpy.dot(plane_thickness,front_down_normal)

    front_down_vertices=numpy.array([v5,v7,v10,v5_extruded,v7_extruded,v10_extruded])
    # print front_down_vertices

    io.add_convex_shape('Front_down',front_down_vertices )
    io.add_object('front_down', [Contactor('Front_down')],
                 translation = trans)
