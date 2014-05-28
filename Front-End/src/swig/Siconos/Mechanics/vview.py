#!/usr/bin/env python
import sys
import os
import vtk
from vtk.util import numpy_support
from math import atan2, pi
import bisect
from numpy.linalg import norm
import numpy
import random
import h5py
import getopt
from contextlib import contextmanager
import tempfile
import numpy as np

from Siconos.Mechanics.ContactDetection.Bullet import btQuaternion
from Siconos.Mechanics.ContactDetection.Bullet.BulletWrap import __mul__ as mul
from Siconos.Mechanics.ContactDetection.Bullet.BulletWrap import quatRotate
#
# a context manager for a *named* temp file
#
@contextmanager
def tmpfile():
    (_, tfilename) = tempfile.mkstemp()
    fid = open(tfilename, 'w')
    yield (fid, tfilename)
    fid.close()
    os.remove(tfilename)

def usage():
    print '{0}: Usage'.format(sys.argv[0])
    print """
    {0} [--help] [tmin=<float value>] [tmax=<float value>]
        [--cf-scale=<float value>] [hdf5 file]
    """


def add_compatiblity_methods(obj):

    if hasattr(obj, 'SetInput'):
        obj.SetInputData = obj.SetInput

    if hasattr(obj, 'AddInput'):
        obj.AddInputData = obj.AddInput

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'dat', 'tmin=', 'tmax=',
                                    'cf-scale='])
except getopt.GetoptError, err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)

min_time = None
max_time = None
scale_factor = 1
for o, a in opts:

    if o == '--help':
        usage()
        exit(0)

    elif o == '--tmin':
        min_time = float(a)

    elif o == '--tmax':
        max_time = float(a)

    elif o == '--cf-scale':
        scale_factor = float(a)

if len(args)>0:
    io_filename = args[0]
else:
    io_filename = 'io.hdf5'

print min_time, max_time
if len(args) == 1:
    filename = args[0]
elif len(args) > 1:
    usage()
    exit(1)


def random_color():
    r = random.uniform(0.1, 0.9)
    g = random.uniform(0.1, 0.9)
    b = random.uniform(0.1, 0.9)
    return r, g, b


def axis_angle(q):
    w, v = q[0], q[1:]
    nv = norm(v)
    theta = 2 * atan2(nv, w)

    if nv != 0.:
        v = [iv / nv for iv in v]
    else:
        v = [0., 0., 0.]
    return v, theta

axis_anglev = numpy.vectorize(axis_angle)
transforms = dict()
contactors = dict()
offsets = dict()

vtkmath = vtk.vtkMath()

class Quaternion():

    def __init__(self, *args):
        self._data = vtk.vtkQuaternion[float](*args)

    def __mul__(self, q):
        r = Quaternion()
        vtkmath.MultiplyQuaternion(self._data, q._data, r._data)
        return r

    def __getitem__(self, i):
        return self._data[i]

    def conjugate(self):
        r = Quaternion((self[0], self[1], self[2], self[3]))
        r._data.Conjugate()
        return r

    def rotate(self, v):
        pv = Quaternion((0, v[0], v[1], v[2]))
        rv = self * pv * self.conjugate()
        assert(rv[0] == 0)
        return [rv[1], rv[2], rv[3]]

    def axisAngle(self):
        r = [0,0,0]
        a = self._data.GetRotationAngleAndAxis(r)
        return r, a


def set_position(instance, q0, q1, q2, q3, q4, q5, q6):


    q = Quaternion((q3, q4, q5, q6))

    for transform, offset in zip(transforms[instance], offsets[instance]):

        p = q.rotate(offset[0])

        r = q * Quaternion(offset[1])

        transform.Identity()
        transform.Translate(q0 + p[0], q1 + p[1], q2 + p[2])

        axis, angle = r.axisAngle()


#        axis, angle = axis_angle((q3, q4, q5, q6))

        transform.RotateWXYZ(angle * 180./pi,
                             axis[0],
                             axis[1],
                             axis[2])

set_positionv = numpy.vectorize(set_position)


def usage():
    print """{0}
    """.format(sys.argv[0])

    #ref_filename = args[0]
    #bind_filename = args[1]
    #pos_filename = args[2]
    #cf_filename = args[3]

ref_filename = 'ref.txt'
cfg_filename = 'shapes.cfg'
bind_filename = 'bindings.dat'
dpos_filename = 'dpos.dat'
spos_filename = 'spos.dat'
cf_filename = 'cf.dat'

refs = []
refs_attrs = []
#with open(ref_filename, 'r') as ref_file:
#    for line in ref_file:
##        line_tokens = shlex.split(line)
#        refs.append(line_tokens[0])
#        refs_attrs.append([float(x) for x in line_tokens[1:]])

shape = dict()

#with open(bind_filename, 'r') as bind_file:
#    for line in bind_file:
#        obj_id, obj_shape = shlex.split(line)
#        shape[int(obj_id)] = int(obj_shape)

pos = dict()
instances = dict()

import numpy


def load():

    inFile = h5py.File(io_filename, 'r')
    spos_data = inFile['data']['static']
    dpos_data = inFile['data']['dynamic']
    try:
        cf_data = inFile['data']['cf'][:].copy()
    except:
        cf_data = None

    solv_data = inFile['data']['solv']

    return spos_data, dpos_data, cf_data, solv_data

spos_data, dpos_data, cf_data, solv_data = load()
#def contact_point_reader():
#    global dos
#    dos.GetOutput().GetFieldData()


#contact_field = vtk.vtkPointData()
#c1 = vtk.vtkPolyData()



#contact_data.SetFieldData(contact_field)

#id_f = numpy.where(cf_data[:,0] == min(cf_data[:,0]))[0]

#cf_provider = vtk.vtkProgrammableSource()
keeper = []

# contact forces provider
class CFprov():

    def __init__(self, data):
        self._data = None
        self._mu_coefs = None
        if data is not None:
            if len(data) > 0:
                self._data = data
                self._mu_coefs = set(self._data[:, 1])
        else:
            self._data = None
            self._mu_coefs = None

        if self._data is not None:
            self._time = min(self._data[:, 0])
        else:
            self._time = 0

        self.cpa_at_time = dict()
        self.cpa = dict()

        self.cpb_at_time = dict()
        self.cpb = dict()

        self.cf_at_time = dict()
        self.cf = dict()

        self.cn_at_time = dict()
        self.cn = dict()

        self._contact_field = dict()
        self._output = dict()

        for mu in self._mu_coefs:
            self._contact_field[mu] = vtk.vtkPointData()
            self._output[mu] = vtk.vtkPolyData()
            self._output[mu].SetFieldData(self._contact_field[mu])


    def method(self):

        if self._data is not None:

            id_f = numpy.where(abs(self._data[:, 0] - self._time) < 1e-15)[0]

            for mu in self._mu_coefs:
                imu = numpy.where(self._data[id_f, 1] == mu)

                self.cpa_at_time[mu] = self._data[id_f[imu], 2:5].copy()
                self.cpb_at_time[mu] = self._data[id_f[imu], 5:8].copy()
                self.cn_at_time[mu] = - self._data[id_f[imu], 8:11].copy()
                self.cf_at_time[mu] = self._data[id_f[imu], 11:14].copy()

                self.cpa[mu] = numpy_support.numpy_to_vtk(self.cpa_at_time[mu])
                self.cpa[mu].SetName('contactPositionsA')

                self.cpb[mu] = numpy_support.numpy_to_vtk(self.cpb_at_time[mu])
                self.cpb[mu].SetName('contactPositionsB')

                self.cn[mu] = numpy_support.numpy_to_vtk(self.cn_at_time[mu])
                self.cn[mu].SetName('contactNormals')

                self.cf[mu] = numpy_support.numpy_to_vtk(self.cf_at_time[mu])
                self.cf[mu].SetName('contactForces')

                self._contact_field[mu].AddArray(self.cpa[mu])
                self._contact_field[mu].AddArray(self.cpb[mu])
                self._contact_field[mu].AddArray(self.cn[mu])
                self._contact_field[mu].AddArray(self.cf[mu])


        else:
            pass

        #self._output.Update()

cf_prov = CFprov(cf_data)


contact_posa = dict()
contact_posb = dict()
contact_pos_force = dict()
contact_pos_norm = dict()

for mu in cf_prov._mu_coefs:

    contact_posa[mu] = vtk.vtkDataObjectToDataSetFilter()
    contact_posb[mu] = vtk.vtkDataObjectToDataSetFilter()

    add_compatiblity_methods(contact_posa[mu])
    add_compatiblity_methods(contact_posb[mu])

    contact_pos_force[mu] = vtk.vtkFieldDataToAttributeDataFilter()
    contact_pos_norm[mu] = vtk.vtkFieldDataToAttributeDataFilter()

    contact_posa[mu].SetDataSetTypeToPolyData()
    contact_posa[mu].SetPointComponent(0, "contactPositionsA", 0)
    contact_posa[mu].SetPointComponent(1, "contactPositionsA", 1)
    contact_posa[mu].SetPointComponent(2, "contactPositionsA", 2)

    contact_posb[mu].SetDataSetTypeToPolyData()
    contact_posb[mu].SetPointComponent(0, "contactPositionsB", 0)
    contact_posb[mu].SetPointComponent(1, "contactPositionsB", 1)
    contact_posb[mu].SetPointComponent(2, "contactPositionsB", 2)

    contact_pos_force[mu].SetInputConnection(contact_posa[mu].GetOutputPort())
    contact_pos_force[mu].SetInputFieldToDataObjectField()
    contact_pos_force[mu].SetOutputAttributeDataToPointData()
    contact_pos_force[mu].SetVectorComponent(0, "contactForces", 0)
    contact_pos_force[mu].SetVectorComponent(1, "contactForces", 1)
    contact_pos_force[mu].SetVectorComponent(2, "contactForces", 2)

    contact_pos_norm[mu].SetInputConnection(contact_posa[mu].GetOutputPort())
    contact_pos_norm[mu].SetInputFieldToDataObjectField()
    contact_pos_norm[mu].SetOutputAttributeDataToPointData()
    contact_pos_norm[mu].SetVectorComponent(0, "contactNormals", 0)
    contact_pos_norm[mu].SetVectorComponent(1, "contactNormals", 1)
    contact_pos_norm[mu].SetVectorComponent(2, "contactNormals", 2)

times = list(set(dpos_data[:, 0]))
times.sort()

ndyna = len(numpy.where(dpos_data[:, 0] == times[0]))

if len(spos_data) > 0:
    nstatic = len(numpy.where(spos_data[:, 0] == times[0]))
#    instances = set(dpos_data[:, 1]).union(set(spos_data[:, 1]))
else:
    nstatic = 0
#    instances = set(dpos_data[:, 1])



cf_prov._time = min(times[:])

cf_prov.method()

cone = dict()
cone_glyph = dict()
cmapper = dict()
cactor = dict()
arrow = dict()
cylinder = dict()
sphere = dict()
arrow_glyph = dict()
gmapper = dict()
gactor = dict()
ctransform = dict()
cylinder_glyph = dict()
clmapper = dict()
sphere_glypha = dict()
sphere_glyphb = dict()
smappera = dict()
smapperb = dict()
sactora = dict()
sactorb = dict()
clactor = dict()

transform = vtk.vtkTransform()
transform.Translate(-0.5, 0., 0.)

for mu in cf_prov._mu_coefs:
    contact_posa[mu].SetInputData(cf_prov._output[mu])
    contact_posa[mu].Update()
    contact_posb[mu].SetInputData(cf_prov._output[mu])
    contact_posb[mu].Update()

    contact_pos_force[mu].Update()
    contact_pos_norm[mu].Update()

    cone[mu] = vtk.vtkConeSource()
    cone[mu].SetResolution(40)

    cone[mu].SetRadius(mu)  # one coef!!

    cone_glyph[mu] = vtk.vtkGlyph3D()
    cone_glyph[mu].SetSourceTransform(transform)

    cone_glyph[mu].SetInputConnection(contact_pos_norm[mu].GetOutputPort())
    cone_glyph[mu].SetSourceConnection(cone[mu].GetOutputPort())
    #cone_glyph.ScalingOn()
    #cone_glyph.SetScaleModeToScaleByVector()
    #cone_glyph.SetRange(0, 100)
    #cone_glyph.ClampingOn()
    cone_glyph[mu]._scale_fact = 1
    cone_glyph[mu].SetScaleFactor(cone_glyph[mu]._scale_fact * scale_factor)
    cone_glyph[mu].SetVectorModeToUseVector()

    cone_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
    cone_glyph[mu].OrientOn()

    cmapper[mu] = vtk.vtkPolyDataMapper()
    cmapper[mu].SetInputConnection(cone_glyph[mu].GetOutputPort())

    cactor[mu] = vtk.vtkActor()
    cactor[mu].GetProperty().SetOpacity(0.4)
    cactor[mu].GetProperty().SetColor(0, 0, 1)
    cactor[mu].SetMapper(cmapper[mu])


    arrow[mu] = vtk.vtkArrowSource()
    arrow[mu].SetTipResolution(40)
    arrow[mu].SetShaftResolution(40)

    cylinder[mu] = vtk.vtkCylinderSource()
    cylinder[mu].SetRadius(.01)
    cylinder[mu].SetHeight(1)

    sphere[mu] = vtk.vtkSphereSource()

    # 1. scale = (scalar value of that particular data index);
    # 2. denominator = Range[1] - Range[0];
    # 3. scale = (scale < Range[0] ? Range[0] : (scale > Range[1] ? Range[1] : scale));
    # 4. scale = (scale - Range[0]) / denominator;
    # 5. scale *= scaleFactor;

    arrow_glyph[mu] = vtk.vtkGlyph3D()
    arrow_glyph[mu].SetInputConnection(contact_pos_force[mu].GetOutputPort())
    arrow_glyph[mu].SetSourceConnection(arrow[mu].GetOutputPort())
    arrow_glyph[mu].ScalingOn()
    arrow_glyph[mu].SetScaleModeToScaleByVector()
    arrow_glyph[mu].SetRange(0, .01)
    arrow_glyph[mu].ClampingOn()
    arrow_glyph[mu]._scale_fact = 5
    arrow_glyph[mu].SetScaleFactor(arrow_glyph[mu]._scale_fact * scale_factor)
    arrow_glyph[mu].SetVectorModeToUseVector()

    arrow_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactForces')
    arrow_glyph[mu].SetInputArrayToProcess(3, 0, 0, 0, 'contactForces')
    arrow_glyph[mu].OrientOn()

    gmapper[mu] = vtk.vtkPolyDataMapper()
    gmapper[mu].SetInputConnection(arrow_glyph[mu].GetOutputPort())
    gmapper[mu].SetScalarModeToUsePointFieldData()
    gmapper[mu].SetColorModeToMapScalars()
    gmapper[mu].ScalarVisibilityOn()
    gmapper[mu].SelectColorArray('contactForces')
    #gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

    gactor[mu] = vtk.vtkActor()
    gactor[mu].SetMapper(gmapper[mu])


    ctransform[mu] = vtk.vtkTransform()
    ctransform[mu].Translate(-0.5, 0, 0)
    ctransform[mu].RotateWXYZ(90, 0, 0, 1)
    cylinder_glyph[mu] = vtk.vtkGlyph3D()
    cylinder_glyph[mu].SetSourceTransform(ctransform[mu])

    cylinder_glyph[mu].SetInputConnection(contact_pos_norm[mu].GetOutputPort())
    cylinder_glyph[mu].SetSourceConnection(cylinder[mu].GetOutputPort())
    cylinder_glyph[mu].SetVectorModeToUseVector()

    cylinder_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
    cylinder_glyph[mu].OrientOn()
    cylinder_glyph[mu]._scale_fact = 1
    cylinder_glyph[mu].SetScaleFactor(cylinder_glyph[mu]._scale_fact * scale_factor)

    clmapper[mu] = vtk.vtkPolyDataMapper()
    clmapper[mu].SetInputConnection(cylinder_glyph[mu].GetOutputPort())


    sphere_glypha[mu] = vtk.vtkGlyph3D()
    sphere_glypha[mu].SetInputConnection(contact_posa[mu].GetOutputPort())
    sphere_glypha[mu].SetSourceConnection(sphere[mu].GetOutputPort())
    sphere_glypha[mu].ScalingOn()
    #sphere_glypha[mu].SetScaleModeToScaleByVector()
    #sphere_glypha[mu].SetRange(-0.5, 2)
    #sphere_glypha[mu].ClampingOn()
    sphere_glypha[mu]._scale_fact = .1
    sphere_glypha[mu].SetScaleFactor(sphere_glypha[mu]._scale_fact * scale_factor)
    #sphere_glypha[mu].SetVectorModeToUseVector()

    sphere_glyphb[mu] = vtk.vtkGlyph3D()
    sphere_glyphb[mu].SetInputConnection(contact_posb[mu].GetOutputPort())
    sphere_glyphb[mu].SetSourceConnection(sphere[mu].GetOutputPort())
    sphere_glyphb[mu].ScalingOn()
    #sphere_glyphb[mu].SetScaleModeToScaleByVector()
    #sphere_glyphb[mu].SetRange(-0.5, 2)
    #sphere_glyphb[mu].ClampingOn()
    sphere_glyphb[mu]._scale_fact = .1
    sphere_glyphb[mu].SetScaleFactor(sphere_glyphb[mu]._scale_fact * scale_factor)
    #sphere_glyphb[mu].SetVectorModeToUseVector()

    #sphere_glyphb[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
    #sphere_glyph.OrientOn()

    smappera[mu] = vtk.vtkPolyDataMapper()
    smappera[mu].SetInputConnection(sphere_glypha[mu].GetOutputPort())
    smapperb[mu] = vtk.vtkPolyDataMapper()
    smapperb[mu].SetInputConnection(sphere_glyphb[mu].GetOutputPort())

    #cmapper.SetScalarModeToUsePointFieldData()
    #cmapper.SetColorModeToMapScalars()
    #cmapper.ScalarVisibilityOn()
    #cmapper.SelectColorArray('contactNormals')
    #gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())


    clactor[mu] = vtk.vtkActor()
    #cactor.GetProperty().SetOpacity(0.4)
    clactor[mu].GetProperty().SetColor(1, 0, 0)
    clactor[mu].SetMapper(clmapper[mu])

    sactora[mu] = vtk.vtkActor()
    sactora[mu].GetProperty().SetColor(1, 0, 0)
    sactora[mu].SetMapper(smappera[mu])

    sactorb[mu] = vtk.vtkActor()
    sactorb[mu].GetProperty().SetColor(0, 1, 0)
    sactorb[mu].SetMapper(smapperb[mu])

#with open(pos_filename, 'r') as pos_file:
#    for line in pos_file:
#        data = shlex.split(line)
#        time = float(data[0])
#        instance = int(data[1])
#        q0 = float(data[2])
#        q1 = float(data[3])
#        q2 = float(data[4])
#        q3 = float(data[5])
#        q4 = float(data[6])
#        q5 = float(data[7])
#        q6 = float(data[8])
#        instances.add(instance)
#        if time in pos:
#            pos[time].append([instance, q0, q1, q2, q3, q4, q5, q6])
#        else:
#            pos[time] = [[instance, q0, q1, q2, q3, q4, q5, q6]]

renderer = vtk.vtkRenderer()
renderer_window = vtk.vtkRenderWindow()
interactor_renderer = vtk.vtkRenderWindowInteractor()

#camera = vtk.vtkCamera()
#camera.SetViewUp(0, 0, -1)
#camera.SetPosition(-221, 40, 204)
#camera.SetFocalPoint(0, 0, 0)
#camera.ComputeViewPlaneNormal()
#camera.SetRoll(180.0)
#camera.Azimuth(80.0)

#renderer.SetBackground(0.85, 0.85, 0.85)
#renderer.SetActiveCamera(camera)
#renderer.ResetCamera()

readers = dict()
mappers = dict()
actors = []

with h5py.File(io_filename, 'r') as io:

    for shape_name in io['data']['ref']:

        if '.vtp' in shape_name:
            with tmpfile() as tmpf:
                tmpf.write(io['data']['ref'][shape_name][:])
                reader = vtk.vtkXMLPolyDataReader()
                reader.SetFileName(tmpf)
                reader.Update()
                readers[shape_name] = reader
        elif 'primitive' not in io['data']['ref'][shape_name].attrs:
            # a convex shape
            points = vtk.vtkPoints()
            convex = vtk.vtkConvexPointSet()
            data = io['data']['ref'][shape_name][:]
            convex.GetPointIds().SetNumberOfIds(data.shape[0])
            for id_, vertice in enumerate(io['data']['ref'][shape_name][:]):
                points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
                convex.GetPointIds().SetId(id_, id_)
            convex_grid = vtk.vtkUnstructuredGrid()
            convex_grid.Allocate(1, 1)
            convex_grid.InsertNextCell(convex.GetCellType(), convex.GetPointIds())
            convex_grid.SetPoints(points)

            source = vtk.vtkProgrammableSource()

            def sourceExec():
                output = source.GetUnstructuredGridOutput()
                output.Allocate(1, 1)
                output.InsertNextCell(convex.GetCellType(), convex.GetPointIds())
                output.SetPoints(points)

            source.SetExecuteMethod(sourceExec)

            readers[shape_name] = source
            mapper = vtk.vtkDataSetMapper()
            add_compatiblity_methods(mapper)
            mapper.SetInputData(convex_grid)
            mappers[shape_name] = mapper

        else:
            primitive = io['data']['ref'][shape_name].attrs['primitive']
            attrs = io['data']['ref'][shape_name][:][0]
            if primitive == 'Sphere':
                source = vtk.vtkSphereSource()
                source.SetRadius(attrs[0])

            elif primitive == 'Cone':
                source = vtk.vtkConeSource()
                source.SetRadius(attrs[0])
                source.SetHeight(attrs[1])
                source.SetResolution(15)
                source.SetDirection(0, 1, 0) # needed

            elif primitive == 'Cylinder':
                source = vtk.vtkCylinderSource()
                source.SetResolution(15)
                source.SetRadius(attrs[0])
                source.SetHeight(attrs[1])
                #           source.SetDirection(0,1,0)

            elif primitive == 'Box':
                source = vtk.vtkCubeSource()
                source.SetXLength(attrs[0])
                source.SetYLength(attrs[1])
                source.SetZLength(attrs[2])

            elif primitive == 'Capsule':
                sphere1 = vtk.vtkSphereSource()
                sphere1.SetRadius(attrs[0])
                sphere1.SetCenter(0, attrs[1] / 2, 0)
                sphere1.SetThetaResolution(15)
                sphere1.SetPhiResolution(15)
                sphere1.Update()

                sphere2 = vtk.vtkSphereSource()
                sphere2.SetRadius(attrs[0])
                sphere2.SetCenter(0, -attrs[1] / 2, 0)
                sphere2.SetThetaResolution(15)
                sphere2.SetPhiResolution(15)
                sphere2.Update()

                cylinder = vtk.vtkCylinderSource()
                cylinder.SetRadius(attrs[0])
                cylinder.SetHeight(attrs[1])
                cylinder.SetResolution(15)
                cylinder.Update()

                data = vtk.vtkMultiBlockDataSet()
                data.SetNumberOfBlocks(3)
                data.SetBlock(0, sphere1.GetOutput())
                data.SetBlock(1, sphere2.GetOutput())
                data.SetBlock(2, cylinder.GetOutput())
                source = vtk.vtkMultiBlockDataGroupFilter()
                add_compatiblity_methods(source)
                source.AddInputData(data)

            readers[shape_name] = source
            mapper = vtk.vtkCompositePolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            mappers[shape_name] = mapper

with h5py.File(io_filename, 'r') as io:
    for instance_name in io['data']['input']:
        instance = int(io['data']['input'][instance_name].attrs['id'])
        contactors[instance] = []
        transforms[instance] = []
        offsets[instance] = []
        for contactor_instance_name in io['data']['input'][instance_name]:
            contactor_name = io['data']['input'][instance_name][contactor_instance_name].attrs['name']
            contactors[instance].append(contactor_name)

            actor = vtk.vtkActor()
            if io['data']['input'][instance_name].attrs['mass'] > 0:
                actor.GetProperty().SetOpacity(0.7)

            actor.GetProperty().SetColor(random_color())
            actor.SetMapper(mappers[contactor_name])
            actors.append(actor)
            renderer.AddActor(actor)
            transform = vtk.vtkTransform()
            actor.SetUserTransform(transform)
            transforms[instance].append(transform)
            offsets[instance].append((io['data']['input'][instance_name][contactor_instance_name].attrs['position'],
                                      io['data']['input'][instance_name][contactor_instance_name].attrs['orientation']))



for mu in cf_prov._mu_coefs:
    renderer.AddActor(gactor[mu])
    renderer.AddActor(cactor[mu])

    renderer.AddActor(clactor[mu])
    renderer.AddActor(sactora[mu])
    renderer.AddActor(sactorb[mu])

import imp
try:
    imp.load_source('myview', 'myview.py')
    import myview
    this_view = myview.MyView(renderer)
except IOError as e:
    pass

id_t0 = numpy.where(dpos_data[:, 0] == min(dpos_data[:, 0]))

if nstatic > 0:
    pos_data = numpy.concatenate((spos_data, dpos_data))
else:
    pos_data = dpos_data[:].copy()

set_positionv(pos_data[id_t0, 1], pos_data[id_t0, 2], pos_data[id_t0, 3],
              pos_data[id_t0, 4], pos_data[id_t0, 5], pos_data[id_t0, 6],
              pos_data[id_t0, 7], pos_data[id_t0, 8])

renderer_window.AddRenderer(renderer)
interactor_renderer.SetRenderWindow(renderer_window)
interactor_renderer.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

# http://www.itk.org/Wiki/VTK/Depth_Peeling ?

# #Use a render window with alpha bits (as initial value is 0 (false) ):
renderer_window.SetAlphaBitPlanes(1)

# # Force to not pick a framebuffer with a multisample buffer ( as initial value is 8):
renderer_window.SetMultiSamples(0)

# # Choose to use depth peeling (if supported) (initial value is 0 (false) )
renderer.SetUseDepthPeeling(1)

# # Set depth peeling parameters.
renderer.SetMaximumNumberOfPeels(100)

# # Set the occlusion ratio (initial value is 0.0, exact image)
renderer.SetOcclusionRatio(0.1)

# callback maker for scale manipulation
def make_scale_observer(glyphs):

    def scale_observer(obj, event):
        slider_repres = obj.GetRepresentation()
        scale_at_pos = slider_repres.GetValue()
        for glyph in glyphs:
            glyph.SetScaleFactor(scale_at_pos * glyph._scale_fact)

    return scale_observer


# make a slider widget and its representation
def make_slider(title, observer, interactor,
                startvalue, minvalue, maxvalue, cx1, cy1, cx2, cy2):
    slider_repres = vtk.vtkSliderRepresentation2D()
    slider_repres.SetMinimumValue(minvalue - (maxvalue - minvalue) / 100)
    slider_repres.SetMaximumValue(maxvalue + (maxvalue - minvalue) / 100)
    slider_repres.SetValue(startvalue)
    slider_repres.SetTitleText(title)
    slider_repres.GetPoint1Coordinate().\
        SetCoordinateSystemToNormalizedDisplay()
    slider_repres.GetPoint1Coordinate().SetValue(cx1, cy1)
    slider_repres.GetPoint2Coordinate().\
        SetCoordinateSystemToNormalizedDisplay()
    slider_repres.GetPoint2Coordinate().SetValue(cx2, cy2)

    slider_repres.SetSliderLength(0.02)
    slider_repres.SetSliderWidth(0.03)
    slider_repres.SetEndCapLength(0.01)
    slider_repres.SetEndCapWidth(0.03)
    slider_repres.SetTubeWidth(0.005)
    slider_repres.SetLabelFormat('%f')
    slider_repres.SetTitleHeight(0.02)
    slider_repres.SetLabelHeight(0.02)

    slider_widget = vtk.vtkSliderWidget()
    slider_widget.SetInteractor(interactor)
    slider_widget.SetRepresentation(slider_repres)
    slider_widget.KeyPressActivationOff()
    slider_widget.SetAnimationModeToAnimate()
    slider_widget.SetEnabled(True)
    slider_widget.AddObserver('InteractionEvent', observer)

    return slider_widget, slider_repres


image_maker = vtk.vtkWindowToImageFilter()
image_maker.SetInput(renderer_window)
writer = vtk.vtkPNGWriter()
writer.SetInputConnection(image_maker.GetOutputPort())

class InputObserver():

    def __init__(self, times, slider_repres):
        self._stimes = set(times)
        self._opacity = 1.0
        self._time_step = (max(self._stimes) - min(self._stimes)) \
                           / len(self._stimes)
        self._time = min(times)
        self._slider_repres = slider_repres
        self._current_id = vtk.vtkIdTypeArray()
        self._renderer = renderer
        self._renderer_window = renderer_window
        self._times = times
        self._image_counter = 0

    def update(self):
        index = bisect.bisect_left(self._times, self._time)
        index = max(0, index)
        index = min(index, len(self._times)-1)
        cf_prov._time = self._times[index]
        cf_prov.method()

#        contact_posa.Update()

#        contact_posb.Update()

        #contact_pos_force.Update()
        # arrow_glyph.Update()
        #gmapper.Update()

        id_t = numpy.where(pos_data[:, 0] == self._times[index])
        set_positionv(pos_data[id_t, 1], pos_data[id_t, 2], pos_data[id_t, 3],
                      pos_data[id_t, 4],
                      pos_data[id_t, 5], pos_data[id_t, 6], pos_data[id_t, 7],
                      pos_data[id_t, 8])

        self._slider_repres.SetValue(self._time)
        renderer_window.Render()

        self._current_id.SetNumberOfValues(1)
        self._current_id.SetValue(0, index)

        self._iter_plot.SetSelection(self._current_id)
        self._prec_plot.SetSelection(self._current_id)
        self._iter_plot_view.Update()
        self._prec_plot_view.Update()
        self._iter_plot_view.GetRenderer().GetRenderWindow().Render()
        self._prec_plot_view.GetRenderer().GetRenderWindow().Render()

    def object_pos(self, id_):
        index = bisect.bisect_left(self._times, self._time)
        index = max(0, index)
        index = min(index, len(self._times)-1)

        id_t = numpy.where(pos_data[:, 0] == self._times[index])
        return (pos_data[id_t[0][id_], 2], pos_data[id_t[0][id_], 3], pos_data[id_t[0][id_], 4])

    def set_opacity(self):
        for instance, actor in zip(instances, actors):
            if instance >= 0:
                actor.GetProperty().SetOpacity(self._opacity)

    def key(self, obj, event):
        key = obj.GetKeySym()
        print 'key', key

        if key == 'r':
            spos_data, dpos_data, cf_data, solv_data = load()
            cf_prov = CFprov(cf_data)
            times = list(set(dpos_data[:, 0]))
            times.sort()
            ndyna = len(numpy.where(dpos_data[:, 0] == times[0]))
            if len(spos_data) > 0:
                nstatic = len(numpy.where(spos_data[:, 0] == times[0]))
                instances = set(dpos_data[:, 1]).union(set(spos_data[:, 1]))
            else:
                nstatic = 0
                instances = set(dpos_data[:, 1])
            cf_prov._time = min(times[:])
            cf_prov.method()
            contact_posa.SetInputData(cf_prov._output)
            contact_posa.Update()
            contact_posb.SetInputData(cf_prov._output)
            contact_posb.Update()
            id_t0 = numpy.where(dpos_data[:, 0] == min(dpos_data[:, 0]))
            contact_pos_force.Update()
            contact_pos_norm.Update()
            if nstatic > 0:
                pos_data = numpy.concatenate((spos_data, dpos_data))
            else:
                pos_data = dpos_data[:].copy()
            min_time = times[0]

            max_time = times[len(times) - 1]

            self._slider_repres.SetMinimumValue(min_time)
            self._slider_repres.SetMaximumValue(max_time)
            self.update()

        if key == 'p':
            self._image_counter += 1
            image_maker.Update()
            writer.SetFileName('vview-{0}.png'.format(self._image_counter))
            writer.Write()

        if key == 'Up':
                self._time_step = self._time_step * 2.
                self._time += self._time_step

        if key == 'Down':
                self._time_step = self._time_step / 2.
                self._time -= self._time_step

        if key == 'Left':
                self._time -= self._time_step

        if key == 'Right':
                self._time += self._time_step

        if key == 't':
                self._opacity -= .1
                self.set_opacity()

        if key == 'T':
                self._opacity += .1
                self.set_opacity()

        if key == 'c':
                print 'camera position:',self._renderer.GetActiveCamera().GetPosition()
                print 'camera focal point', self._renderer.GetActiveCamera().GetFocalPoint()
                print 'camera clipping plane', self._renderer.GetActiveCamera().GetClippingRange()


        if key == 'C':
                this_view.action(self)
        self.update()

    def time(self, obj, event):
        slider_repres = obj.GetRepresentation()
        self._time = slider_repres.GetValue()
        self.update()

        # observer on 2D chart
    def iter_plot_observer(self, obj, event):

        if self._iter_plot.GetSelection() is not None:
            # just one selection at the moment!
            if self._iter_plot.GetSelection().GetMaxId() >= 0:
                self._time = self._times[self._iter_plot.GetSelection().GetValue(0)]
                # -> recompute index ...
                self.update()

    def prec_plot_observer(self, obj, event):
        if self._prec_plot.GetSelection() is not None:
            # just one selection at the moment!
            if self._prec_plot.GetSelection().GetMaxId() >= 0:
                self._time = self._times[self._prec_plot.GetSelection().GetValue(0)]
                # -> recompute index ...
                self.update()

slider_repres = vtk.vtkSliderRepresentation2D()

if min_time is None:
    min_time = times[0]

if max_time is None:
    max_time = times[len(times) - 1]

slider_repres.SetMinimumValue(min_time)
slider_repres.SetMaximumValue(max_time)
slider_repres.SetValue(min_time)
slider_repres.SetTitleText("time")
slider_repres.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
slider_repres.GetPoint1Coordinate().SetValue(0.4, 0.9)
slider_repres.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
slider_repres.GetPoint2Coordinate().SetValue(0.9, 0.9)

slider_repres.SetSliderLength(0.02)
slider_repres.SetSliderWidth(0.03)
slider_repres.SetEndCapLength(0.01)
slider_repres.SetEndCapWidth(0.03)
slider_repres.SetTubeWidth(0.005)
slider_repres.SetLabelFormat("%3.4lf")
slider_repres.SetTitleHeight(0.02)
slider_repres.SetLabelHeight(0.02)

slider_widget = vtk.vtkSliderWidget()
slider_widget.SetInteractor(interactor_renderer)
slider_widget.SetRepresentation(slider_repres)
slider_widget.KeyPressActivationOff()
slider_widget.SetAnimationModeToAnimate()
slider_widget.SetEnabled(True)

input_observer = InputObserver(times, slider_repres)


slider_widget.AddObserver("InteractionEvent", input_observer.time)

interactor_renderer.AddObserver('KeyPressEvent', input_observer.key)
interactor_renderer.AddObserver('KeyPressEvent', input_observer.key)

# Create a vtkLight, and set the light parameters.
light = vtk.vtkLight()
light.SetFocalPoint(0, 0, 0)
light.SetPosition(0, 0, 500)
#light.SetLightTypeToHeadlight()
renderer.AddLight(light)

hlight = vtk.vtkLight()
hlight.SetFocalPoint(0, 0, 0)
#hlight.SetPosition(0, 0, 500)
hlight.SetLightTypeToHeadlight()
renderer.AddLight(hlight)



#import numpy as np
#import matplotlib.pyplot as plt

#plt.plot(solv_data[:,0], solv_data[:,1])
#plt.show()


# Warning! numpy support offer a view on numpy array
# the numpy array must not be garbage collected!
nxtime = solv_data[:, 0].copy()
nxiters = solv_data[:, 1].copy()
nprecs = solv_data[:, 2].copy()
xtime = numpy_support.numpy_to_vtk(nxtime)
xiters = numpy_support.numpy_to_vtk(nxiters)
xprecs = numpy_support.numpy_to_vtk(nprecs)

xtime.SetName('time')
xiters.SetName('iterations')
xprecs.SetName('precisions')

table = vtk.vtkTable()
table.AddColumn(xtime)
table.AddColumn(xiters)
table.AddColumn(xprecs)
#table.Dump()

tview_iter = vtk.vtkContextView()
tview_prec = vtk.vtkContextView()

chart_iter = vtk.vtkChartXY()
chart_prec = vtk.vtkChartXY()
tview_iter.GetScene().AddItem(chart_iter)
tview_prec.GetScene().AddItem(chart_prec)
iter_plot = chart_iter.AddPlot(vtk.vtkChart.LINE)
iter_plot.SetLabel('Solver iterations')
iter_plot.GetXAxis().SetTitle('time')
iter_plot.GetYAxis().SetTitle('iterations')

prec_plot = chart_prec.AddPlot(vtk.vtkChart.LINE)
prec_plot.SetLabel('Solver precisions')
prec_plot.GetXAxis().SetTitle('time')
prec_plot.GetYAxis().SetTitle('precisions')

add_compatiblity_methods(iter_plot)
add_compatiblity_methods(prec_plot)

iter_plot.SetInputData(table, 'time', 'iterations')
prec_plot.SetInputData(table, 'time', 'precisions')
iter_plot.SetWidth(5.0)
prec_plot.SetWidth(5.0)
iter_plot.SetColor(0, 255, 0, 255)
prec_plot.SetColor(0, 255, 0, 255)

input_observer._iter_plot = iter_plot
input_observer._prec_plot = prec_plot
input_observer._iter_plot_view = tview_iter
input_observer._prec_plot_view = tview_prec

tview_iter.GetInteractor().AddObserver('RightButtonReleaseEvent',
                                       input_observer.iter_plot_observer)

tview_prec.GetInteractor().AddObserver('RightButtonReleaseEvent',
                                       input_observer.prec_plot_observer)


#screen_size = renderer_window.GetScreenSize()
renderer_window.SetSize(600, 600)
tview_iter.GetRenderer().GetRenderWindow().SetSize(600, 200)
tview_prec.GetRenderer().GetRenderWindow().SetSize(600, 200)


tview_iter.GetInteractor().Initialize()
#tview_iter.GetInteractor().Start()
tview_iter.GetRenderer().SetBackground(.9, .9, .9)
tview_iter.GetRenderer().Render()

tview_prec.GetInteractor().Initialize()
#tview_prec.GetInteractor().Start()
tview_prec.GetRenderer().SetBackground(.9, .9, .9)
tview_prec.GetRenderer().Render()

slwsc, slrepsc = make_slider('Scale',
                             make_scale_observer([cone_glyph, cylinder_glyph, sphere_glypha, sphere_glyphb, arrow_glyph]
                                                 ),
                                                 interactor_renderer,
                                                 scale_factor, scale_factor - scale_factor / 2,
                                                 scale_factor + scale_factor / 2,
                                                 0.01, 0.01, 0.01, 0.7)



renderer_window.Render()
interactor_renderer.Initialize()
interactor_renderer.Start()
