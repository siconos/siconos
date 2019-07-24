#!/usr/bin/env @PYTHON_EXECUTABLE@
"""
Description: Export a Siconos mechanics-IO HDF5 file in VTK format.
"""

# Lighter imports before command line parsing
from __future__ import print_function
import sys
import os
import getopt

#
# a replacement for vview --vtk-export
#

def usage(long=False):
    print(__doc__); print()
    print('Usage:  {0} [--help] [--version] [--ascii] <HDF5>'
          .format(os.path.split(sys.argv[0])[1]))
    if long:
        print()
        print("""Options:
        --help     display this message
        --version  display version information
        --ascii    export file in ascii format
        """)
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help','version','ascii'])
except getopt.GetoptError as err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)

ascii_mode = False

for o, a in opts:
    if o == '--help':
        usage(long=True)
        exit(0)
    if o == '--version':
        print('{0} @SICONOS_VERSION@'.format(os.path.split(sys.argv[0])[1]))
        exit(0)
    if o in ('--ascii'):
        ascii_mode = True

min_time = None
max_time = None
cf_scale_factor = 1
normalcone_ratio = 1
time_scale_factor = 1
vtk_export_mode = True

if len(args) > 0:
    io_filename = args[0]

else:
    usage()
    exit(1)

# Heavier imports after command line parsing
import vtk
from vtk.util import numpy_support
from math import atan2, pi
import bisect
from numpy.linalg import norm
import numpy
import random
from siconos.io.mechanics_hdf5 import MechanicsHdf5

# attach velocity
# contact points and associated forces are embedded in on a PolyData source

class UnstructuredGridSource(vtk.vtkProgrammableSource):

    def GetOutputPort(self):
        # 3: UnstructuredGridOutput for vtkProgrammableSource
        return vtk.vtkProgrammableSource.GetOutputPort(self, 3)


class ConvexSource(UnstructuredGridSource):

    def __init__(self, convex, points):
        self._convex = convex
        self._points = points
        self.SetExecuteMethod(self.method)

    def method(self):
        output = self.GetUnstructuredGridOutput()
        output.Allocate(1, 1)
        output.InsertNextCell(
            convex.GetCellType(), self._convex.GetPointIds())
        output.SetPoints(self._points)

def add_compatiblity_methods(obj):
    """
    Add missing methods in previous VTK versions.
    """

    if hasattr(obj, 'SetInput'):
        obj.SetInputData = obj.SetInput

    if hasattr(obj, 'AddInput'):
        obj.AddInputData = obj.AddInput

transforms = dict()
transformers = dict()
data_connectors_v = dict()
data_connectors_t = dict()
data_connectors_d = dict()

big_data_source = vtk.vtkMultiBlockDataGroupFilter()
add_compatiblity_methods(big_data_source)

big_data_writer = vtk.vtkXMLMultiBlockDataWriter()
add_compatiblity_methods(big_data_writer)

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
        # assert(rv[0] == 0)
        return [rv[1], rv[2], rv[3]]

    def axisAngle(self):
        r = [0, 0, 0]
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

        transform.RotateWXYZ(angle * 180. / pi,
                             axis[0],
                             axis[1],
                             axis[2])

set_positionv = numpy.vectorize(set_position)


def build_set_velocity(dico):

    def set_velocity(instance, v0, v1, v2, v3, v4, v5):

        if instance in dico:
            dico[instance]._data[:] = [v0, v1, v2, v3, v4, v5]
            dico[instance]._connector.Update()

    set_velocityv = numpy.vectorize(set_velocity)
    return set_velocityv

def build_set_translation(dico):

    def set_translation(instance, x0, x1, x2 ):

        if instance in dico:
            dico[instance]._data[:] = [x0, x1, x2]
            dico[instance]._connector.Update()

    set_translationv = numpy.vectorize(set_translation)
    return set_translationv

def build_set_displacement(dico):

    def set_displacement(instance, x0, x1, x2 ):

        if instance in dico:
            dico[instance]._data[:] = [x0, x1, x2]
            dico[instance]._connector.Update()

    set_displacementv = numpy.vectorize(set_displacement)
    return set_displacementv

def step_reader(step_string):

    from OCC.StlAPI import StlAPI_Writer
    from OCC.STEPControl import STEPControl_Reader
    from OCC.BRep import BRep_Builder
    from OCC.TopoDS import TopoDS_Compound
    from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)

    stl_writer = StlAPI_Writer()
    stl_writer.SetASCIIMode(True)

    with io.tmpfile(contents=io.shapes()[shape_name][:][0]) as tmpfile:
        step_reader = STEPControl_Reader()

        status = step_reader.ReadFile(tmpfile[1])

        if status == IFSelect_RetDone:  # check status
            failsonly = False
            step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
            step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)

            ok = step_reader.TransferRoot(1)
            nbs = step_reader.NbShapes()

            l = []
            for i in range(1, nbs + 1):
                shape = step_reader.Shape(i)

                builder.Add(comp, shape)

            with io.tmpfile(suffix='.stl') as tmpf:
                    stl_writer.Write(comp, tmpf[1])
                    tmpf[0].flush()

                    reader = vtk.vtkSTLReader()
                    reader.SetFileName(tmpf[1])
                    reader.Update()

                    return reader


def brep_reader(brep_string, indx):

    from OCC.StlAPI import StlAPI_Writer

    from OCC.BRepTools import BRepTools_ShapeSet
    shape_set = BRepTools_ShapeSet()
    shape_set.ReadFromString(brep_string)
    shape = shape_set.Shape(shape_set.NbShapes())
    location = shape_set.Locations().Location(indx)
    shape.Location(location)

    stl_writer = StlAPI_Writer()

    with io.tmpfile(suffix='.stl') as tmpf:
        stl_writer.Write(shape, tmpf[1])
        tmpf[0].flush()

        reader = vtk.vtkSTLReader()
        reader.SetFileName(tmpf[1])
        reader.Update()

        return reader


refs = []
refs_attrs = []
shape = dict()

pos = dict()
instances = dict()

with MechanicsHdf5(io_filename=io_filename, mode='r') as io:

    def load():

        ispos_data = io.static_data()
        idpos_data = io.dynamic_data()
        ivelo_data = io.velocities_data()
        icf_data = io.contact_forces_data()[:]

        isolv_data = io.solver_data()

        return ispos_data, idpos_data, ivelo_data, icf_data, isolv_data

    spos_data, dpos_data, velo_data, cf_data, solv_data = load()

    class DataConnector():

        def __init__(self, instance, data_name='velocity', data_size=6):

            self._instance = instance
            self._data_name = data_name
            self._data_size = data_size
            self._connector = vtk.vtkProgrammableFilter()
            self._connector.SetExecuteMethod(self.method)
            self._data = numpy.zeros(data_size)
            self._vtk_data = vtk.vtkFloatArray()
            self._vtk_data.SetName(data_name)
            self._vtk_data.SetNumberOfComponents(data_size)
            self._vtk_data.SetNumberOfTuples(1)

        def method(self):
            input = self._connector.GetInput()
            output = self._connector.GetOutput()
            output.ShallowCopy(input)

            if output.GetFieldData().GetArray(self._data_name) is None:
                output.GetFieldData().AddArray(self._vtk_data)

            data = self._data

            data_t = tuple(data[0:self._data_size])
            output.GetFieldData().GetArray(self._data_name).SetTuple(
                0, data_t)

    # contact forces provider
    class ContactInfoSource():

        def __init__(self, data):
            self._data = None

            if data is not None:
                if len(data) > 0:
                    self._data = data

            else:
                self._data = None

            if self._data is not None:
                self._time = min(self._data[:, 0])
            else:
                self._time = 0

            self._contact_source_a = vtk.vtkProgrammableSource()
            self._contact_source_b = vtk.vtkProgrammableSource()

            self._contact_source_a.SetExecuteMethod(self.method)
            self._contact_source_b.SetExecuteMethod(self.method)

        def method(self):

            # multiblock += contact points
            output_a = self._contact_source_a.GetPolyDataOutput()
            output_b = self._contact_source_b.GetPolyDataOutput()

            id_f = numpy.where(
                abs(self._data[:, 0] - self._time) < 1e-15)[0]

            self.cpa_export = self._data[
                id_f, 2:5].copy()

            self.cpb_export = self._data[
                id_f, 5:8].copy()

            self.cn_export = self._data[
                id_f, 8:11].copy()

            self.cf_export = self._data[
                id_f, 11:14].copy()

            self.cpa_ = numpy_support.numpy_to_vtk(
                self.cpa_export)
            self.cpa_.SetName('contact_positions_A')

            self.cpb_ = numpy_support.numpy_to_vtk(
                self.cpb_export)
            self.cpb_.SetName('contact_positions_B')

            self.cn_ = numpy_support.numpy_to_vtk(
                self.cn_export)
            self.cn_.SetName('contact_normals')

            self.cf_ = numpy_support.numpy_to_vtk(
                self.cf_export)
            self.cf_.SetName('contact_forces')

            output_a.Allocate(len(self.cpa_export), 1)

            cpa_points = vtk.vtkPoints()
            cpa_points.SetNumberOfPoints(len(self.cpa_export))
            cpa_points.SetData(self.cpa_)
            output_a.SetPoints(cpa_points)

            # normal and forces are attached to A points
            output_a.GetPointData().AddArray(self.cn_)
            output_a.GetPointData().AddArray(self.cf_)

            output_b.Allocate(len(self.cpb_export), 1)

            cpb_points = vtk.vtkPoints()
            cpb_points.SetNumberOfPoints(len(self.cpb_export))
            cpb_points.SetData(self.cpb_)

            output_b.SetPoints(cpb_points)

    #  Step 2
    #
    #
    readers = dict()
    vtk_reader = {'vtp': vtk.vtkXMLPolyDataReader,
                  'stl': vtk.vtkSTLReader}

    for shape_name in io.shapes():

        shape_type = io.shapes()[shape_name].attrs['type']

        if shape_type in ['vtp', 'stl']:
            with io.tmpfile() as tmpf:
                tmpf[0].write(str(io.shapes()[shape_name][:][0]))
                tmpf[0].flush()
                reader = vtk_reader[shape_type]()
                reader.SetFileName(tmpf[1])
                reader.Update()
                readers[shape_name] = reader

        elif shape_type in ['brep']:
            # try to find an associated shape
            if 'associated_shape' in io.shapes()[shape_name].attrs:
                associated_shape = \
                    io.shapes()[shape_name].\
                    attrs['associated_shape']

            else:
                if 'brep' in io.shapes()[shape_name].attrs:
                    brep = io.shapes()[shape_name].attrs['brep']
                else:
                    brep = shape_name

                reader = brep_reader(str(io.shapes()[brep][:][0]),
                                     io.shapes()[brep].attrs['occ_indx'])
                readers[shape_name] = reader

        elif shape_type in ['stp', 'step']:
            # try to find an associated shape
            if 'associated_shape' in io.shapes()[shape_name].attrs:
                associated_shape = \
                    io.shapes()[shape_name].\
                    attrs['associated_shape']
                # delayed

            else:
                reader = step_reader(str(io.shapes()[shape_name][:]))

                readers[shape_name] = reader

        elif shape_type == 'convex':
            # a convex shape
            points = vtk.vtkPoints()
            convex = vtk.vtkConvexPointSet()
            data = io.shapes()[shape_name][:]
            convex.GetPointIds().SetNumberOfIds(data.shape[0])

            for id_, vertice in enumerate(io.shapes()[shape_name][:]):
                points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
                convex.GetPointIds().SetId(id_, id_)

            readers[shape_name] = ConvexSource(convex, points)

        else:
            assert shape_type == 'primitive'
            primitive = io.shapes()[shape_name].attrs['primitive']
            attrs = io.shapes()[shape_name][:][0]
            if primitive == 'Sphere':
                source = vtk.vtkSphereSource()
                source.SetRadius(attrs[0])

            elif primitive == 'Cone':
                source = vtk.vtkConeSource()
                source.SetRadius(attrs[0])
                source.SetHeight(attrs[1])
                source.SetResolution(15)
                source.SetDirection(0, 1, 0)  # needed

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

    for instance_name in io.instances():

        instance = int(io.instances()[instance_name].attrs['id'])
        contactors[instance] = []
        transforms[instance] = []
        offsets[instance] = []
        for contactor_instance_name in io.instances()[instance_name]:
            contactor_name = io.instances()[instance_name][
                contactor_instance_name].attrs['name']

            contactors[instance].append(contactor_name)

            transform = vtk.vtkTransform()
            transformer = vtk.vtkTransformFilter()

            if contactor_name in readers:
                transformer.SetInputConnection(
                    readers[contactor_name].GetOutputPort())
            else:
                print ('WARNING: cannot find a shape source for instance:',
                       instance)

            transformer.SetTransform(transform)
            transformers[contactor_name] = transformer

            data_connectors_v[instance] = DataConnector(instance)
            data_connectors_v[instance]._connector.SetInputConnection(
                transformer.GetOutputPort())
            data_connectors_v[instance]._connector.Update()
            big_data_source.AddInputConnection(
                data_connectors_v[instance]._connector.GetOutputPort())

            data_connectors_t[instance] = DataConnector(instance, data_name='translation', data_size=3)
            data_connectors_t[instance]._connector.SetInputConnection(
                transformer.GetOutputPort())
            data_connectors_t[instance]._connector.Update()
            big_data_source.AddInputConnection(
                data_connectors_t[instance]._connector.GetOutputPort())

            data_connectors_d[instance] = DataConnector(instance, data_name='displacement', data_size=3)
            data_connectors_d[instance]._connector.SetInputConnection(
                transformer.GetOutputPort())
            data_connectors_d[instance]._connector.Update()
            big_data_source.AddInputConnection(
                data_connectors_d[instance]._connector.GetOutputPort())

            transforms[instance].append(transform)
            offsets[instance].append(
                (io.instances()[
                    instance_name][
                 contactor_instance_name].attrs['translation'],
                    io.instances()[instance_name][contactor_instance_name].attrs['orientation']))

    pos_data = dpos_data[:].copy()
    spos_data = spos_data[:].copy()
    velo_data = velo_data[:].copy()

    set_velocityv = build_set_velocity(data_connectors_v)
    set_translationv = build_set_translation(data_connectors_t)
    set_displacementv = build_set_displacement(data_connectors_d)

    times = list(set(dpos_data[:, 0]))
    times.sort()

    contact_info_source = ContactInfoSource(cf_data)

    pveloa = DataConnector(0)
    pvelob = DataConnector(0)

    pveloa._connector.SetInputConnection(
        contact_info_source._contact_source_a.GetOutputPort())
    pvelob._connector.SetInputConnection(
        contact_info_source._contact_source_a.GetOutputPort())

    big_data_source.AddInputConnection(
        pveloa._connector.GetOutputPort())
    big_data_source.AddInputConnection(
        pvelob._connector.GetOutputPort())

    big_data_writer.SetInputConnection(big_data_source.GetOutputPort())
    ntime = len(times)
    k=0
    packet= int(ntime/100)+1
    for time in times:
        k=k+1
        if (k%packet == 0):
            sys.stdout.write('.')
        index = bisect.bisect_left(times, time)
        index = max(0, index)
        index = min(index, len(times) - 1)

        contact_info_source._time = times[index]

        # fix: should be called by contact_source?
        contact_info_source.method()

        id_t = numpy.where(pos_data[:, 0] == times[index])

        if numpy.shape(spos_data)[0] > 0:
            set_positionv(spos_data[:, 1], spos_data[:, 2],
                          spos_data[:, 3],
                          spos_data[:, 4], spos_data[:, 5],
                          spos_data[:, 6],
                          spos_data[:, 7], spos_data[:, 8])

        set_positionv(
            pos_data[id_t, 1], pos_data[id_t, 2], pos_data[id_t, 3],
            pos_data[id_t, 4], pos_data[id_t, 5], pos_data[id_t, 6],
            pos_data[id_t, 7], pos_data[id_t, 8])

        id_tv = numpy.where(velo_data[:, 0] == times[index])

        set_velocityv(
            velo_data[id_tv, 1],
            velo_data[id_tv, 2],
            velo_data[id_tv, 3],
            velo_data[id_tv, 4],
            velo_data[id_tv, 5],
            velo_data[id_tv, 6],
            velo_data[id_tv, 7])

        set_translationv(
            pos_data[id_t, 1],
            pos_data[id_t, 2],
            pos_data[id_t, 3],
            pos_data[id_t, 4],
        )

        # set_displacementv(
        #     pos_data[id_t, 1],
        #     pos_data[id_t, 2]- pos_data[0, 2],
        #     pos_data[id_t, 3]- pos_data[0, 3],
        #     pos_data[id_t, 4]- pos_data[0, 4]
        # ) # should be w.r.t initial position

        big_data_writer.SetFileName('{0}-{1}.{2}'.format(os.path.splitext(
                                    os.path.basename(io_filename))[0],
            index, big_data_writer.GetDefaultFileExtension()))
        big_data_writer.SetTimeStep(times[index])
        big_data_source.Update()

        if ascii_mode:
            big_data_writer.SetDataModeToAscii()
        big_data_writer.Write()
    print(' ')
