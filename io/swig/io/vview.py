#!/usr/bin/env @PYTHON_EXECUTABLE@
import sys, os, json
import vtk
from vtk.util import numpy_support
from math import pi
import bisect
from numpy.linalg import norm
import numpy
import random

import getopt

from siconos.io.mechanics_io import Quaternion, Hdf5
from siconos.io.mechanics_io import tmpfile as io_tmpfile

## Persistent configuration
config = {'window_size': [600,600]}
config_fn = os.path.join(os.environ['HOME'], '.config', 'siconos_vview.json')
should_save_config = True

def load_configuration():
    global should_save_config
    if os.path.exists(config_fn):
        try:
            config.update(json.load(open(config_fn)))
            should_save_config = True
        except:
            should_save_config = False
            print("Warning: Error loading configuration `{}'".format(config_fn))

def save_configuration():
    try:
        if not os.path.exists(os.path.join(os.environ['HOME'], '.config')):
            os.mkdir(os.path.join(os.environ['HOME'], '.config'))
        json.dump(config, open(config_fn,'w'))
    except:
        print("Error saving configuration `{}'".format(config_fn))

# Load it immediately
load_configuration()

## Print usage information
def usage():
    print('{0}: Usage'.format(sys.argv[0]))
    print("""
    {0} [--help] [tmin=<float value>] [tmax=<float value>]
        [--cf-scale=<float value>] [--no-cf]
        [--advance=<'fps' or float value>] [--fps=float value]
        [--camera=x,y,z] [--lookat=x,y,z] [--up=x,y,z] [--ortho=scale]
        <hdf5 file>
    """)
    print("""
    Options :
      --help
        display this message
     --tmin= value
       set the time lower bound for visualization
     --tmax= value
       set the time upper bound for visualization
     --cf-scale= value  (default : 1.0 )
       rescale the arrow representing the contact forces by the value.
       the normal cone and the contact points are also rescaled
     --no-cf
       do not display contact forces
     --normalcone-ratio = value  (default : 1.0 )
       introduce a ratio between the representation of the contact forces arrows
       the normal cone and the contact points. useful when the contact forces are
       small with respect to the characteristic dimesion
     --advance= value or 'fps'
       automatically advance time during recording (default : don't advance)
     --fps= value
       frames per second of generated video (default 25)
     --camera=x,y,z
       initial position of the camera (default=above looking down)
     --lookat=x,y,z
       initial direction to look (default=center of bounding box)
     --up=x,y,z
       initial up direction of the camera (default=y-axis)
     --ortho=scale
       start in ortho mode with given parallel scale (default=perspective)
    """)


def add_compatiblity_methods(obj):
    """
    Add missing methods in previous VTK versions.
    """

    if hasattr(obj, 'SetInput'):
        obj.SetInputData = obj.SetInput

    if hasattr(obj, 'AddInput'):
        obj.AddInputData = obj.AddInput

## Parse command line
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help', 'dat', 'tmin=', 'tmax=', 'no-cf',
                                    'cf-scale=', 'normalcone-ratio=',
                                    'advance=', 'fps=', 'camera=', 'lookat=',
                                    'up=', 'ortho='])
except getopt.GetoptError as err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)

min_time = None
max_time = None
cf_scale_factor = 1
normalcone_ratio = 1
time_scale_factor = 1
view_cycle = -1
advance_by_time = None
frames_per_second = 25
cf_disable = False
initial_camera = [None] * 4

for o, a in opts:

    if o == '--help':
        usage()
        exit(0)

    elif o == '--tmin':
        min_time = float(a)

    elif o == '--tmax':
        max_time = float(a)

    elif o == '--cf-scale':
        cf_scale_factor = float(a)

    elif o == '--no-cf':
        cf_disable = True

    elif o == '--normalcone-ratio':
        normalcone_ratio = float(a)

    elif o == '--advance':
        if 'fps' in a:
            advance_by_time = eval(a, {'fps': 1.0 / frames_per_second})
        else:
            advance_by_time = float(a)

    elif o == '--fps':
        frames_per_second = int(a)

    elif o == '--camera':
        initial_camera[0] = map(float, a.split(','))

    elif o == '--lookat':
        initial_camera[1] = map(float, a.split(','))

    elif o == '--up':
        initial_camera[2] = map(float, a.split(','))

    elif o == '--ortho':
        initial_camera[3] = float(a)

if frames_per_second == 0:
    frames_per_second = 25

if len(args) > 0:
    io_filename = args[0]

else:
    usage()
    exit(1)

## Utilities

def random_color():
    r = random.uniform(0.1, 0.9)
    g = random.uniform(0.1, 0.9)
    b = random.uniform(0.1, 0.9)
    return r, g, b

transforms = dict()
transformers = dict()

big_data_source = vtk.vtkMultiBlockDataGroupFilter()
add_compatiblity_methods(big_data_source)

big_data_writer = vtk.vtkXMLMultiBlockDataWriter()
add_compatiblity_methods(big_data_writer)
big_data_writer.SetInputConnection(big_data_source.GetOutputPort())

contactors = dict()
offsets = dict()


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

    with io_tmpfile(contents=io.shapes()[shape_name][:][0]) as tmpfile:
        step_reader = STEPControl_Reader()

        status = step_reader.ReadFile(tmpfile[1])

        if status == IFSelect_RetDone:  # check status
            failsonly = False
            step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
            step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
            ok = step_reader.TransferRoot(1)
            nbs = step_reader.NbShapes()

            for i in range(1, nbs + 1):
                shape = step_reader.Shape(i)

                builder.Add(comp, shape)

            with io_tmpfile(suffix='.stl') as tmpf:
                    stl_writer.Write(comp, tmpf[1])
                    tmpf[0].flush()

                    reader = vtk.vtkSTLReader()
                    reader.SetFileName(tmpf[1])
                    reader.Update()

                    return reader


def iges_reader(iges_string):

    from OCC.StlAPI import StlAPI_Writer
    from OCC.IGESControl import IGESControl_Reader
    from OCC.BRep import BRep_Builder
    from OCC.TopoDS import TopoDS_Compound
    from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)

    stl_writer = StlAPI_Writer()
    stl_writer.SetASCIIMode(True)

    with io_tmpfile(contents=io.shapes()[shape_name][:][0]) as tmpfile:
        iges_reader = IGESControl_Reader()

        status = iges_reader.ReadFile(tmpfile[1])

        if status == IFSelect_RetDone:  # check status
            failsonly = False
            iges_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
            iges_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
            ok = iges_reader.TransferRoots()
            nbs = iges_reader.NbShapes()

            for i in range(1, nbs + 1):
                shape = iges_reader.Shape(i)

                builder.Add(comp, shape)

            with io_tmpfile(suffix='.stl') as tmpf:
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

    with io_tmpfile(suffix='.stl') as tmpf:
        stl_writer.Write(shape, tmpf[1])
        tmpf[0].flush()

        reader = vtk.vtkSTLReader()
        reader.SetFileName(tmpf[1])
        reader.Update()

        return reader


## Program starts

# Read file and open VTK interaction window

refs = []
refs_attrs = []
shape = dict()

pos = dict()
instances = dict()

with Hdf5(io_filename=io_filename, mode='r') as io:

    def load():

        ispos_data = io.static_data()
        idpos_data = io.dynamic_data()
        try:
            idom_data = io.domains_data()
        except ValueError:
            idom_data = None

        icf_data = io.contact_forces_data()[:]

        isolv_data = io.solver_data()

        return ispos_data, idpos_data, idom_data, icf_data, isolv_data

    spos_data, dpos_data, dom_data, cf_data, solv_data = load()

    # contact forces provider
    class CFprov():

        def __init__(self, data, dom_data):
            self._data = None
            self._datap = numpy.array(
                [[1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.]])
            self._mu_coefs = []
            if data is not None:
                if len(data) > 0:
                    self._data = data
                    self._mu_coefs = set(self._data[:, 1])
            else:
                self._data = None
                self._mu_coefs = []

            self._dom_data = dom_data

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

            self.dom_at_time = dict()
            self.dom = dict()

            self._contact_field = dict()
            self._output = dict()

            for mu in self._mu_coefs:
                self._contact_field[mu] = vtk.vtkPointData()
                self._output[mu] = vtk.vtkPolyData()
                self._output[mu].SetFieldData(self._contact_field[mu])

        def xmethod(self):

            if self._data is not None:

                id_f = numpy.where(
                    abs(self._data[:, 0] - self._time) < 1e-15)[0]

                dom_id_f = None
                if self._dom_data is not None:
                    dom_id_f = numpy.where(
                        abs(self._dom_data[:, 0] - self._time) < 1e-15)[0]

                for mu in self._mu_coefs:

                    try:
                        imu = numpy.where(
                            abs(self._data[id_f, 1] - mu) < 1e-15)[0]

                        dom_imu = None
                        if dom_id_f is not None:
                            dom_imu = numpy.where(
                                self._dom_data[dom_id_f,-1] == self._data[id_f[imu],-1]
                            )[0]

                        self.cpa_at_time[mu] = self._data[
                            id_f[imu], 2:5]
                        self.cpb_at_time[mu] = self._data[
                            id_f[imu], 5:8]
                        self.cn_at_time[mu] = - self._data[
                            id_f[imu], 8:11]
                        self.cf_at_time[mu] = self._data[
                            id_f[imu], 11:14]

                        self.dom_at_time[mu] = None
                        if dom_imu is not None:
                            self.dom_at_time[mu] = self._dom_data[
                                dom_id_f[imu], 1]

                        self.cpa[mu] = numpy_support.numpy_to_vtk(
                            self.cpa_at_time[mu])
                        self.cpa[mu].SetName('contactPositionsA')

                        self.cpb[mu] = numpy_support.numpy_to_vtk(
                            self.cpb_at_time[mu])
                        self.cpb[mu].SetName('contactPositionsB')

                        self.cn[mu] = numpy_support.numpy_to_vtk(
                            self.cn_at_time[mu])
                        self.cn[mu].SetName('contactNormals')

                        self.cf[mu] = numpy_support.numpy_to_vtk(
                            self.cf_at_time[mu])
                        self.cf[mu].SetName('contactForces')

                        if self.dom_at_time[mu] is not None:
                            self.dom[mu] = numpy_support.numpy_to_vtk(
                                self.dom_at_time[mu])
                            self.dom[mu].SetName('domains')

                        self._contact_field[mu].AddArray(self.cpa[mu])
                        self._contact_field[mu].AddArray(self.cpb[mu])
                        self._contact_field[mu].AddArray(self.cn[mu])
                        self._contact_field[mu].AddArray(self.cf[mu])
                        if self.dom[mu] is not None:
                            self._contact_field[mu].AddArray(self.dom[mu])
                    except:
                        pass

            else:
                pass

            # self._output.Update()

    contact_posa = dict()
    contact_posb = dict()
    contact_pos_force = dict()
    contact_pos_norm = dict()

    def init_contact_pos(mu):

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

        contact_pos_force[mu].SetInputConnection(
            contact_posa[mu].GetOutputPort())
        contact_pos_force[mu].SetInputFieldToDataObjectField()
        contact_pos_force[mu].SetOutputAttributeDataToPointData()
        contact_pos_force[mu].SetVectorComponent(0, "contactForces", 0)
        contact_pos_force[mu].SetVectorComponent(1, "contactForces", 1)
        contact_pos_force[mu].SetVectorComponent(2, "contactForces", 2)

        contact_pos_norm[mu].SetInputConnection(
            contact_posa[mu].GetOutputPort())
        contact_pos_norm[mu].SetInputFieldToDataObjectField()
        contact_pos_norm[mu].SetOutputAttributeDataToPointData()
        contact_pos_norm[mu].SetVectorComponent(0, "contactNormals", 0)
        contact_pos_norm[mu].SetVectorComponent(1, "contactNormals", 1)
        contact_pos_norm[mu].SetVectorComponent(2, "contactNormals", 2)

        if cf_prov.dom_at_time is not None:
            contact_pos_norm[mu].SetScalarComponent(0, "domains", 0)

    cf_prov = None
    if not cf_disable:
        cf_prov = CFprov(cf_data, dom_data)
        for mu in cf_prov._mu_coefs:
            init_contact_pos(mu)

    times = dpos_data[:, 0]

    if (len(times) == 0):
        print('No dynamic data found!  Empty simulation.')

    if cf_prov is not None:
        cf_prov._time = min(times[:])
        cf_prov.xmethod()

    cone = dict()
    cone_glyph = dict()
    cmapper = dict()
    cLUT = dict()
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
    times_of_birth = dict()

    transform = vtk.vtkTransform()
    transform.Translate(-0.5, 0., 0.)

    def init_cf_sources(mu):
        contact_posa[mu].SetInputData(cf_prov._output[mu])
        contact_posa[mu].Update()
        contact_posb[mu].SetInputData(cf_prov._output[mu])
        contact_posb[mu].Update()

        contact_pos_force[mu].Update()
        contact_pos_norm[mu].Update()

        big_data_source.AddInputConnection(contact_posa[mu].GetOutputPort())
        big_data_source.AddInputConnection(contact_posb[mu].GetOutputPort())
        big_data_source.AddInputConnection(
            contact_pos_force[mu].GetOutputPort())
        big_data_source.AddInputConnection(
            contact_pos_norm[mu].GetOutputPort())

        cone[mu] = vtk.vtkConeSource()
        cone[mu].SetResolution(40)

        cone[mu].SetRadius(mu)  # one coef!!

        cone_glyph[mu] = vtk.vtkGlyph3D()
        cone_glyph[mu].SetSourceTransform(transform)

        cone_glyph[mu].SetInputConnection(contact_pos_norm[mu].GetOutputPort())
        cone_glyph[mu].SetSourceConnection(cone[mu].GetOutputPort())

        cone_glyph[mu]._scale_fact = normalcone_ratio
        cone_glyph[mu].SetScaleFactor(
            cone_glyph[mu]._scale_fact *cf_scale_factor)
        cone_glyph[mu].SetVectorModeToUseVector()

        cone_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
        cone_glyph[mu].OrientOn()

        # Don't allow scalar to affect size of glyph
        cone_glyph[mu].SetScaleModeToDataScalingOff()

        # Allow scalar to affect color of glyph
        cone_glyph[mu].SetColorModeToColorByScalar()

        cmapper[mu] = vtk.vtkPolyDataMapper()
        cmapper[mu].SetInputConnection(cone_glyph[mu].GetOutputPort())

        # Random color map, up to 256 domains
        cLUT[mu] = vtk.vtkLookupTable()
        cLUT[mu].SetNumberOfColors(256)
        cLUT[mu].Build()
        for i in range(256):
            cLUT[mu].SetTableValue(i, *random_color())
        cLUT[mu].SetTableRange(0, 255)

        # By default don't allow scalars to have an effect
        cmapper[mu].ScalarVisibilityOff()

        # If domain information is available, we turn on the color
        # table and turn on scalars
        if cf_prov.dom_at_time is not None:
            cmapper[mu].SetLookupTable(cLUT[mu])
            cmapper[mu].SetColorModeToMapScalars()
            cmapper[mu].SetScalarModeToUsePointData()
            cmapper[mu].SetScalarRange(0,255)
            cmapper[mu].ScalarVisibilityOn()

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
        arrow_glyph[mu].SetInputConnection(
            contact_pos_force[mu].GetOutputPort())
        arrow_glyph[mu].SetSourceConnection(arrow[mu].GetOutputPort())
        arrow_glyph[mu].ScalingOn()
        arrow_glyph[mu].SetScaleModeToScaleByVector()
        arrow_glyph[mu].SetRange(0, .01)
        arrow_glyph[mu].ClampingOn()
        arrow_glyph[mu]._scale_fact = 5
        arrow_glyph[mu].SetScaleFactor(
            arrow_glyph[mu]._scale_fact * cf_scale_factor)
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
        # gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

        gactor[mu] = vtk.vtkActor()
        gactor[mu].SetMapper(gmapper[mu])

        ctransform[mu] = vtk.vtkTransform()
        ctransform[mu].Translate(-0.5, 0, 0)
        ctransform[mu].RotateWXYZ(90, 0, 0, 1)
        cylinder_glyph[mu] = vtk.vtkGlyph3D()
        cylinder_glyph[mu].SetSourceTransform(ctransform[mu])

        cylinder_glyph[mu].SetInputConnection(
            contact_pos_norm[mu].GetOutputPort())
        cylinder_glyph[mu].SetSourceConnection(cylinder[mu].GetOutputPort())
        cylinder_glyph[mu].SetVectorModeToUseVector()

        cylinder_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
        cylinder_glyph[mu].OrientOn()
        cylinder_glyph[mu]._scale_fact = normalcone_ratio
        cylinder_glyph[mu].SetScaleFactor(
            cylinder_glyph[mu]._scale_fact * cf_scale_factor)

        clmapper[mu] = vtk.vtkPolyDataMapper()
        clmapper[mu].SetInputConnection(cylinder_glyph[mu].GetOutputPort())

        sphere_glypha[mu] = vtk.vtkGlyph3D()
        sphere_glypha[mu].SetInputConnection(contact_posa[mu].GetOutputPort())
        sphere_glypha[mu].SetSourceConnection(sphere[mu].GetOutputPort())
        sphere_glypha[mu].ScalingOn()
        # sphere_glypha[mu].SetScaleModeToScaleByVector()
        # sphere_glypha[mu].SetRange(-0.5, 2)
        # sphere_glypha[mu].ClampingOn()
        sphere_glypha[mu]._scale_fact = .1 * normalcone_ratio
        sphere_glypha[mu].SetScaleFactor(
            sphere_glypha[mu]._scale_fact * cf_scale_factor)
        # sphere_glypha[mu].SetVectorModeToUseVector()

        sphere_glyphb[mu] = vtk.vtkGlyph3D()
        sphere_glyphb[mu].SetInputConnection(contact_posb[mu].GetOutputPort())
        sphere_glyphb[mu].SetSourceConnection(sphere[mu].GetOutputPort())
        sphere_glyphb[mu].ScalingOn()
        # sphere_glyphb[mu].SetScaleModeToScaleByVector()
        # sphere_glyphb[mu].SetRange(-0.5, 2)
        # sphere_glyphb[mu].ClampingOn()
        sphere_glyphb[mu]._scale_fact = .1 * normalcone_ratio
        sphere_glyphb[mu].SetScaleFactor(
            sphere_glyphb[mu]._scale_fact * cf_scale_factor)
        # sphere_glyphb[mu].SetVectorModeToUseVector()

        # sphere_glyphb[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
        # sphere_glyph.OrientOn()

        smappera[mu] = vtk.vtkPolyDataMapper()
        smappera[mu].SetInputConnection(sphere_glypha[mu].GetOutputPort())
        smapperb[mu] = vtk.vtkPolyDataMapper()
        smapperb[mu].SetInputConnection(sphere_glyphb[mu].GetOutputPort())

        # cmapper.SetScalarModeToUsePointFieldData()
        # cmapper.SetColorModeToMapScalars()
        # cmapper.ScalarVisibilityOn()
        # cmapper.SelectColorArray('contactNormals')
        # gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

        clactor[mu] = vtk.vtkActor()
        # cactor.GetProperty().SetOpacity(0.4)
        clactor[mu].GetProperty().SetColor(1, 0, 0)
        clactor[mu].SetMapper(clmapper[mu])

        sactora[mu] = vtk.vtkActor()
        sactora[mu].GetProperty().SetColor(1, 0, 0)
        sactora[mu].SetMapper(smappera[mu])

        sactorb[mu] = vtk.vtkActor()
        sactorb[mu].GetProperty().SetColor(0, 1, 0)
        sactorb[mu].SetMapper(smapperb[mu])

    if cf_prov is not None:
        [init_cf_sources(mu) for mu in cf_prov._mu_coefs]

    renderer = vtk.vtkRenderer()
    renderer_window = vtk.vtkRenderWindow()
    interactor_renderer = vtk.vtkRenderWindowInteractor()

    readers = dict()
    datasets = dict()
    mappers = dict()
    dynamic_actors = dict()
    static_actors = dict()
    vtk_reader = {'vtp': vtk.vtkXMLPolyDataReader,
                  'stl': vtk.vtkSTLReader}

    for shape_name in io.shapes():

        shape_type = io.shapes()[shape_name].attrs['type']
        scale = None
        if 'scale' in io.shapes()[shape_name].attrs:
            scale = io.shapes()[shape_name].attrs['scale']

        if shape_type in ['vtp', 'stl']:
            with io_tmpfile() as tmpf:
                tmpf[0].write(str(io.shapes()[shape_name][:][0]))
                tmpf[0].flush()
                reader = vtk_reader[shape_type]()
                reader.SetFileName(tmpf[1])
                reader.Update()
                readers[shape_name] = reader

                # a try for smooth rendering but it does not work here
                normals = vtk.vtkPolyDataNormals()
                normals.SetInputConnection(reader.GetOutputPort())
                normals.SetFeatureAngle(60.0)

                mapper = vtk.vtkDataSetMapper()
                add_compatiblity_methods(mapper)
                mapper.SetInputConnection(normals.GetOutputPort())
                mapper.ScalarVisibilityOff()
                # delayed (see the one in brep)
                # note: "lambda : mapper" fails (dynamic scope)
                # and (x for x in [mapper]) is ok.
                mappers[shape_name] = (x for x in [mapper])

        elif shape_type in ['brep']:
            # try to find an associated shape
            if 'associated_shape' in io.shapes()[shape_name].attrs:
                associated_shape = \
                    io.shapes()[shape_name].\
                    attrs['associated_shape']
                # delayed
                mappers[shape_name] = (x for x in
                                       [mappers[associated_shape]()])
            else:
                if 'brep' in io.shapes()[shape_name].attrs:
                    brep = io.shapes()[shape_name].attrs['brep']
                else:
                    brep = shape_name

                reader = brep_reader(str(io.shapes()[brep][:][0]),
                                     io.shapes()[brep].attrs['occ_indx'])
                readers[shape_name] = reader
                mapper = vtk.vtkDataSetMapper()
                add_compatiblity_methods(mapper)
                mapper.SetInputConnection(reader.GetOutputPort())
                mappers[shape_name] = (x for x in [mapper])

        elif shape_type in ['stp', 'step', 'igs', 'iges']:
            # try to find an associated shape
            if 'associated_shape' in io.shapes()[shape_name].attrs:
                associated_shape = \
                    io.shapes()[shape_name].\
                    attrs['associated_shape']
                # delayed
                mappers[shape_name] = (
                    x for x in [mappers[associated_shape]()])

            elif shape_type in ['stp', 'step']:
                reader = step_reader(str(io.shapes()[shape_name][:]))

                readers[shape_name] = reader
                mapper = vtk.vtkDataSetMapper()
                add_compatiblity_methods(mapper)
                mapper.SetInputConnection(reader.GetOutputPort())
                mappers[shape_name] = (x for x in [mapper])
            else:
                assert shape_type in ['igs', 'iges']

                reader = iges_reader(str(io.shapes()[shape_name][:]))

                readers[shape_name] = reader
                mapper = vtk.vtkDataSetMapper()
                add_compatiblity_methods(mapper)
                mapper.SetInputConnection(reader.GetOutputPort())
                mappers[shape_name] = (x for x in [mapper])

        elif shape_type == 'heightmap':
            points = vtk.vtkPoints()
            shape = io.shapes()[shape_name]
            extents = list(shape.attrs['rect']) + [numpy.max(shape) - numpy.min(shape)]

            # Data points are adjusted to center tangentially, but
            # vertical position is left alone; i.e., supports
            # non-zero-centered data.  User must use contactor
            # translation to compensate if desired, or simply adjust
            # data itself to desired origin.
            for x,d in enumerate(shape):
                for y,v in enumerate(d):
                    points.InsertNextPoint(
                        float(x) / (shape.shape[0]-1) * extents[0] - extents[0]/2,
                        float(y) / (shape.shape[1]-1) * extents[1] - extents[1]/2,
                        v)

            polydata = vtk.vtkPolyData()
            polydata.SetPoints(points)
            delaunay = vtk.vtkDelaunay2D()
            delaunay.SetInputData(polydata)
            delaunay.Update()
            datasets[shape_name] = polydata
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(delaunay.GetOutputPort())
            add_compatiblity_methods(mapper)
            mappers[shape_name] = None
            mappers[shape_name] = (x for x in [mapper])

        elif shape_type == 'convex':
            # a convex shape
            points = vtk.vtkPoints()
            convex = vtk.vtkConvexPointSet()
            data = io.shapes()[shape_name][:]
            convex.GetPointIds().SetNumberOfIds(data.shape[0])
            for id_, vertice in enumerate(io.shapes()[shape_name][:]):
                points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
                convex.GetPointIds().SetId(id_, id_)
            convex_grid = vtk.vtkUnstructuredGrid()
            convex_grid.Allocate(1, 1)
            convex_grid.InsertNextCell(
                convex.GetCellType(), convex.GetPointIds())
            convex_grid.SetPoints(points)

            # not a source!
            datasets[shape_name] = convex_grid

            mapper = vtk.vtkDataSetMapper()
            add_compatiblity_methods(mapper)
            mapper.SetInputData(convex_grid)
            # delayed
            mappers[shape_name] = (x for x in [mapper])

        else:
            assert shape_type == 'primitive'
            primitive = io.shapes()[shape_name].attrs['primitive']
            attrs = io.shapes()[shape_name][:][0]
            if primitive == 'Sphere':
                source = vtk.vtkSphereSource()
                source.SetRadius(attrs[0])
                source.SetThetaResolution(15)
                source.SetPhiResolution(15)

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
            mapper = vtk.vtkCompositePolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            mappers[shape_name] = (x for x in [mapper])

    fixed_mappers = dict()
    for shape_name in io.shapes():
        if shape_name not in fixed_mappers:
            fixed_mappers[shape_name] = mappers[shape_name].next()

    for instance_name in io.instances():

        instance = int(io.instances()[instance_name].attrs['id'])
        contactors[instance] = []
        transforms[instance] = []
        offsets[instance] = []
        times_of_birth[instance] = io.instances()[instance_name].\
                                   attrs.get('time_of_birth',-1)

        if instance >= 0:
            dynamic_actors[instance] = list()
        else:
            static_actors[instance] = list()

        for contactor_instance_name in io.instances()[instance_name]:
            contactor_shape_name = io.instances()[instance_name][
                contactor_instance_name].attrs['name']
            contactors[instance].append(contactor_shape_name)

            actor = vtk.vtkActor()
            if io.instances()[instance_name].attrs.get('mass', 0) > 0:
                # objects that may move
                dynamic_actors[instance].append(actor)

                actor.GetProperty().SetOpacity(
                    config.get('dynamic_opacity', 0.7))
            else:
                # objects that are not supposed to move
                static_actors[instance].append(actor)

                actor.GetProperty().SetOpacity(
                    config.get('static_opacity', 1.0))

            actor.GetProperty().SetColor(random_color())
            actor.SetMapper(fixed_mappers[contactor_shape_name])

            renderer.AddActor(actor)

            transform = vtk.vtkTransform()
            transformer = vtk.vtkTransformFilter()

            if contactor_shape_name in readers:
                transformer.SetInputConnection(
                    readers[contactor_shape_name].GetOutputPort())
            else:
                transformer.SetInputData(datasets[contactor_shape_name])

            if 'scale' in io.shapes()[contactor_shape_name].attrs:
                scale = io.shapes()[contactor_shape_name].attrs['scale']
                scale_transform = vtk.vtkTransform()
                scale_transform.Scale(scale, scale, scale)
                scale_transform.SetInput(transform)
                transformer.SetTransform(scale_transform)
                actor.SetUserTransform(scale_transform)
            else:
                transformer.SetTransform(transform)
                actor.SetUserTransform(transform)

            transformers[contactor_shape_name] = transformer

            big_data_source.AddInputConnection(transformer.GetOutputPort())

            transforms[instance].append(transform)

            if 'center_of_mass' in io.instances()[instance_name].attrs:
                center_of_mass = io.instances()[instance_name].\
                                 attrs['center_of_mass'].astype(float)
            else:
                center_of_mass = [0., 0., 0.]

            offsets[instance].append(
                (numpy.subtract(io.instances()[
                    instance_name][
                        contactor_instance_name].attrs['translation'].astype(float),
                                center_of_mass),
                    io.instances()[instance_name][contactor_instance_name].\
                 attrs['orientation'].astype(float)))

    pos_data = dpos_data[:]
    spos_data = spos_data[:]

    # this sets the position for all transforms associated to an instance
    def set_position_i(instance, q0, q1, q2, q3, q4, q5, q6):
        #if (numpy.any(numpy.isnan([q0, q1, q2, q3, q4, q5, q6]))
        #    or numpy.any(numpy.isinf([q0, q1, q2, q3, q4, q5, q6]))):
        #print('Bad position for', instance, q0, q1, q2, q3, q4, q5, q6)
        #return

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

    # the numpy vectorization is ok on column vectors for each args
    set_position_v = numpy.vectorize(set_position_i)

    def set_position(data):
        set_position_v(data[:, 1],
                       data[:, 2],
                       data[:, 3],
                       data[:, 4],
                       data[:, 5],
                       data[:, 6],
                       data[:, 7],
                       data[:, 8])

    # to be removed if ok
    def set_positionv_old(x):
        for y in x.reshape(-1, 8):
            set_position(*(y.reshape(8)))

    # set visibility for all actors associated to a dynamic instance
    def set_actors_viz(instance, time):
        if times_of_birth[instance] <= time:
            for actor in dynamic_actors[instance]:
                actor.VisibilityOn()
        else:
            for actor in dynamic_actors[instance]:
                actor.VisibilityOff()

    # here the numpy vectorization is used with a column vector and a
    # scalar for the time arg
    set_actors_vizzz = numpy.vectorize(set_actors_viz)

    def set_dynamic_actors_visibility(time):
        set_actors_vizzz(dynamic_actors.keys(), time)


    # to be removed if ok
    def set_dynamic_actors_visibility_old(id_t=None):
        for instance, actor in dynamic.actors.items() + static_actors.items():
            # Instance is a static object
            visible = instance < 0
            # Instance is in the current timestep
            if id_t:
                visible = visible or instance in pos_data[id_t, 1]
            # Instance has time of birth <= 0
            else:
                tob = [io.instances()[k].attrs['time_of_birth']
                       for k in io.instances()
                       if io.instances()[k].attrs['id'] == instance][0]
                visible = visible or tob <= 0
            if visible:
                [actor.VisibilityOn() for actor in actorlist]
            else:
                [actor.VisibilityOff() for actor in actorlist]

    if cf_prov is not None:
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

    time0 = None
    try:
        # Positions at first time step
        time0 = min(dpos_data[:, 0])
        id_t0 = numpy.where(dpos_data[:, 0] == time0)
        pos_t0 = pos_data[id_t0, 0:9]
    except ValueError:
        # this is for the case simulation hass not been ran and
        # time does not exists
        time0 = 0
        id_t0 = None
        pos_t0 = numpy.array([
            numpy.hstack(([0.,
                           io.instances()[k].attrs['id']]
                          ,io.instances()[k].attrs['translation']
                          ,io.instances()[k].attrs['orientation']))
            for k in io.instances()
            if io.instances()[k].attrs['id'] >= 0])

    if numpy.shape(spos_data)[0] > 0:
        set_position(spos_data)
        print(spos_data.shape)
        # static objects are always visible
        for instance, actors in static_actors.items():
            for actor in actors:
                 actor.VisibilityOn()

    set_position(*pos_t0)

    set_dynamic_actors_visibility(time0)

    renderer_window.AddRenderer(renderer)
    interactor_renderer.SetRenderWindow(renderer_window)
    interactor_renderer.GetInteractorStyle(
    ).SetCurrentStyleToTrackballCamera()

    # http://www.itk.org/Wiki/VTK/Depth_Peeling ?

    # Use a render window with alpha bits (as initial value is 0 (false) ):
    renderer_window.SetAlphaBitPlanes(1)

    # Force to not pick a framebuffer with a multisample buffer ( as initial
    # value is 8):
    renderer_window.SetMultiSamples(0)

    # Choose to use depth peeling (if supported) (initial value is 0
    # (false) )
    renderer.SetUseDepthPeeling(1)

    # Set depth peeling parameters.
    renderer.SetMaximumNumberOfPeels(100)

    # Set the occlusion ratio (initial value is 0.0, exact image)
    renderer.SetOcclusionRatio(0.1)

    # Set the initial camera position and orientation if specified
    if initial_camera[0] is not None:
        renderer.GetActiveCamera().SetPosition(*initial_camera[0])
    if initial_camera[1] is not None:
        renderer.GetActiveCamera().SetFocalPoint(*initial_camera[1])
    if initial_camera[2] is not None:
        renderer.GetActiveCamera().SetViewUp(*initial_camera[2])
    if initial_camera[3] is not None:
        renderer.GetActiveCamera().ParallelProjectionOn()
        renderer.GetActiveCamera().SetParallelScale(initial_camera[3])

    # callback maker for scale manipulation
    def make_scale_observer(glyphs):

        def scale_observer(obj, event):
            slider_repres = obj.GetRepresentation()
            scale_at_pos = slider_repres.GetValue()
            for glyph in glyphs:
                for k in glyph:
                    glyph[k].SetScaleFactor(
                        scale_at_pos * glyph[k]._scale_fact)

        return scale_observer

    # callback maker for time scale manipulation
    def make_time_scale_observer(time_slider_repres, time_observer):

        delta_time = max_time - min_time

        def time_scale_observer(obj, event):
            slider_repres = obj.GetRepresentation()
            time_scale_at_pos = 1. - slider_repres.GetValue()

            current_time = time_observer._time

            shift = (current_time - min_time) / delta_time

            xmin_time = min_time + time_scale_at_pos / 2. * delta_time
            xmax_time = max_time - time_scale_at_pos / 2. * delta_time

            xdelta_time = xmax_time - xmin_time

            new_mintime = max(min_time, current_time - xdelta_time)
            new_maxtime = min(max_time, current_time + xdelta_time)

            time_slider_repres.SetMinimumValue(new_mintime)
            time_slider_repres.SetMaximumValue(new_maxtime)

        return time_scale_observer

    # make a slider widget and its representation
    def make_slider(title, observer, interactor,
                    startvalue, minvalue, maxvalue, cx1, cy1, cx2, cy2):
        slider_repres = vtk.vtkSliderRepresentation2D()
        slider_repres.SetMinimumValue(
            minvalue - (maxvalue - minvalue) / 100)
        slider_repres.SetMaximumValue(
            maxvalue + (maxvalue - minvalue) / 100)
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

    recorder = vtk.vtkOggTheoraWriter()
    recorder.SetQuality(2)
    recorder.SetRate(frames_per_second)
    recorder.SetFileName(os.path.splitext(io_filename)[0]+'.avi')
    recorder.SetInputConnection(image_maker.GetOutputPort())

    writer = vtk.vtkPNGWriter()
    writer.SetInputConnection(image_maker.GetOutputPort())

    class InputObserver():

        def __init__(self, times=None, slider_repres=None):
            self._opacity = 1.0
            self._current_id = vtk.vtkIdTypeArray()
            self._renderer = renderer
            self._renderer_window = renderer_window
            self._image_counter = 0
            self._view_cycle = -1
            self._recording = False
            self._times = None

            if times is None or len(times)==0:
                return
            self._times = times
            self._stimes = set(times)
            self._time_step = (max(self._stimes) - min(self._stimes)) \
                / len(self._stimes)
            self._time = min(times)
            if slider_repres is None:
                return
            self._slider_repres = slider_repres

        def update(self):
            global cf_prov
            if self._times is None:
                renderer_window.Render()
                return
            index = bisect.bisect_left(self._times, self._time)
            index = max(0, index)
            index = min(index, len(self._times) - 1)
            if cf_prov is not None:
                cf_prov._time = self._times[index]
                cf_prov.xmethod()

    #        contact_posa.Update()

    #        contact_posb.Update()

            # contact_pos_force.Update()
            # arrow_glyph.Update()
            # gmapper.Update()

            # set_positionv(spos_data[:, 1:9])

            set_dynamic_actors_visibility(self._times[index])

            id_t = numpy.where(pos_data[:, 0] == self._times[index])
            set_position(*pos_data[id_t, :])

            self._slider_repres.SetValue(self._time)

            self._current_id.SetNumberOfValues(1)
            self._current_id.SetValue(0, index)

            self._iter_plot.SetSelection(self._current_id)
            self._prec_plot.SetSelection(self._current_id)

            renderer_window.Render()

        def object_pos(self, id_):
            index = bisect.bisect_left(self._times, self._time)
            index = max(0, index)
            index = min(index, len(self._times) - 1)

            id_t = numpy.where(pos_data[:, 0] == self._times[index])
            return (pos_data[id_t[0][id_], 2], pos_data[id_t[0][id_], 3], pos_data[id_t[0][id_], 4])

        def set_opacity(self):
            for instance, actors in dynamic_actors.items():
                for actor in actors:
                    actor.GetProperty().SetOpacity(self._opacity)

        def key(self, obj, event):
            global cf_prov
            key = obj.GetKeySym()
            print('key', key)

            if key == 'r':
                spos_data, dpos_data, dom_data, cf_data, solv_data = load()
                if not cf_disable:
                    cf_prov = CFprov(cf_data, dom_data)
                times = list(set(dpos_data[:, 0]))
                times.sort()

                if len(spos_data) > 0:
                    instances = set(dpos_data[:, 1]).union(
                        set(spos_data[:, 1]))
                else:
                    instances = set(dpos_data[:, 1])

                if cf_prov is not None:
                    cf_prov._time = min(times[:])
                    cf_prov.xmethod()
                    for mu in cf_prov._mu_coefs:
                        contact_posa[mu].SetInputData(cf_prov._output)
                        contact_posa[mu].Update()
                        contact_posb[mu].SetInputData(cf_prov._output)
                        contact_posb[mu].Update()
                        contact_pos_force[mu].Update()
                        contact_pos_norm[mu].Update()

                id_t0 = numpy.where(
                    dpos_data[:, 0] == time0)

                pos_data = dpos_data[:]
                min_time = times[0]
                set_dynamic_actors_visibility(time0)

                max_time = times[len(times) - 1]

                self._slider_repres.SetMinimumValue(min_time)
                self._slider_repres.SetMaximumValue(max_time)
                self.update()

            if key == 'p':
                self._image_counter += 1
                image_maker.Update()
                writer.SetFileName(
                    'vview-{0}.png'.format(self._image_counter))
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
                    print('camera position:', self._renderer.GetActiveCamera().GetPosition())
                    print('camera focal point', self._renderer.GetActiveCamera().GetFocalPoint())
                    print('camera clipping plane', self._renderer.GetActiveCamera().GetClippingRange())
                    print('camera up vector', self._renderer.GetActiveCamera().GetViewUp())
                    if self._renderer.GetActiveCamera().GetParallelProjection() != 0:
                        print('camera parallel scale', self._renderer.GetActiveCamera().GetParallelScale())

            if key == 'o':
                    self._renderer.GetActiveCamera().SetParallelProjection(
                        1 - self._renderer.GetActiveCamera().GetParallelProjection())

            if key == 'v':
                    # Cycle through some useful views
                    dist = norm(self._renderer.GetActiveCamera().GetPosition())
                    # dist2 = norm([numpy.sqrt(dist**2)/3]*2)
                    d3 = norm([numpy.sqrt(dist**2) / 3] * 3)
                    self._view_cycle += 1
                    if self._view_cycle == 0:
                            print('Left')
                            self._renderer.GetActiveCamera().SetPosition(
                                dist, 0, 0)
                            self._renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
                            self._renderer.GetActiveCamera().SetViewUp(0, 0, 1)
                    elif self._view_cycle == 1:
                            print('Right')
                            self._renderer.GetActiveCamera().SetPosition(
                                0, dist, 0)
                            self._renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
                            self._renderer.GetActiveCamera().SetViewUp(0, 0, 1)
                    elif self._view_cycle == 2:
                            print('Top')
                            self._renderer.GetActiveCamera().SetPosition(
                                0, 0, dist)
                            self._renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
                            self._renderer.GetActiveCamera().SetViewUp(1, 0, 0)
                    else:  # Corner
                            print('Corner')
                            self._renderer.GetActiveCamera().SetPosition(
                                d3, d3, d3)
                            self._renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
                            self._renderer.GetActiveCamera().SetViewUp(
                                -1, -1, 1)
                            self._view_cycle = -1
                    self._renderer.ResetCameraClippingRange()

            if key == 's':
                if not self._recording:
                    recorder.Start()
                    self._recording = True

            if key == 'e':
                if self._recording:
                    self._recording = False
                    recorder.End()

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
                    self._time = self._times[
                        self._iter_plot.GetSelection().GetValue(0)]
                    # -> recompute index ...
                    self.update()

        def prec_plot_observer(self, obj, event):
            if self._prec_plot.GetSelection() is not None:
                # just one selection at the moment!
                if self._prec_plot.GetSelection().GetMaxId() >= 0:
                    self._time = self._times[
                        self._prec_plot.GetSelection().GetValue(0)]
                    # -> recompute index ...
                    self.update()

        def recorder_observer(self, obj, event):
            if self._recording:
                if advance_by_time is not None:
                    self._time += advance_by_time
                    slwsc.SetEnabled(False)  # Scale slider
                    xslwsc.SetEnabled(False)  # Time scale slider
                    # slider_widget.SetEnabled(False) # Time slider
                    # widget.SetEnabled(False) # Axis widget
                self.update()
                image_maker.Modified()
                recorder.Write()
                if advance_by_time is not None:
                    slwsc.SetEnabled(True)
                    xslwsc.SetEnabled(True)
                    # slider_widget.SetEnabled(True)
                    # widget.SetEnabled(True) # Axis widget

                if self._time >= max(self._times):
                    self._recording = False
                    recorder.End()

    if len(times) > 0:
        slider_repres = vtk.vtkSliderRepresentation2D()

        if min_time is None:
            min_time = times[0]
        if max_time is None:
            max_time = times[len(times) - 1]

        slider_repres.SetMinimumValue(min_time)
        slider_repres.SetMaximumValue(max_time)
        slider_repres.SetValue(min_time)
        slider_repres.SetTitleText("time")
        slider_repres.GetPoint1Coordinate(
        ).SetCoordinateSystemToNormalizedDisplay()
        slider_repres.GetPoint1Coordinate().SetValue(0.4, 0.9)
        slider_repres.GetPoint2Coordinate(
        ).SetCoordinateSystemToNormalizedDisplay()
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
    else:
        input_observer = InputObserver()
    interactor_renderer.AddObserver('KeyPressEvent', input_observer.key)

    # Create a vtkLight, and set the light parameters.
    light = vtk.vtkLight()
    light.SetFocalPoint(0, 0, 0)
    light.SetPosition(0, 0, 500)
    # light.SetLightTypeToHeadlight()
    renderer.AddLight(light)

    hlight = vtk.vtkLight()
    hlight.SetFocalPoint(0, 0, 0)
    # hlight.SetPosition(0, 0, 500)
    hlight.SetLightTypeToHeadlight()
    renderer.AddLight(hlight)

    # Warning! numpy support offer a view on numpy array
    # the numpy array must not be garbage collected!
    nxtime = solv_data[:, 0]
    nxiters = solv_data[:, 1]
    nprecs = solv_data[:, 2]
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
    # table.Dump()

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

    # screen_size = renderer_window.GetScreenSize()
    renderer_window.SetSize(*config['window_size'])
    renderer_window.SetWindowName('vview: ' + io_filename)
    tview_iter.GetRenderer().GetRenderWindow().SetSize(600, 200)
    tview_prec.GetRenderer().GetRenderWindow().SetSize(600, 200)

    tview_iter.GetInteractor().Initialize()
    # tview_iter.GetInteractor().Start()
    tview_iter.GetRenderer().SetBackground(.9, .9, .9)
    tview_iter.GetRenderer().Render()

    tview_prec.GetInteractor().Initialize()
    # tview_prec.GetInteractor().Start()
    tview_prec.GetRenderer().SetBackground(.9, .9, .9)
    tview_prec.GetRenderer().Render()

    if io.contact_forces_data().shape[0] > 0:
        slwsc, slrepsc = make_slider('Scale',
                                 make_scale_observer([cone_glyph, cylinder_glyph, sphere_glypha, sphere_glyphb, arrow_glyph]
                                                     ),
                                 interactor_renderer,
                                 cf_scale_factor, cf_scale_factor -
                                 cf_scale_factor / 2,
                                 cf_scale_factor + cf_scale_factor / 2,
                                 0.01, 0.01, 0.01, 0.7)

    if len(times) > 0:
        xslwsc, xslrepsc = make_slider('Time scale',
                                   make_time_scale_observer(
                                       slider_repres, input_observer),
                                   interactor_renderer,
                                   time_scale_factor, time_scale_factor -
                                   time_scale_factor / 2,
                                   time_scale_factor + time_scale_factor / 2,
                                   0.1, 0.9, 0.3, 0.9)

    renderer_window.Render()
    interactor_renderer.Initialize()

    # display coordinates axes
    axes = vtk.vtkAxesActor()
    axes.SetTotalLength(1.0, 1.0, 1.0)
    widget = vtk.vtkOrientationMarkerWidget()
    # widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 )
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(interactor_renderer)
    # widget.SetViewport( 0.0, 0.0, 40.0, 40.0 );
    widget.SetEnabled(True)
    widget.InteractiveOn()

    timer_id = interactor_renderer.CreateRepeatingTimer(40)   # 25 fps
    interactor_renderer.AddObserver(
        'TimerEvent', input_observer.recorder_observer)

    interactor_renderer.Start()

## Finalize on quit

# Update configuration and save it
config['window_size'] = renderer_window.GetSize()
if should_save_config:
    save_configuration()
