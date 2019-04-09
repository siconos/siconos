# Mechanics IO

from __future__ import print_function

import os
import sys

from math import cos, sin, asin, atan2
from scipy import constants

import numpy as np
import h5py
import bisect
import time
import pickle

import tempfile
from contextlib import contextmanager

# Siconos Mechanics imports
from siconos.mechanics.collision.tools import Contactor, Volume, Shape

### Constants

joint_points_axes = {
    'KneeJointR': (1, 0),
    'PivotJointR': (1, 1),
    'PrismaticJointR': (0, 1),
    'CylindricalJointR': (1, 1),
    'FixedJointR': (0, 0),
}

### Utility functions

def floatv(v):
    return [float(x) for x in v]


def arguments():
    """Returns tuple containing dictionary of calling function's
    named arguments and a list of calling function's unnamed
    positional arguments.
    """
    from inspect import getargvalues, stack
    posname, kwname, args = getargvalues(stack()[1][0])[-3:]
    posargs = args.pop(posname, [])
    args.update(args.pop(kwname, []))
    return args, posargs


def check_points_axes(name, joint_class, points, axes):
    def check(x, idx):
        def er():
            n = joint_points_axes[joint_class][idx]
            raise ValueError('{} ({}) expects {} {} (got {})'
                             .format(joint_class, name, n,
                                     ['point','points','axis','axes'][idx*2+1*(n!=1)],
                                     x))
        if np.shape(x)==(0,) or np.shape(x)==():
            num = 0
        else:
            if len(np.shape(x))!=2 or np.shape(x)[1]!=3: er()
            num = np.shape(x)[0]
        if (joint_class in joint_points_axes
            and joint_points_axes[joint_class][idx] != num): er()
    check(points, 0)
    check(axes, 1)

@contextmanager
def tmpfile(suffix='', prefix='siconos_io', contents=None,
            debug=False):
    """
    A context manager for a named temporary file.
    """
    (_, tfilename) = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    fid = open(tfilename, 'w')
    if contents is not None:
        fid.write(contents)
        fid.flush()

    class TmpFile:

        def __init__(self, fid, name):
            self.fid = fid
            self.name = name

        def __getitem__(self, n):
            if n == 0:
                return self.fid
            elif n == 1:
                return self.name
            else:
                raise IndexError

    r = TmpFile(fid, tfilename)

    yield r
    fid.close()
    if not debug:
        os.remove(tfilename)


def warn(msg):
    sys.stderr.write('{0}: {1}'.format(sys.argv[0], msg))


def object_id(obj):
    """returns an unique object identifier"""
    return obj.__hash__()


def group(h, name, must_exist=True):
    try:
        return h[name]
    except KeyError:
        if must_exist:
            return h.create_group(name)
        else:
            try:
                return h.create_group(name)
            except ValueError:
                # could not create group, return None
                # (file is probably in read-only mode)
                return None

def data(h, name, nbcolumns, use_compression=False):
    try:
        return h[name]
    except KeyError:
        comp = use_compression and nbcolumns > 0
        return h.create_dataset(name, (0, nbcolumns),
                                maxshape=(None, nbcolumns),
                                chunks=[None,(4000,nbcolumns)][comp],
                                compression=[None,'gzip'][comp],
                                compression_opts=[None,9][comp])


def add_line(dataset, line):
    dataset.resize(dataset.shape[0] + 1, 0)
    dataset[dataset.shape[0] - 1, :] = line



#
# misc fixes
#
# fix ctr.'name' in old hdf5 files
#
def upgrade_io_format(filename):

    with MechanicsHdf5(filename, mode='a') as io:

        for instance_name in io.instances():
            for contactor_instance_name in io.instances()[instance_name]:
                contactor = io.instances()[instance_name][
                    contactor_instance_name]
                if 'name' in contactor.attrs:
                    warn("""
contactor {0} attribute 'name': renamed in 'shape_name'
                    """)
                    contactor.attrs['shape_name'] = contactor['name']
                    del contactor['name']


def str_of_file(filename):
    with open(filename, 'r') as f:
        return str(f.read())


def file_of_str(filename, string):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

    with open(filename, "w") as f:
        f.write(string)


#
# fix orientation -> rotation ?
#
def quaternion_get(orientation):
    """
    Get quaternion from orientation
    """
    if len(orientation) == 2:
            # axis + angle
            axis=orientation[0]
            assert len(axis) == 3
            angle=orientation[1]
            assert type(angle) is float
            n=sin(angle / 2.) / np.linalg.norm(axis)

            ori=[cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
    else:
        assert(len(orientation)==4)
        # a given quaternion
        ori=orientation
    return ori

def quaternion_multiply(q1, q0):
    w0, x0, y0, z0 = q0
    w1, x1, y1, z1 = q1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)



def phi(q0, q1, q2, q3):
    """
    Euler angle phi from quaternion.
    """
    return atan2(2*(q0*q1+q2*q3), 1-2*(q1*q1+q2*q2))


def theta(q0, q1, q2, q3):
    """
    Euler angle theta from quaternion.
    """
    return asin(2*(q0*q2-q3*q1))


def psi(q0, q1, q2, q3):
    """
    Euler angle psi from quaternion.
    """
    return atan2(2*(q0*q3+q1*q2), 1-2*(q2*q2+q3*q3))

# vectorized versions
phiv = np.vectorize(phi)
thetav = np.vectorize(theta)
psiv = np.vectorize(psi)


#
# inertia
#
def compute_inertia_and_center_of_mass(shapes, io=None):
    """
    Compute inertia from a list of Shapes.
    """
    from OCC.GProp import GProp_GProps
    from OCC.BRepGProp import brepgprop_VolumeProperties
    from OCC.gp import gp_Ax1, gp_Dir
    from siconos.mechanics import occ

    system = GProp_GProps()

    for shape in shapes:

        iprops = GProp_GProps()

        if shape.data is None:
            if io is not None:
                shape.data = io._shape.get(shape.shape_name, new_instance=True)
            else:
                warn('cannot get shape {0}'.format(shape.shape_name))
                return None

        iishape = shape.data

        ishape = occ.OccContactShape(iishape).data()
        # the shape relative displacement
        occ.occ_move(ishape, list(shape.translation) +
                     list(shape.orientation))

        brepgprop_VolumeProperties(iishape, iprops)

        density = None

        if hasattr(shape, 'mass') and shape.mass is not None:
            density = shape.mass / iprops.Mass()

        elif shape.parameters is not None and \
             hasattr(shape.parameters, 'density'):
            density = shape.parameters.density
            #print('shape.parameters.density:', shape.parameters.density)
        else:
            density = 1.

        assert density is not None
        # print("shape", shape.shape_name)
        # print('density:', density)
        # print('iprops.Mass():', iprops.Mass())

        system.Add(iprops, density)


    mass=  system.Mass()
    assert (system.Mass() > 0.)

    computed_com = system.CentreOfMass()

    gp_mat= system.MatrixOfInertia()
    inertia_matrix = np.zeros((3,3))
    for i in range(0,3):
        for j in range(0,3):
            inertia_matrix[i,j]=  gp_mat.Value(i+1,j+1)

    I1 =  system.MomentOfInertia(
        gp_Ax1(computed_com, gp_Dir(1, 0, 0)))
    I2 =  system.MomentOfInertia(
        gp_Ax1(computed_com, gp_Dir(0, 1, 0)))
    I3 = system.MomentOfInertia(
        gp_Ax1(computed_com, gp_Dir(0, 0, 1)))

    inertia = [I1, I2, I3]
    center_of_mass = np.array([computed_com.Coord(1),
                               computed_com.Coord(2),
                               computed_com.Coord(3)])

    return mass, center_of_mass, inertia, inertia_matrix


def occ_topo_list(shape):
    """ return the edges & faces from `shape`

    :param shape: a TopoDS_Shape
    :return: a list of edges and faces
    """

    from OCC.TopAbs import TopAbs_FACE
    from OCC.TopAbs import TopAbs_EDGE
    from OCC.TopExp import TopExp_Explorer
    from OCC.TopoDS import topods_Face, topods_Edge


    topExp = TopExp_Explorer()
    topExp.Init(shape, TopAbs_FACE)
    faces = []
    edges = []

    while topExp.More():
        face = topods_Face(topExp.Current())
        faces.append(face)
        topExp.Next()

    topExp.Init(shape, TopAbs_EDGE)

    while topExp.More():
        edge = topods_Edge(topExp.Current())
        edges.append(edge)
        topExp.Next()

    return faces, edges


def occ_load_file(filename):
    """
    load in pythonocc a igs or step file

    :param filename: a filename with extension
    :return: a topods_shape
    """

    from OCC.STEPControl import STEPControl_Reader
    from OCC.IGESControl import IGESControl_Reader
    from OCC.BRep import BRep_Builder
    from OCC.TopoDS import TopoDS_Compound
    from OCC.IFSelect import IFSelect_RetDone,\
    IFSelect_ItemsByEntity

    reader_switch = {'stp': STEPControl_Reader,
                     'step': STEPControl_Reader,
                     'igs': IGESControl_Reader,
                     'iges': IGESControl_Reader}

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)

    reader = reader_switch[os.path.splitext(filename)[1][1:].lower()]()

    status = reader.ReadFile(filename)

    if status == IFSelect_RetDone:  # check status
        failsonly = False
        reader.PrintCheckLoad(
            failsonly, IFSelect_ItemsByEntity)
        reader.PrintCheckTransfer(
            failsonly, IFSelect_ItemsByEntity)

        ok = reader.TransferRoots()
        nbs = reader.NbShapes()

        for i in range(1, nbs + 1):
            shape = reader.Shape(i)
            builder.Add(comp, shape)

    return comp


def topods_shape_reader(shape, deflection=0.001):

    from OCC.StlAPI import StlAPI_Writer
    from OCC.BRepMesh import BRepMesh_IncrementalMesh

    import vtk

    stl_writer = StlAPI_Writer()

    with tmpfile(suffix='.stl') as tmpf:
        mesh = BRepMesh_IncrementalMesh(shape, deflection)
        mesh.Perform()
        assert mesh.IsDone()
        stl_writer.SetASCIIMode(False)
        stl_writer.Write(shape, tmpf[1])
        tmpf[0].flush()

        reader = vtk.vtkSTLReader()
        reader.SetFileName(tmpf[1])
        reader.Update()

        return reader


def brep_reader(brep_string, indx):

    from OCC.StlAPI import StlAPI_Writer
    from OCC.BRepTools import BRepTools_ShapeSet
    import vtk

    shape_set = BRepTools_ShapeSet()
    shape_set.ReadFromString(brep_string)
    shape = shape_set.Shape(shape_set.NbShapes())
    location = shape_set.Locations().Location(indx)
    shape.Location(location)

    stl_writer = StlAPI_Writer()

    with tmpfile(suffix='.stl') as tmpf:
        stl_writer.Write(shape, tmpf[1])
        tmpf[0].flush()

        reader = vtk.vtkSTLReader()
        reader.SetFileName(tmpf[1])
        reader.Update()

        return reader


class MechanicsHdf5(object):

    """a MechanicsHdf5 context manager prepares a simulation description
    to be executed by MechanicsRunner.
    """

    def __init__(self, io_filename=None, mode='w',
                 use_compression=False, output_domains=False, verbose=True):
        if io_filename is None:
            self._io_filename = '{0}.hdf5'.format(
                os.path.splitext(os.path.basename(sys.argv[0]))[0])
        else:
            self._io_filename = io_filename
        self._mode = mode
        self._static_data = None
        self._velocities_data = None
        self._dynamic_data = None
        self._cf_data = None
        self._domain_data = None
        self._solv_data = None
        self._input = None
        self._nslaws_data = None
        self._nslaws = dict()
        self._out = None
        self._data = None
        self._ref = None
        self._permanent_interactions = None
        self._joints = None
        self._boundary_conditions = None
        self._plugins = None
        self._external_functions = None
        self._number_of_shapes = 0
        self._number_of_permanent_interactions = 0
        self._number_of_dynamic_objects = 0
        self._number_of_static_objects = 0
        self._use_compression = use_compression
        self._should_output_domains = output_domains
        self._verbose = verbose

    def __enter__(self):
        self._out = h5py.File(self._io_filename, self._mode)
        self._data = group(self._out, 'data')
        self._ref = group(self._data, 'ref')
        self._permanent_interactions = group(self._data, 'permanent_interactions',
                                             must_exist=False)
        self._joints = group(self._data, 'joints', must_exist=False)
        self._plugins = group(self._data, 'plugins', must_exist=False)
        self._external_functions = group(self._data, 'external_functions', must_exist=False)
        try:
            self._boundary_conditions = group(self._data, 'boundary_conditions',
                                              must_exist=(self._mode=='w'))
        except Exception as e :
            print('Warning -  group(self._data, boundary_conditions ) : ',  e)
        self._static_data = data(self._data, 'static', 9,
                                 use_compression = self._use_compression)
        self._velocities_data = data(self._data, 'velocities', 8,
                                     use_compression = self._use_compression)
        self._dynamic_data = data(self._data, 'dynamic', 9,
                                  use_compression = self._use_compression)
        self._cf_data = data(self._data, 'cf', 26,
                             use_compression = self._use_compression)
        if self._should_output_domains or 'domain' in self._data:
            self._domain_data = data(self._data, 'domain', 3,
                                     use_compression = self._use_compression)
        self._solv_data = data(self._data, 'solv', 4,
                               use_compression = self._use_compression)
        self._input = group(self._data, 'input')

        self._nslaws_data = group(self._data, 'nslaws')
        return self

    def __exit__(self, type_, value, traceback):
        self._out.close()

    def print_verbose(self, *args, **kwargs):
            if self._verbose:
                print('[io.mechanics]', *args, **kwargs)


# hdf5 structure

    def shapes(self):
        """
        Shapes : parameterized primitives or user defined
                 (convex set or meshes)
        """
        return self._ref

    def permanent_interactions(self):
        """
        Permanent interactions.
        """
        return self._permanent_interactions

    def static_data(self):
        """
        Coordinates and orientations of static objects.
        """
        return self._static_data

    def dynamic_data(self):
        """
        Coordinates and orientations of dynamic objects.
        """
        return self._dynamic_data

    def velocities_data(self):
        """
        Velocities of dynamic objects
        """
        return self._velocities_data

    def contact_forces_data(self):
        """
        Contact points information.
        """
        return self._cf_data

    def domains_data(self):
        """
        Contact point domain information.
        """
        return self._domain_data

    def solver_data(self):
        """
        Solver output
        """
        return self._solv_data

    def instances(self):
        """
        Scene objects.
        """
        return self._input

    def nonsmooth_laws(self):
        """
        Non smooth laws between group of contactors.
        """
        return self._nslaws_data

    def joints(self):
        """
        Joints between dynamic objects or between an object and the scenery.
        """
        return self._joints

    def boundary_conditions(self):
        """
        Boundary conditions applied to  dynamic objects
        """
        return self._boundary_conditions

    def add_plugin_source(self, name, filename):
        """
        Add C source plugin
        """

        if name not in self._plugins:
            plugin_src = self._plugins.create_dataset(name, (1,),
                                                      dtype=h5py.new_vlen(str))
            plugin_src[:] = str_of_file(filename)
            plugin_src.attrs['filename'] = filename

    def add_external_function(self, name, body_name, function_name,
                              plugin_name, plugin_function_name):

        if name not in self._external_functions:
            ext_fun = group(self._external_functions, name)
            ext_fun.attrs['body_name'] = body_name
            ext_fun.attrs['function_name'] = function_name
            ext_fun.attrs['plugin_name'] = plugin_name
            ext_fun.attrs['plugin_function_name'] = plugin_function_name

    def add_external_bc_function(self, name, body_name, bc_indices,
                              plugin_name, plugin_function_name):

        if name not in self._external_functions:
            ext_fun = group(self._external_functions, name)
            ext_fun.attrs['body_name'] = body_name
            ext_fun.attrs['plugin_name'] = plugin_name
            ext_fun.attrs['plugin_function_name'] = plugin_function_name
            ext_fun.attrs['bc_indices'] = bc_indices

    def add_mesh_from_string(self, name, shape_data, scale=None,
                             insideMargin=None, outsideMargin=None):
        """
        Add a mesh shape from a string.
        Accepted format : mesh encoded in VTK .vtp format
        """

        if name not in self._ref:

            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape[:] = shape_data
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'vtp'
            if scale is not None:
                shape.attrs['scale'] = scale
            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            self._number_of_shapes += 1

    def add_mesh_from_file(self, name, filename, scale=None,
                           insideMargin=None, outsideMargin=None):
        """
        Add a mesh shape from a file.
        Accepted format : .stl or mesh encoded in VTK .vtp format
        """

        import vtk

        if filename[0] != os.path.sep:
            filename = os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0],
                                    filename)
        if name not in self._ref:

            if os.path.splitext(filename)[-1][1:] == 'stl':
                reader = vtk.vtkSTLReader()
                reader.SetFileName(filename)
                reader.Update()

                if reader.GetErrorCode() != 0:
                    print('vtkSTLReader error', reader.GetErrorCode())
                    sys.exit(1)

                with tmpfile() as tmpf:
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetInputData(reader.GetOutput())
                    writer.SetFileName(tmpf[1])
                    writer.Write()

                    shape_data = str_of_file(tmpf[1])

            else:
                assert os.path.splitext(filename)[-1][1:] == 'vtp'
                shape_data = str_of_file(filename)

            self.add_mesh_from_string(name, shape_data, scale=scale,
                                   insideMargin=insideMargin,
                                   outsideMargin=outsideMargin)

    def add_height_map(self, name, heightmap, rectangle,
                       insideMargin=None, outsideMargin=None):
        """
        Add a heightmap represented as a SiconosMatrix
        """
        assert(heightmap.shape[0] >= 2)
        assert(heightmap.shape[1] >= 2)
        if name not in self._ref:
            shape = self._ref.create_dataset(name, data=heightmap)
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'heightmap'

            # measurements of the heightfield, i.e. length of sides of
            # the rectangle where heightmap will be placed -- height
            # is represented by heightmap values
            assert(len(rectangle)==2)
            shape.attrs['rect'] = rectangle # tuple (length x, length y)

            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            self._number_of_shapes += 1

    def add_brep_from_string(self, name, shape_data):
        """
        Add a brep contained in a string.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            if type(shape_data) == str:
                # raw str
                shape[:] = shape_data
            else:
                # __getstate__ as with pythonocc
                shape[:] = shape_data[0]
                shape.attrs['occ_indx'] = shape_data[1]

            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'brep'

            self._number_of_shapes += 1

    def add_occ_shape(self, name, occ_shape):
        """
        Add an OpenCascade TopoDS_Shape.
        """

        if name not in self._ref:

            from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs

            # step format is used for the storage.
            step_writer = STEPControl_Writer()

            step_writer.Transfer(occ_shape, STEPControl_AsIs)

            shape_data = None

            with tmpfile() as tmpf:

                status = step_writer.Write(tmpf[1])

                tmpf[0].flush()
                shape_data = str_of_file(tmpf[1])

                shape = self._ref.create_dataset(name, (1,),
                                                 dtype=h5py.new_vlen(str))
                shape[:] = shape_data
                shape.attrs['id'] = self._number_of_shapes
                shape.attrs['type'] = 'step'
                self._number_of_shapes += 1

    def add_shape_data_from_file(self, name, filename):
        """
        Add shape data from a file.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape[:] = str_of_file(filename)
            shape.attrs['id'] = self._number_of_shapes
            try:
                shape.attrs['type'] = os.path.splitext(filename)[1][1:]
            except:
                shape.attrs['type'] = 'unknown'

            self._number_of_shapes += 1

    def add_interaction(self, name, body1_name, contactor1_name=None,
                        body2_name=None, contactor2_name=None,
                        distance_calculator='cadmbtb',
                        offset1=0.0, offset2=0.0):
        """
        Add permanent interactions between two objects contactors.
        """
        if name not in self.permanent_interactions():
            pinter = self.permanent_interactions().\
                      create_dataset(name, (1,),
                                     dtype=h5py.new_vlen(str))
            pinter.attrs['id'] = self._number_of_permanent_interactions
            pinter.attrs['type'] = 'permanent_interaction'
            pinter.attrs['body1_name'] = body1_name
            pinter.attrs['body2_name'] = body2_name
            if contactor1_name is not None:
                pinter.attrs['contactor1_name'] = contactor1_name
            if contactor2_name is not None:
                pinter.attrs['contactor2_name'] = contactor2_name
            pinter.attrs['distance_calculator'] = distance_calculator
            pinter.attrs['offset1'] = offset1
            pinter.attrs['offset2'] = offset2

            self._number_of_permanent_interactions += 1

    def add_convex_shape(self, name, points,
                         insideMargin=None, outsideMargin=None):
        """
        Add a convex shape defined by a list of points.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name,
                                             (np.shape(points)[0],
                                              np.shape(points)[1]))
            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            shape[:] = points[:]
            shape.attrs['type'] = 'convex'
            shape.attrs['id'] = self._number_of_shapes
            self._number_of_shapes += 1

    def add_primitive_shape(self, name, primitive, params,
                            insideMargin=None, outsideMargin=None):
        """
        Add a primitive shape.
        """
        if name not in self._ref:
            shape=self._ref.create_dataset(name, (1, len(params)))
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'primitive'
            shape.attrs['primitive'] = primitive
            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            shape[:] = params
            self._number_of_shapes += 1

    def add_object(self, name, shapes,
                   translation =[0, 0, 0],
                   orientation=[1, 0, 0, 0],
                   velocity=[0, 0, 0, 0, 0, 0],
                   use_volume_centroid_as_initial_translation=False,
                   mass=None, center_of_mass=[0, 0, 0],
                   inertia=None, time_of_birth=-1, time_of_death=-1,
                   allow_self_collide=False):
        """Add an object with associated shapes as a list of Volume or
        Contactor objects. Contact detection and processing is
        defined by the Contactor objects. The Volume objects are used for
        the computation of inertia and center of mass if not provided.

        Each Contactor and Volume object may have a relative
        translation and a relative orientation expressed in the bodyframe
        coordinates.

        Parameters
        ----------
        name: string
            The name of the object.

        shapes: iterable
            The list of associated Contactor or Volume objects.

        translation: array_like of length 3
            Initial translation of the object (mandatory)

        velocity: array_like of length 6
            Initial velocity of the object.
            The components are those of the translation velocity along
            x, y and z axis and the rotation velocity around x, y and
            z axis.  The default velocity is [0, 0, 0, 0, 0, 0].

        mass: float
            The mass of the object, if it is zero the object is defined as
            a static object involved only in contact detection.
            The default value is zero.

        center_of_mass: array_like of length 3
            The position of the center of mass expressed in the body frame
            coordinates.

        inertia: array_like of length 3 or 3x3 matrix.
            The principal moments of inertia (array of length 3) or
            a full 3x3 inertia matrix

        use_volume_centroid_as_initial_translation: boolean.
            if True and if a Volume is given is the list of shape, the position of
            the volume centroid is used as initial translation.
        """
        # print(arguments())
        ori = quaternion_get(orientation)

        assert (len(translation)==3)
        assert (len(ori)==4)
        is_center_of_mass_computed=False
        if name not in self._input:

            com_translation = [0., 0., 0.]

            if (inertia is None) or (mass is None):
                if any(map(lambda s: isinstance(s,Volume), shapes)):

                    # a computed inertia and center of mass
                    # occ only
                    volumes = filter(lambda s: isinstance(s, Volume),
                                     shapes)

                    computed_mass, com, computed_inertia, computed_inertia_matrix = compute_inertia_and_center_of_mass(volumes, self)
                    print('{0}: computed mass from Volume'.format(name))
                    print('{0}: computed center of mass:'.format(name),
                          com[0],
                          com[1],
                          com[2])
                    print('{0}: computed mass:'.format(name),
                          computed_mass)
                    print('{0}: computed inertia:'.format(name),
                          computed_inertia[0], computed_inertia[1], computed_inertia[2])
                    print('{0}: computed inertia matrix:'.format(name),
                          computed_inertia_matrix)
                    is_center_of_mass_computed=True
                    if mass is None:
                        mass=computed_mass

                    if inertia is None:
                        inertia=computed_inertia_matrix

            obj =group(self._input, name)


            if use_volume_centroid_as_initial_translation and is_center_of_mass_computed:
                translation=com
                for s in shapes:
                    s.translation =  s.translation - com





            if time_of_birth >= 0:
                obj.attrs['time_of_birth']=time_of_birth
            if time_of_death >= 0:
                obj.attrs['time_of_death']=time_of_death

            if mass is not None: obj.attrs['mass']=mass
            obj.attrs['translation']=translation
            obj.attrs['orientation']=ori
            obj.attrs['velocity']=velocity
            obj.attrs['center_of_mass']=center_of_mass

            if inertia is not None:
                obj.attrs['inertia']=inertia
            if allow_self_collide is not None:
                obj.attrs['allow_self_collide']=allow_self_collide

            contactors = shapes

            for num, ctor in enumerate(contactors):

                if ctor.instance_name is not None:
                    # a specified name
                    instance_name = ctor.instance_name
                else:
                    # the default name for contactor
                    instance_name = '{0}-{1}'.format(ctor.shape_name, num)

                dat = data(obj, instance_name, 0,
                           use_compression=self._use_compression)

                dat.attrs['instance_name'] = instance_name
                dat.attrs['shape_name'] = ctor.shape_name
                if hasattr(ctor, 'group'):
                    dat.attrs['group'] = ctor.group

                if hasattr(ctor, 'parameters') and \
                        ctor.parameters is not None:
                    ## we add np.void to manage writing string in hdf5 files see http://docs.h5py.org/en/latest/strings.html
                    dat.attrs['parameters'] = np.void(pickle.dumps(ctor.parameters))

                if hasattr(ctor, 'contact_type') and \
                   ctor.contact_type is not None:
                    dat.attrs['type'] = ctor.contact_type

                if hasattr(ctor, 'contact_index') and \
                   ctor.contact_index is not None:
                    dat.attrs['contact_index'] = ctor.contact_index

                dat.attrs['translation'] = ctor.translation
                dat.attrs['orientation'] = quaternion_get(ctor.orientation)

            if mass is None or mass == 0:
                obj.attrs['id']=- (self._number_of_static_objects + 1)
                self._number_of_static_objects += 1

            else:
                obj.attrs['id']=(self._number_of_dynamic_objects + 1)
                self._number_of_dynamic_objects += 1

            return obj

    def add_Newton_impact_rolling_friction_nsl(self, name, mu, mu_r, e=0, collision_group1=0,
                                   collision_group2=0):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactFrictionNSL are supported.
        name is an user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiants.

        """
        if name not in self._nslaws_data:
            nslaw=self._nslaws_data.create_dataset(name, (0,))
            nslaw.attrs['type']='NewtonImpactRollingFrictionNSL'
            nslaw.attrs['mu']=mu
            nslaw.attrs['mu_r']=mu_r
            nslaw.attrs['e']=e
            nslaw.attrs['gid1']=collision_group1
            nslaw.attrs['gid2']=collision_group2

    def add_Newton_impact_friction_nsl(self, name, mu, e=0, collision_group1=0,
                                   collision_group2=0):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactFrictionNSL are supported.
        name is an user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiants.

        """
        if name not in self._nslaws_data:
            nslaw=self._nslaws_data.create_dataset(name, (0,))
            nslaw.attrs['type']='NewtonImpactFrictionNSL'
            nslaw.attrs['mu']=mu
            nslaw.attrs['e']=e
            nslaw.attrs['gid1']=collision_group1
            nslaw.attrs['gid2']=collision_group2

    # Note, default groups are -1 here, indicating not to add them to
    # the nslaw lookup table for contacts, since 1D impacts are
    # useless in this case.  They are however useful for joint stops.
    def add_Newton_impact_nsl(self, name, e=0, collision_group1=-1,
                              collision_group2=-1):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactNSL are supported.
        name is a user identifier and must be unique,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiers.

        As opposed to add_Newton_impact_friction_nsl, the default groups are
        -1, making the NSL unassociated with point contacts.  It can
        by used for joint stops however.
        """
        if name not in self._nslaws_data:
            nslaw=self._nslaws_data.create_dataset(name, (0,))
            nslaw.attrs['type']='NewtonImpactNSL'
            nslaw.attrs['e']=e
            nslaw.attrs['gid1']=collision_group1
            nslaw.attrs['gid2']=collision_group2

    # Note, default groups are -1 here, indicating not to add them to
    # the nslaw lookup table for contacts, since 1D impacts are
    # useless in this case.  They are however useful for joint friction.
    def add_relay_nsl(self, name, lb, ub, size=1, collision_group1=-1,
                      collision_group2=-1):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactNSL are supported.
        name is a user identifier and must be unique,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiers.

        As opposed to add_Newton_impact_friction_nsl, the default groups are
        -1, making the NSL unassociated with point contacts.  It can
        by used for joint stops however.
        """
        if name not in self._nslaws_data:
            nslaw=self._nslaws_data.create_dataset(name, (0,))
            nslaw.attrs['type']='RelayNSL'
            nslaw.attrs['size']=size
            nslaw.attrs['lb']=lb
            nslaw.attrs['ub']=ub
            nslaw.attrs['gid1']=collision_group1
            nslaw.attrs['gid2']=collision_group2

    def add_joint(self, name, object1, object2=None,
                 points=[[0, 0, 0]], axes=[[0, 1, 0]],
                 joint_class='PivotJointR', absolute=None,
                 allow_self_collide=None, nslaws=None, stops=None,
                 friction=None,coupled=None,references=None):
        """
        add a joint between two objects
        """
        if name in self.joints():
            raise ValueError('Joint {} already in simulation!'.format(name))
        else:
            joint=self.joints().create_dataset(name, (0,))
            joint.attrs['object1']=object1
            if object2 is not None:
                joint.attrs['object2']=object2
            joint.attrs['type']=joint_class
            check_points_axes(name, joint_class, points, axes)
            if points is not None:
                joint.attrs['points']=points
            if axes is not None:
                joint.attrs['axes']=axes
            if absolute in [True, False]:
                joint.attrs['absolute']=absolute

            if allow_self_collide in [True, False]:
                joint.attrs['allow_self_collide']=allow_self_collide
            if nslaws is not None:
                # either name of one nslaw, or a list of names same length as stops
                joint.attrs['nslaws'] = np.array(nslaws, dtype='S')
            if stops is not None:
                joint.attrs['stops'] = stops # must be a table of [[axis,pos,dir]..]
            if friction is not None:
                # must be an NSL name (e.g.  RelayNSL), or list of same
                joint.attrs['friction'] = np.array(friction, dtype='S')
            if coupled is not None:
                # must be a list of tuples of two integers (DoF
                # indexes) and a float (ratio)
                for c in coupled:
                    assert(len(c)==3)
                joint.attrs['coupled'] = np.array(coupled)
            if references is not None:
                # must be a list of two joint names and one DS name
                assert(len(references)==2 or len(references)==3)
                joint.attrs['references'] = np.array(references, dtype='S')

    def add_boundary_condition(self, name, object1, indices=None, bc_class='HarmonicBC',
                               v=None, a=None, b=None, omega=None, phi=None):
        """
        add boundarycondition to the object object1

        implementation only works for HarmonicBC for the moment
        """
        if name not in self.boundary_conditions():
            boundary_condition=self.boundary_conditions().create_dataset(name, (0,))
            boundary_condition.attrs['object1']=object1
            boundary_condition.attrs['indices']=indices
            boundary_condition.attrs['type']=bc_class
            if bc_class == 'HarmonicBC' :
                boundary_condition.attrs['a']= a
                boundary_condition.attrs['b']= b
                boundary_condition.attrs['omega']= omega
                boundary_condition.attrs['phi']= phi
            elif bc_class == 'BoundaryCondition' :
                boundary_condition.attrs['v']= v
            elif bc_class == 'FixedBC' :
                pass # nothing to do
            else:
                raise NotImplementedError
