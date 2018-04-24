#!/usr/bin/env @PYTHON_EXECUTABLE@
"""
Description: Viewer and exporter for Siconos mechanics-IO HDF5 files based on VTK.
"""

# Lighter imports before command line parsing
from __future__ import print_function
import sys
import os
import json
import getopt
import math

# Exports from this module
__all__ = ['VView', 'VViewOptions', 'VExportOptions', 'VViewConfig']

if hasattr(math, 'inf'):
    infinity = math.inf
else:
    infinity = float('inf')

## Persistent configuration
class VViewConfig(dict):
    def __init__(self, d={'background_color' : [0., 0. , 0.],
                          'window_size': [600,600]}, filename=None):
        super(self.__class__, self).__init__(d)
        self.should_save_config = True
        if filename is not None:
            self.filename = filename
        else:
            self.filename = os.path.join(os.environ['HOME'], '.config',
                                         'siconos_vview.json')

    def load_configuration(self):
        if os.path.exists(self.filename):
            try:
                self.update(json.load(open(self.filename)))
                self.should_save_config = True
            except:
                self.should_save_config = False
                print("Warning: Error loading configuration `{}'".format(self.filename))

    def save_configuration(self, force=False):
        if not force and not self.should_save_config:
            return
        try:
            if not os.path.exists(os.path.join(os.environ['HOME'], '.config')):
                os.mkdir(os.path.join(os.environ['HOME'], '.config'))
            json.dump(self, open(self.filename,'w'))
        except:
            print("Error saving configuration `{}'".format(self.filename))


class VViewOptions(object):
    def __init__(self):
        self.min_time = None
        self.max_time = None
        self.cf_scale_factor = 1
        self.normalcone_ratio = 1
        self.time_scale_factor = 1
        self.advance_by_time = None
        self.frames_per_second = 25
        self.cf_disable = False
        self.initial_camera = [None] * 4
        self.visible_mode = 'all'

    ## Print usage information
    def usage(self, long=False):
        print(__doc__); print()
        print('Usage: {0} [OPTION]... <HDF5>'
              .format(os.path.split(sys.argv[0])[1]))
        print()
        if not long:
            print("""[--help] [--tmin=<float value>] [--tmax=<float value>]
            [--cf-scale=<float value>] [--no-cf] [--normalcone-ratio = <float value>]
            [--advance=<'fps' or float value>] [--fps=float value]
            [--camera=x,y,z] [--lookat=x,y,z] [--up=x,y,z] [--ortho=scale]
            [--visible=all,avatars,contactors]
            """)
        else:
            print("""Options:
     --help
       display this message
     --version
       display version information
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
     --visible=all
       all: view all contactors and avatars
       avatars: view only avatar if an avatar is defined (for each object)
       contactors: ignore avatars, view only contactors
         where avatars are contactors with collision_group=-1
    """)

    def parse(self):
        ## Parse command line
        try:
            opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                           ['help', 'version',
                                            'dat', 'tmin=', 'tmax=', 'no-cf',
                                            'cf-scale=', 'normalcone-ratio=',
                                            'advance=', 'fps=', 'camera=', 'lookat=',
                                            'up=', 'ortho=', 'visible='])
            self.configure(opts, args)
        except getopt.GetoptError as err:
            sys.stderr.write('{0}\n'.format(str(err)))
            self.usage()
            exit(2)

    def configure(self, opts, args):
        for o, a in opts:

            if o == '--help':
                self.usage(long=True)
                exit(0)

            elif o == '--version':
                print('{0} @SICONOS_VERSION@'.format(os.path.split(sys.argv[0])[1]))
                exit(0)

            elif o == '--tmin':
                self.min_time = float(a)

            elif o == '--tmax':
                self.max_time = float(a)

            elif o == '--cf-scale':
                self.cf_scale_factor = float(a)

            elif o == '--no-cf':
                self.cf_disable = True

            elif o == '--normalcone-ratio':
                self.normalcone_ratio = float(a)

            elif o == '--advance':
                if 'fps' in a:
                    self.advance_by_time = \
                        eval(a, {'fps': 1.0 / self.frames_per_second})
                else:
                    self.advance_by_time = float(a)

            elif o == '--fps':
                self.frames_per_second = int(a)

            elif o == '--camera':
                self.initial_camera[0] = map(float, a.split(','))

            elif o == '--lookat':
                self.initial_camera[1] = map(float, a.split(','))

            elif o == '--up':
                self.initial_camera[2] = map(float, a.split(','))

            elif o == '--ortho':
                self.initial_camera[3] = float(a)

            elif o == '--visible':
                self.visible_mode = a

        if self.frames_per_second == 0:
            self.frames_per_second = 25

        if len(args) > 0:
            self.io_filename = args[0]

        else:
            self.usage()
            exit(1)


class VExportOptions(VViewOptions):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ascii_mode = False

    def usage(self, long=False):
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

    def parse(self):
        ## Parse command line
        try:
            opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                           ['help','version','ascii'])
            self.configure(opts, args)
        except getopt.GetoptError as err:
                sys.stderr.write('{0}\n'.format(str(err)))
                self.usage()
                exit(2)

    def configure(self, opts, args):
        for o, a in opts:
            if o == '--help':
                self.usage(long=True)
                exit(0)
            if o == '--version':
                print('{0} @SICONOS_VERSION@'.format(os.path.split(sys.argv[0])[1]))
                exit(0)
            if o in ('--ascii'):
                self.ascii_mode = True

        if len(args) > 0:
            self.io_filename = args[0]

        else:
            self.usage()
            exit(1)


## Utilities

def add_compatiblity_methods(obj):
    """
    Add missing methods in previous VTK versions.
    """

    if hasattr(obj, 'SetInput'):
        obj.SetInputData = obj.SetInput

    if hasattr(obj, 'AddInput'):
        obj.AddInputData = obj.AddInput

def random_color():
    r = random.uniform(0.1, 0.9)
    g = random.uniform(0.1, 0.9)
    b = random.uniform(0.1, 0.9)
    return r, g, b


class Quaternion():
    def __init__(self, *args):
        import vtk
        self._vtkmath = vtk.vtkMath()
        self._data = vtk.vtkQuaternion[float](*args)

    def __mul__(self, q):
        r = Quaternion()
        self._vtkmath.MultiplyQuaternion(self._data, q._data, r._data)
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

        self.dom_at_time = [dict(),None][dom_data is None]
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

                    self._contact_field[mu].AddArray(self.cpa[mu])
                    self._contact_field[mu].AddArray(self.cpb[mu])
                    self._contact_field[mu].AddArray(self.cn[mu])
                    self._contact_field[mu].AddArray(self.cf[mu])

                    if dom_imu is not None:
                        self.dom_at_time[mu] = self._dom_data[
                            dom_id_f[imu], 1]
                        self.dom[mu] = numpy_support.numpy_to_vtk(
                            self.dom_at_time[mu])
                        self.dom[mu].SetName('domains')
                        self._contact_field[mu].AddArray(self.dom[mu])

                except KeyError:
                    pass

        else:
            pass

        # self._output.Update()


class InputObserver():

    def __init__(self, vview, times=None, slider_repres=None):
        self.vview = vview
        self._opacity = 1.0
        self._opacity_static = 1.0
        self._opacity_contact = 0.4
        self._current_id = vtk.vtkIdTypeArray()
        self._renderer = vview.renderer
        self._renderer_window = vview.renderer_window
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
        if self._times is None:
            self.vview.renderer_window.Render()
            return
        index = bisect.bisect_left(self._times, self._time)
        index = max(0, index)
        index = min(index, len(self._times) - 1)
        if self.vview.cf_prov is not None:
            self.vview.cf_prov._time = self._times[index]
            self.vview.cf_prov.xmethod()

        self.vview.set_dynamic_actors_visibility(self._times[index])

        id_t = numpy.where(self.vview.pos_data[:, 0] == self._times[index])
        self.vview.set_position(*self.vview.pos_data[id_t, :])

        self._slider_repres.SetValue(self._time)

        self._current_id.SetNumberOfValues(1)
        self._current_id.SetValue(0, index)

        self.vview.iter_plot.SetSelection(self._current_id)
        self.vview.prec_plot.SetSelection(self._current_id)

        self.vview.renderer_window.Render()

    def object_pos(self, id_):
        index = bisect.bisect_left(self._times, self._time)
        index = max(0, index)
        index = min(index, len(self._times) - 1)

        id_t = numpy.where(pos_data[:, 0] == self._times[index])
        return (pos_data[id_t[0][id_], 2], pos_data[id_t[0][id_], 3], pos_data[id_t[0][id_], 4])

    def set_opacity(self):
        for instance, actors in self.vview.dynamic_actors.items():
            for actor,_,_ in actors:
                actor.GetProperty().SetOpacity(self._opacity)

    def set_opacity_static(self):
        for instance, actors in self.vview.static_actors.items():
            for actor,_,_ in actors:
                actor.GetProperty().SetOpacity(self._opacity_static)

    def set_opacity_contact(self):
        for mu in self.vview.cf_prov._mu_coefs:
            self.vview.cactor[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.gactor[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.clactor[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.sactora[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.sactorb[mu].GetProperty().SetOpacity(self._opacity_contact)
                    

    def key(self, obj, event):
        key = obj.GetKeySym()
        print('key', key)

        if key == 'r':
            self.vview.reload()
            self._slider_repres.SetMinimumValue(self.vview.min_time)
            self._slider_repres.SetMaximumValue(self.vview.max_time)
            self.update()

        if key == 'p':
            self._image_counter += 1
            self.vview.image_maker.Update()
            self.vview.writer.SetFileName(
                'vview-{0}.png'.format(self._image_counter))
            self.vview.writer.Write()

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
            print('Decrease the opacity of bodies')
            self._opacity -= .1
            self.set_opacity()

        if key == 'T':
            print('Increase the opacity of bodies')
            self._opacity += .1
            self.set_opacity()
            
        if key == 'y':
            print('Decrease the opacity of static bodies')
            self._opacity_static -= .1
            self.set_opacity_static()

        if key == 'Y':
            print('Increase the opacity of static bodies')
            self._opacity_static += .1
            self.set_opacity_static()
            
        if key == 'u':
            print('Decrease the opacity of contact elements')
            self._opacity_contact -= .1
            self.set_opacity_contact()

        if key == 'U':
            print('Increase the opacity of contact elements')
            self._opacity_contact += .1
            self.set_opacity_contact()
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
            self.toggle_recording(True)

        if key == 'e':
            # Note 'e' has the effect to also "end" the program due to
            # default behaviour of vtkInteractorStyle, see class
            # documentation.
            self.toggle_recording(False)

        if key == 'C':
                this_view.action(self)
        self.update()

    def time(self, obj, event):
        slider_repres = obj.GetRepresentation()
        self._time = slider_repres.GetValue()
        self.update()

    def toggle_recording(self, recording):
        if recording and not self._recording:
            fps = 25
            self._timer_id = (self.vview.interactor_renderer
                              .CreateRepeatingTimer(1000//fps))
            self._recording = True
            self.vview.recorder.Start()
        elif self._recording and not recording:
            self.vview.interactor_renderer.DestroyTimer(self._timer_id)
            self._timer_id = None
            self._recording = False
            self.vview.recorder.End()

        # observer on 2D chart
    def iter_plot_observer(self, obj, event):

        if self.vview.iter_plot.GetSelection() is not None:
            # just one selection at the moment!
            if self.vview.iter_plot.GetSelection().GetMaxId() >= 0:
                self._time = self._times[
                    self.vview.iter_plot.GetSelection().GetValue(0)]
                # -> recompute index ...
                self.update()

    def prec_plot_observer(self, obj, event):
        if self.vview.prec_plot.GetSelection() is not None:
            # just one selection at the moment!
            if self.vview.prec_plot.GetSelection().GetMaxId() >= 0:
                self._time = self._times[
                    self.vview.prec_plot.GetSelection().GetValue(0)]
                # -> recompute index ...
                self.update()

    def recorder_observer(self, obj, event):
        if self._recording:
            if self.vview.opts.advance_by_time is not None:
                self._time += self.vview.opts.advance_by_time
                self.vview.slwsc.SetEnabled(False)  # Scale slider
                self.vview.xslwsc.SetEnabled(False)  # Time scale slider
                # slider_widget.SetEnabled(False) # Time slider (TODO video options)
                # widget.SetEnabled(False) # Axis widget
            self.update()
            self.vview.image_maker.Modified()
            self.vview.recorder.Write()
            if self.vview.opts.advance_by_time is not None:
                self.vview.slwsc.SetEnabled(True)
                self.vview.xslwsc.SetEnabled(True)
                # slider_widget.SetEnabled(True)
                # widget.SetEnabled(True) # Axis widget

            # End video if done
            if self._time >= max(self._times):
                self.toggle_recording(False)


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


# attach velocity
# contact points and associated forces are embedded in on a PolyData source
# (deferred class creation due to needing vtk to be imported)

def makeConvexSourceClass():
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
                self._convex.GetCellType(), self._convex.GetPointIds())
            output.SetPoints(self._points)
    return ConvexSource

# Read file and open VTK interaction window

class VView(object):
    def __init__(self, io, options, config=None):
        self.opts = options
        self.config = [config,VViewConfig()][config is None]
        self.gui_initialized = False
        self.io = io
        self.refs = []
        self.refs_attrs = []
        self.shape = dict()
        self.pos = dict()
        self.instances = dict()

        (self.spos_data, self.dpos_data, self.dom_data,
         self.cf_data, self.solv_data, self.velo_data) = self.load()

        self.contact_posa = dict()
        self.contact_posb = dict()
        self.contact_pos_force = dict()
        self.contact_pos_norm = dict()

        self.cone = dict()
        self.cone_glyph = dict()
        self.cmapper = dict()
        self.cLUT = dict()
        self.cactor = dict()
        self.arrow = dict()
        self.cylinder = dict()
        self.sphere = dict()
        self.arrow_glyph = dict()
        self.gmapper = dict()
        self.gactor = dict()
        self.ctransform = dict()
        self.cylinder_glyph = dict()
        self.clmapper = dict()
        self.sphere_glypha = dict()
        self.sphere_glyphb = dict()
        self.smappera = dict()
        self.smapperb = dict()
        self.sactora = dict()
        self.sactorb = dict()
        self.clactor = dict()
        self.data_connectors_v = dict()
        self.data_connectors_t = dict()
        self.data_connectors_d = dict()
        self.times_of_birth = dict()
        self.times_of_death = dict()
        self.min_time = self.opts.min_time
        self.max_time = self.opts.max_time

        self.transforms = dict()
        self.transformers = dict()

        self.offsets = dict()

    def load(self):

        ispos_data = self.io.static_data()
        idpos_data = self.io.dynamic_data()
        try:
            idom_data = self.io.domains_data()
        except ValueError:
            idom_data = None

        icf_data = self.io.contact_forces_data()[:]

        isolv_data = self.io.solver_data()
        ivelo_data = self.io.velocities_data()

        return ispos_data, idpos_data, idom_data, icf_data, isolv_data, ivelo_data

    def reload(self):
        (self.spos_data, self.dpos_data, self.dom_data,
         self.cf_data, self.solv_data, self.velo_data) = self.load()
        if not self.opts.cf_disable:
            self.cf_prov = CFprov(self.cf_data, self.dom_data)
        times = list(set(self.dpos_data[:, 0]))
        times.sort()

        if len(self.spos_data) > 0:
            self.instances = set(self.dpos_data[:, 1]).union(
                set(self.spos_data[:, 1]))
        else:
            self.instances = set(self.dpos_data[:, 1])

        if self.cf_prov is not None:
            self.cf_prov._time = min(times[:])
            self.cf_prov.xmethod()
            for mu in self.cf_prov._mu_coefs:
                self.contact_posa[mu].SetInputData(self.cf_prov._output[mu])
                self.contact_posa[mu].Update()
                self.contact_posb[mu].SetInputData(self.cf_prov._output[mu])
                self.contact_posb[mu].Update()
                self.contact_pos_force[mu].Update()
                self.contact_pos_norm[mu].Update()

        self.id_t0 = numpy.where(
            self.dpos_data[:, 0] == self.time0)

        self.pos_data = self.dpos_data[:]
        self.min_time = times[0]
        self.set_dynamic_actors_visibility(self.time0)

        self.max_time = times[len(times) - 1]

    def init_contact_pos(self, mu):

        self.contact_posa[mu] = vtk.vtkDataObjectToDataSetFilter()
        self.contact_posb[mu] = vtk.vtkDataObjectToDataSetFilter()

        add_compatiblity_methods(self.contact_posa[mu])
        add_compatiblity_methods(self.contact_posb[mu])

        self.contact_pos_force[mu] = vtk.vtkFieldDataToAttributeDataFilter()
        self.contact_pos_norm[mu] = vtk.vtkFieldDataToAttributeDataFilter()

        self.contact_posa[mu].SetDataSetTypeToPolyData()
        self.contact_posa[mu].SetPointComponent(0, "contactPositionsA", 0)
        self.contact_posa[mu].SetPointComponent(1, "contactPositionsA", 1)
        self.contact_posa[mu].SetPointComponent(2, "contactPositionsA", 2)

        self.contact_posb[mu].SetDataSetTypeToPolyData()
        self.contact_posb[mu].SetPointComponent(0, "contactPositionsB", 0)
        self.contact_posb[mu].SetPointComponent(1, "contactPositionsB", 1)
        self.contact_posb[mu].SetPointComponent(2, "contactPositionsB", 2)

        self.contact_pos_force[mu].SetInputConnection(
            self.contact_posa[mu].GetOutputPort())
        self.contact_pos_force[mu].SetInputFieldToDataObjectField()
        self.contact_pos_force[mu].SetOutputAttributeDataToPointData()
        self.contact_pos_force[mu].SetVectorComponent(0, "contactForces", 0)
        self.contact_pos_force[mu].SetVectorComponent(1, "contactForces", 1)
        self.contact_pos_force[mu].SetVectorComponent(2, "contactForces", 2)

        self.contact_pos_norm[mu].SetInputConnection(
            self.contact_posa[mu].GetOutputPort())
        self.contact_pos_norm[mu].SetInputFieldToDataObjectField()
        self.contact_pos_norm[mu].SetOutputAttributeDataToPointData()
        self.contact_pos_norm[mu].SetVectorComponent(0, "contactNormals", 0)
        self.contact_pos_norm[mu].SetVectorComponent(1, "contactNormals", 1)
        self.contact_pos_norm[mu].SetVectorComponent(2, "contactNormals", 2)

        if self.cf_prov.dom_at_time is not None:
            self.contact_pos_norm[mu].SetScalarComponent(0, "domains", 0)

    def init_cf_sources(self, mu, transform):
        self.contact_posa[mu].SetInputData(self.cf_prov._output[mu])
        self.contact_posa[mu].Update()
        self.contact_posb[mu].SetInputData(self.cf_prov._output[mu])
        self.contact_posb[mu].Update()

        self.contact_pos_force[mu].Update()
        self.contact_pos_norm[mu].Update()

        self.big_data_source.AddInputConnection(self.contact_posa[mu].GetOutputPort())
        self.big_data_source.AddInputConnection(self.contact_posb[mu].GetOutputPort())
        self.big_data_source.AddInputConnection(
            self.contact_pos_force[mu].GetOutputPort())
        self.big_data_source.AddInputConnection(
            self.contact_pos_norm[mu].GetOutputPort())

        self.cone[mu] = vtk.vtkConeSource()
        self.cone[mu].SetResolution(40)

        self.cone[mu].SetRadius(mu)  # one coef!!

        self.cone_glyph[mu] = vtk.vtkGlyph3D()
        self.cone_glyph[mu].SetSourceTransform(transform)

        self.cone_glyph[mu].SetInputConnection(self.contact_pos_norm[mu].GetOutputPort())
        self.cone_glyph[mu].SetSourceConnection(self.cone[mu].GetOutputPort())

        self.cone_glyph[mu]._scale_fact = self.opts.normalcone_ratio
        self.cone_glyph[mu].SetScaleFactor(
            self.cone_glyph[mu]._scale_fact *self.opts.cf_scale_factor)
        self.cone_glyph[mu].SetVectorModeToUseVector()

        self.cone_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
        self.cone_glyph[mu].OrientOn()

        # Don't allow scalar to affect size of glyph
        self.cone_glyph[mu].SetScaleModeToDataScalingOff()

        # Allow scalar to affect color of glyph
        self.cone_glyph[mu].SetColorModeToColorByScalar()

        self.cmapper[mu] = vtk.vtkPolyDataMapper()
        self.cmapper[mu].SetInputConnection(self.cone_glyph[mu].GetOutputPort())

        # Random color map, up to 256 domains
        self.cLUT[mu] = vtk.vtkLookupTable()
        self.cLUT[mu].SetNumberOfColors(256)
        self.cLUT[mu].Build()
        for i in range(256):
            self.cLUT[mu].SetTableValue(i, *random_color())
        self.cLUT[mu].SetTableRange(0, 255)

        # By default don't allow scalars to have an effect
        self.cmapper[mu].ScalarVisibilityOff()

        # If domain information is available, we turn on the color
        # table and turn on scalars
        if self.cf_prov.dom_at_time is not None:
            self.cmapper[mu].SetLookupTable(self.cLUT[mu])
            self.cmapper[mu].SetColorModeToMapScalars()
            self.cmapper[mu].SetScalarModeToUsePointData()
            self.cmapper[mu].SetScalarRange(0,255)
            self.cmapper[mu].ScalarVisibilityOn()

        self.cactor[mu] = vtk.vtkActor()
        
        self.cactor[mu].GetProperty().SetOpacity(self.config.get('contact_opacity', 0.4))
        self.cactor[mu].GetProperty().SetColor(0, 0, 1)
        self.cactor[mu].SetMapper(self.cmapper[mu])

        self.arrow[mu] = vtk.vtkArrowSource()
        self.arrow[mu].SetTipResolution(40)
        self.arrow[mu].SetShaftResolution(40)

        self.cylinder[mu] = vtk.vtkCylinderSource()
        self.cylinder[mu].SetRadius(.01)
        self.cylinder[mu].SetHeight(1)

        self.sphere[mu] = vtk.vtkSphereSource()

        # 1. scale = (scalar value of that particular data index);
        # 2. denominator = Range[1] - Range[0];
        # 3. scale = (scale < Range[0] ? Range[0] : (scale > Range[1] ? Range[1] : scale));
        # 4. scale = (scale - Range[0]) / denominator;
        # 5. scale *= scaleFactor;

        self.arrow_glyph[mu] = vtk.vtkGlyph3D()
        self.arrow_glyph[mu].SetInputConnection(
            self.contact_pos_force[mu].GetOutputPort())
        self.arrow_glyph[mu].SetSourceConnection(self.arrow[mu].GetOutputPort())
        self.arrow_glyph[mu].ScalingOn()
        self.arrow_glyph[mu].SetScaleModeToScaleByVector()
        self.arrow_glyph[mu].SetRange(0, .01)
        self.arrow_glyph[mu].ClampingOn()
        self.arrow_glyph[mu]._scale_fact = 5
        self.arrow_glyph[mu].SetScaleFactor(
            self.arrow_glyph[mu]._scale_fact * self.opts.cf_scale_factor)
        self.arrow_glyph[mu].SetVectorModeToUseVector()

        self.arrow_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactForces')
        self.arrow_glyph[mu].SetInputArrayToProcess(3, 0, 0, 0, 'contactForces')
        self.arrow_glyph[mu].OrientOn()

        self.gmapper[mu] = vtk.vtkPolyDataMapper()
        self.gmapper[mu].SetInputConnection(self.arrow_glyph[mu].GetOutputPort())
        self.gmapper[mu].SetScalarModeToUsePointFieldData()
        self.gmapper[mu].SetColorModeToMapScalars()
        self.gmapper[mu].ScalarVisibilityOn()
        self.gmapper[mu].SelectColorArray('contactForces')
        # gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

        self.gactor[mu] = vtk.vtkActor()
        self.gactor[mu].SetMapper(self.gmapper[mu])

        self.ctransform[mu] = vtk.vtkTransform()
        self.ctransform[mu].Translate(-0.5, 0, 0)
        self.ctransform[mu].RotateWXYZ(90, 0, 0, 1)
        self.cylinder_glyph[mu] = vtk.vtkGlyph3D()
        self.cylinder_glyph[mu].SetSourceTransform(self.ctransform[mu])

        self.cylinder_glyph[mu].SetInputConnection(
            self.contact_pos_norm[mu].GetOutputPort())
        self.cylinder_glyph[mu].SetSourceConnection(self.cylinder[mu].GetOutputPort())
        self.cylinder_glyph[mu].SetVectorModeToUseVector()

        self.cylinder_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
        self.cylinder_glyph[mu].OrientOn()
        self.cylinder_glyph[mu]._scale_fact = self.opts.normalcone_ratio
        self.cylinder_glyph[mu].SetScaleFactor(
            self.cylinder_glyph[mu]._scale_fact * self.opts.cf_scale_factor)

        self.clmapper[mu] = vtk.vtkPolyDataMapper()
        self.clmapper[mu].SetInputConnection(self.cylinder_glyph[mu].GetOutputPort())

        self.sphere_glypha[mu] = vtk.vtkGlyph3D()
        self.sphere_glypha[mu].SetInputConnection(self.contact_posa[mu].GetOutputPort())
        self.sphere_glypha[mu].SetSourceConnection(self.sphere[mu].GetOutputPort())
        self.sphere_glypha[mu].ScalingOn()
        # self.sphere_glypha[mu].SetScaleModeToScaleByVector()
        # self.sphere_glypha[mu].SetRange(-0.5, 2)
        # self.sphere_glypha[mu].ClampingOn()
        self.sphere_glypha[mu]._scale_fact = .1 * self.opts.normalcone_ratio
        self.sphere_glypha[mu].SetScaleFactor(
            self.sphere_glypha[mu]._scale_fact * self.opts.cf_scale_factor)
        # self.sphere_glypha[mu].SetVectorModeToUseVector()

        self.sphere_glyphb[mu] = vtk.vtkGlyph3D()
        self.sphere_glyphb[mu].SetInputConnection(self.contact_posb[mu].GetOutputPort())
        self.sphere_glyphb[mu].SetSourceConnection(self.sphere[mu].GetOutputPort())
        self.sphere_glyphb[mu].ScalingOn()
        # self.sphere_glyphb[mu].SetScaleModeToScaleByVector()
        # self.sphere_glyphb[mu].SetRange(-0.5, 2)
        # self.sphere_glyphb[mu].ClampingOn()
        self.sphere_glyphb[mu]._scale_fact = .1 * self.opts.normalcone_ratio
        self.sphere_glyphb[mu].SetScaleFactor(
            self.sphere_glyphb[mu]._scale_fact * self.opts.cf_scale_factor)
        # self.sphere_glyphb[mu].SetVectorModeToUseVector()

        # self.sphere_glyphb[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contactNormals')
        # self.sphere_glyph.OrientOn()

        self.smappera[mu] = vtk.vtkPolyDataMapper()
        self.smappera[mu].SetInputConnection(self.sphere_glypha[mu].GetOutputPort())
        self.smapperb[mu] = vtk.vtkPolyDataMapper()
        self.smapperb[mu].SetInputConnection(self.sphere_glyphb[mu].GetOutputPort())

        # self.cmapper.SetScalarModeToUsePointFieldData()
        # self.cmapper.SetColorModeToMapScalars()
        # self.cmapper.ScalarVisibilityOn()
        # self.cmapper.SelectColorArray('contactNormals')
        # self.gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contactForces').GetRange())

        self.clactor[mu] = vtk.vtkActor()
        # cactor.GetProperty().SetOpacity(0.4)
        self.clactor[mu].GetProperty().SetColor(1, 0, 0)
        self.clactor[mu].SetMapper(self.clmapper[mu])

        self.sactora[mu] = vtk.vtkActor()
        self.sactora[mu].GetProperty().SetColor(1, 0, 0)
        self.sactora[mu].SetMapper(self.smappera[mu])

        self.sactorb[mu] = vtk.vtkActor()
        self.sactorb[mu].GetProperty().SetColor(0, 1, 0)
        self.sactorb[mu].SetMapper(self.smapperb[mu])


    def init_shape(self, shape_name):

        shape_type = (self.io.shapes()[shape_name].attrs['type'])
        try:
            # work-around h5py unicode bug
            # https://github.com/h5py/h5py/issues/379
            shape_type  = shape_type.decode('utf-8')
        except AttributeError:
            pass
        scale = None
        if 'scale' in self.io.shapes()[shape_name].attrs:
            scale = self.io.shapes()[shape_name].attrs['scale']

        ConvexSource = makeConvexSourceClass()

        if shape_type in ['vtp', 'stl']:
            with io_tmpfile() as tmpf:
                tmpf[0].write(str(self.io.shapes()[shape_name][:][0]))
                tmpf[0].flush()
                reader = self.vtk_reader[shape_type]()
                reader.SetFileName(tmpf[1])
                reader.Update()
                self.readers[shape_name] = reader

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
                self.mappers[shape_name] = (x for x in [mapper])

        elif shape_type in ['brep']:
            # try to find an associated shape
            if 'associated_shape' in self.io.shapes()[shape_name].attrs:
                associated_shape = \
                    self.io.shapes()[shape_name].\
                    attrs['associated_shape']
                # delayed
                self.mappers[shape_name] = (x for x in
                                       [mappers[associated_shape]()])
            else:
                if 'brep' in self.io.shapes()[shape_name].attrs:
                    brep = self.io.shapes()[shape_name].attrs['brep']
                else:
                    brep = shape_name

                reader = brep_reader(str(self.io.shapes()[brep][:][0]),
                                     self.io.shapes()[brep].attrs['occ_indx'])
                self.readers[shape_name] = reader
                mapper = vtk.vtkDataSetMapper()
                add_compatiblity_methods(mapper)
                mapper.SetInputConnection(reader.GetOutputPort())
                self.mappers[shape_name] = (x for x in [mapper])

        elif shape_type in ['stp', 'step', 'igs', 'iges']:
            # try to find an associated shape
            if 'associated_shape' in self.io.shapes()[shape_name].attrs:
                associated_shape = \
                    self.io.shapes()[shape_name].\
                    attrs['associated_shape']
                # delayed
                self.mappers[shape_name] = (
                    x for x in [mappers[associated_shape]()])

            elif shape_type in ['stp', 'step', 'igs', 'iges']:
                with io_tmpfile(
                        debug=True,
                        suffix='.{0}'.format(shape_type),
                        contents=str(self.io.shapes()[shape_name][:][0])) as tmpf:
                    shape = occ_load_file(tmpf[1])

                    # whole shape
                    reader = topods_shape_reader(shape)
                    self.readers[shape_name] = reader
                    mapper = vtk.vtkDataSetMapper()
                    add_compatiblity_methods(mapper)
                    mapper.SetInputConnection(reader.GetOutputPort())
                    self.mappers[shape_name] = (x for x in [mapper])

                    # subparts
                    faces, edges = occ_topo_list(shape)
                    for i, f in enumerate(faces):
                        shape_indx = ('Face', shape_name, i)
                        reader = topods_shape_reader(f)
                        self.readers[shape_indx] = reader
                        mapper = vtk.vtkDataSetMapper()
                        add_compatiblity_methods(mapper)
                        mapper.SetInputConnection(reader.GetOutputPort())
                        self.mappers[shape_indx] = (x for x in [mapper])

                    for i, e in enumerate(edges):
                        shape_indx = ('Edge', shape_name, i)
                        reader = topods_shape_reader(e)
                        self.readers[shape_indx] = reader
                        mapper = vtk.vtkDataSetMapper()
                        add_compatiblity_methods(mapper)
                        mapper.SetInputConnection(reader.GetOutputPort())
                        self.mappers[shape_indx] = (x for x in [mapper])

        elif shape_type == 'heightmap':
            points = vtk.vtkPoints()
            shape = self.io.shapes()[shape_name]
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
            self.datasets[shape_name] = polydata
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(delaunay.GetOutputPort())
            add_compatiblity_methods(mapper)
            self.mappers[shape_name] = None
            self.mappers[shape_name] = (x for x in [mapper])

        elif shape_type == 'convex':
            # a convex shape
            points = vtk.vtkPoints()
            convex = vtk.vtkConvexPointSet()
            data = self.io.shapes()[shape_name][:]
            convex.GetPointIds().SetNumberOfIds(data.shape[0])

            for id_, vertice in enumerate(self.io.shapes()[shape_name][:]):
                points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
                convex.GetPointIds().SetId(id_, id_)

            source = ConvexSource(convex, points)
            self.readers[shape_name] = source

            # not a source!
            self.datasets[shape_name] = source.GetUnstructuredGridOutput()

            mapper = vtk.vtkDataSetMapper()
            add_compatiblity_methods(mapper)
            mapper.SetInputData(source.GetUnstructuredGridOutput())
            self.mappers[shape_name] = (x for x in [mapper])

        else:
            assert shape_type == 'primitive'
            primitive = self.io.shapes()[shape_name].attrs['primitive']
            attrs = self.io.shapes()[shape_name][:][0]
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

            self.readers[shape_name] = source
            mapper = vtk.vtkCompositePolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            self.mappers[shape_name] = (x for x in [mapper])

    def init_shapes(self):
        for shape_name in self.io.shapes():
            self.init_shape(shape_name)

        for shape_name in self.mappers.keys():
            if shape_name not in self.unfrozen_mappers:
                self.unfrozen_mappers[shape_name] = next(self.mappers[shape_name])

    def init_contactor(self, contactor_instance_name, instance, instid):
        contactor = instance[contactor_instance_name]
        contact_shape_indx = None
        if 'shape_name' not in contactor.attrs:
            print("Warning: old format: ctr.name must be ctr.shape_name for contact {0}".format(contactor_instance_name))
            shape_attr_name='name'
        else:
            shape_attr_name='shape_name'

        if 'group' in  contactor.attrs:
            collision_group = contactor.attrs['group']
        else:
            collision_group = -1
        if 'type' in contactor.attrs:
            contact_type = contactor.attrs['type']
            contact_index = contactor.attrs['contact_index']
            contact_shape_indx = (contact_type, contactor.attrs[shape_attr_name],
                                  contact_index)
        else:
            contact_shape_indx = contactor.attrs[shape_attr_name]
            try:
                # work-around h5py unicode bug
                # https://github.com/h5py/h5py/issues/379
                contact_shape_indx  = contact_shape_indx.decode('utf-8')
            except AttributeError:
                pass


        actor = vtk.vtkActor()
        if instance.attrs.get('mass', 0) > 0:
            # objects that may move
            self.dynamic_actors[instid].append((actor, contact_shape_indx,
                                                collision_group))

            actor.GetProperty().SetOpacity(
                self.config.get('dynamic_opacity', 0.7))
        else:
            # objects that are not supposed to move
            self.static_actors[instid].append((actor, contact_shape_indx,
                                               collision_group))

            actor.GetProperty().SetOpacity(
                self.config.get('static_opacity', 1.0))

        actor.GetProperty().SetColor(random_color())

        actor.SetMapper(self.unfrozen_mappers[contact_shape_indx])

        self.renderer.AddActor(actor)

        transform = vtk.vtkTransform()
        transformer = vtk.vtkTransformFilter()

        if contact_shape_indx in self.readers:
            transformer.SetInputConnection(
                self.readers[contact_shape_indx].GetOutputPort())
        else:
            transformer.SetInputData(self.datasets[contact_shape_indx])

        if isinstance(contact_shape_indx, tuple):
            contact_shape_name = contact_shape_indx[1]
        else:
            contact_shape_name = contact_shape_indx

        if 'scale' in self.io.shapes()[contact_shape_name].attrs:
            scale = self.io.shapes()[contact_shape_name].attrs['scale']
            scale_transform = vtk.vtkTransform()
            scale_transform.Scale(scale, scale, scale)
            scale_transform.SetInput(transform)
            transformer.SetTransform(scale_transform)
            actor.SetUserTransform(scale_transform)
        else:
            transformer.SetTransform(transform)
            actor.SetUserTransform(transform)

        self.transformers[contact_shape_indx] = transformer

        self.big_data_source.AddInputConnection(transformer.GetOutputPort())

        self.transforms[instid].append(transform)

        if 'center_of_mass' in instance.attrs:
            center_of_mass = instance.\
                             attrs['center_of_mass'].astype(float)
        else:
            center_of_mass = [0., 0., 0.]

        self.offsets[instid].append(
            (numpy.subtract(contactor.attrs['translation'].astype(float),
                            center_of_mass),
             contactor.attrs['orientation'].astype(float)))


        self.data_connectors_v[instid] = DataConnector(instid)
        self.data_connectors_v[instid]._connector.SetInputConnection(
            transformer.GetOutputPort())
        self.data_connectors_v[instid]._connector.Update()
        self.big_data_source.AddInputConnection(
            self.data_connectors_v[instid]._connector.GetOutputPort())

        self.data_connectors_t[instid] = DataConnector(instid,
                                                       data_name='translation',
                                                       data_size=3)
        self.data_connectors_t[instid]._connector.SetInputConnection(
            transformer.GetOutputPort())
        self.data_connectors_t[instid]._connector.Update()
        self.big_data_source.AddInputConnection(
            self.data_connectors_t[instid]._connector.GetOutputPort())

        self.data_connectors_d[instid] = DataConnector(instid,
                                                       data_name='displacement',
                                                       data_size=3)
        self.data_connectors_d[instid]._connector.SetInputConnection(
            transformer.GetOutputPort())
        self.data_connectors_d[instid]._connector.Update()
        self.big_data_source.AddInputConnection(
            self.data_connectors_d[instid]._connector.GetOutputPort())

    def init_instance(self, instance_name):
        instance = self.io.instances()[instance_name]
        instid = int(instance.attrs['id'])
        self.transforms[instid] = []
        self.offsets[instid] = []
        if 'time_of_birth' in instance.attrs:
            self.times_of_birth[instid] = instance.attrs['time_of_birth']
        if 'time_of_death' in instance.attrs:
            self.times_of_death[instid] = instance.attrs['time_of_death']

        if instid >= 0:
            self.dynamic_actors[instid] = list()
        else:
            self.static_actors[instid] = list()

        for contactor_instance_name in instance:
            self.init_contactor(contactor_instance_name, instance,  instid)

    def init_instances(self):
        for instance_name in self.io.instances():
            self.init_instance(instance_name)

    # this sets the position for all transforms associated to an instance
    def set_position_i(self, instance, q0, q1, q2, q3, q4, q5, q6):
        if (numpy.any(numpy.isnan([q0, q1, q2, q3, q4, q5, q6]))
           or numpy.any(numpy.isinf([q0, q1, q2, q3, q4, q5, q6]))):
            print('Bad position for object number', int(instance),' :',  q0, q1, q2, q3, q4, q5, q6)
            return

        q = Quaternion((q3, q4, q5, q6))

        for transform, offset in zip(self.transforms[instance],
                                     self.offsets[instance]):

            p = q.rotate(offset[0])

            r = q * Quaternion(offset[1])

            transform.Identity()
            transform.Translate(q0 + p[0], q1 + p[1], q2 + p[2])

            axis, angle = r.axisAngle()

            transform.RotateWXYZ(angle * 180. / pi,
                                 axis[0],
                                 axis[1],
                                 axis[2])

    def set_position(self, data):
        self.set_position_v(data[:, 1],
                            data[:, 2],
                            data[:, 3],
                            data[:, 4],
                            data[:, 5],
                            data[:, 6],
                            data[:, 7],
                            data[:, 8])

    def build_set_functions(self, dcv=None, dct=None, dcd=None):
        if dcv is None: dcv = self.data_connectors_v
        if dct is None: dct = self.data_connectors_t
        if dcd is None: dcd = self.data_connectors_d

        # the numpy vectorization is ok on column vectors for each args
        self.set_position_v = numpy.vectorize(self.set_position_i)

        # here the numpy vectorization is used with a column vector and a
        # scalar for the time arg
        self.set_visibility_v = numpy.vectorize(self.set_dynamic_instance_visibility)

        def set_velocity(instance, v0, v1, v2, v3, v4, v5):
            if instance in dcv:
                dcv[instance]._data[:] = [v0, v1, v2, v3, v4, v5]
                dcv[instance]._connector.Update()

        if dcv is not None:
            self.set_velocity_v = numpy.vectorize(set_velocity)

        def set_translation(instance, x0, x1, x2 ):
            if instance in dct:
                dct[instance]._data[:] = [x0, x1, x2]
                dct[instance]._connector.Update()

        if dct is not None:
            self.set_translation_v = numpy.vectorize(set_translation)

        def set_displacement(instance, x0, x1, x2 ):
            if instance in dcd:
                dcd[instance]._data[:] = [x0, x1, x2]
                dcd[instance]._connector.Update()

        if dcd is not None:
            self.set_displacement_v = numpy.vectorize(set_displacement)

    # set visibility for all actors associated to a dynamic instance
    def set_dynamic_instance_visibility(self, instance, time):
        tob = self.times_of_birth.get(instance, -1)
        tod = self.times_of_death.get(instance, infinity)
        has_avatar = False
        if self.opts.visible_mode=='avatars' or self.opts.visible_mode=='contactors':
            for actor, index, group in self.dynamic_actors[instance]:
                if group==-1:
                    has_avatar = True
                    break
        if (tob <= time and tod >= time):
            for actor, index, group in self.dynamic_actors[instance]:
                if not has_avatar or visible_mode=='all':
                    actor.VisibilityOn()
                elif visible_mode=='avatars' and group==-1 and has_avatar:
                    actor.VisibilityOn()
                elif visible_mode=='contactors' and group!=-1 and has_avatar:
                    actor.VisibilityOn()
                else:
                    actor.VisibilityOff()
        else:
            for actor, index, group in self.dynamic_actors[instance]:
                actor.VisibilityOff()

    def set_dynamic_actors_visibility(self, time):
        self.set_visibility_v(list(self.dynamic_actors.keys()), time)

    # callback maker for scale manipulation
    def make_scale_observer(self, glyphs):

        def scale_observer(obj, event):
            slider_repres = obj.GetRepresentation()
            scale_at_pos = slider_repres.GetValue()
            for glyph in glyphs:
                for k in glyph:
                    glyph[k].SetScaleFactor(
                        scale_at_pos * glyph[k]._scale_fact)

        return scale_observer

    # callback maker for time scale manipulation
    def make_time_scale_observer(self, time_slider_repres, time_observer):

        delta_time = self.max_time - self.min_time

        def time_scale_observer(obj, event):
            slider_repres = obj.GetRepresentation()
            time_scale_at_pos = 1. - slider_repres.GetValue()

            current_time = time_observer._time

            shift = (current_time - self.min_time) / delta_time

            xmin_time = self.min_time + time_scale_at_pos / 2. * delta_time
            xmax_time = self.max_time - time_scale_at_pos / 2. * delta_time

            xdelta_time = xmax_time - xmin_time

            new_mintime = max(self.min_time, current_time - xdelta_time)
            new_maxtime = min(self.max_time, current_time + xdelta_time)

            time_slider_repres.SetMinimumValue(new_mintime)
            time_slider_repres.SetMaximumValue(new_maxtime)

        return time_scale_observer

    # make a slider widget and its representation
    def make_slider(self, title, observer, interactor,
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

        background_color = self.config.get('background_color', [.0,.0,.0])
        reverse_background_color =numpy.ones(3) - background_color

        if (numpy.linalg.norm(background_color-reverse_background_color) < 0.2):
            reverse_background_color = numpy.ones(3)
        slider_repres.GetSliderProperty().SetColor(*reverse_background_color)
        slider_repres.GetTitleProperty().SetColor(*reverse_background_color);
        slider_repres.GetLabelProperty().SetColor(*reverse_background_color);
        slider_repres.GetTubeProperty().SetColor(*reverse_background_color);
        slider_repres.GetCapProperty().SetColor(*reverse_background_color);

        slider_widget = vtk.vtkSliderWidget()
        slider_widget.SetInteractor(interactor)
        slider_widget.SetRepresentation(slider_repres)
        slider_widget.KeyPressActivationOff()
        slider_widget.SetAnimationModeToAnimate()
        slider_widget.SetEnabled(True)
        slider_widget.AddObserver('InteractionEvent', observer)

        return slider_widget, slider_repres

    def setup_initial_position(self):
        self.time0 = None
        try:
            # Positions at first time step
            self.time0 = min(self.dpos_data[:, 0])
            self.id_t0 = numpy.where(self.dpos_data[:, 0] == self.time0)
            self.pos_t0 = self.pos_data[self.id_t0, 0:9]
        except ValueError:
            # this is for the case simulation hass not been ran and
            # time does not exists
            self.time0 = 0
            self.id_t0 = None
            self.pos_t0 = numpy.array([
                numpy.hstack(([0.,
                               self.io.instances()[k].attrs['id']]
                              ,self.io.instances()[k].attrs['translation']
                              ,self.io.instances()[k].attrs['orientation']))
                for k in self.io.instances()
                if self.io.instances()[k].attrs['id'] >= 0])

        if numpy.shape(self.spos_data)[0] > 0:
            self.set_position(self.spos_data)
            # static objects are always visible
            for instance, actors in self.static_actors.items():
                for actor,_,_ in actors:
                     actor.VisibilityOn()

        self.set_position(*self.pos_t0)

        self.set_dynamic_actors_visibility(self.time0)

    def setup_vtk_renderer(self):
        self.renderer_window.AddRenderer(self.renderer)
        self.interactor_renderer.SetRenderWindow(self.renderer_window)
        self.interactor_renderer.GetInteractorStyle(
        ).SetCurrentStyleToTrackballCamera()

        # http://www.itk.org/Wiki/VTK/Depth_Peeling ?

        # Use a render window with alpha bits (as initial value is 0 (false) ):
        self.renderer_window.SetAlphaBitPlanes(1)

        # Force to not pick a framebuffer with a multisample buffer ( as initial
        # value is 8):
        self.renderer_window.SetMultiSamples(0)

        # Choose to use depth peeling (if supported) (initial value is 0
        # (false) )
        self.renderer.SetUseDepthPeeling(1)

        # Set depth peeling parameters.
        self.renderer.SetMaximumNumberOfPeels(100)

        # Set the occlusion ratio (initial value is 0.0, exact image)
        self.renderer.SetOcclusionRatio(0.1)


        
        # Set the initial camera position and orientation if specified
        if self.opts.initial_camera[0] is not None:
            self.renderer.GetActiveCamera().SetPosition(*self.opts.initial_camera[0])
        if self.opts.initial_camera[1] is not None:
            self.renderer.GetActiveCamera().SetFocalPoint(*self.opts.initial_camera[1])
        if self.opts.initial_camera[2] is not None:
            self.renderer.GetActiveCamera().SetViewUp(*self.opts.initial_camera[2])
        if self.opts.initial_camera[3] is not None:
            self.renderer.GetActiveCamera().ParallelProjectionOn()
            self.renderer.GetActiveCamera().SetParallelScale(
                self.opts.initial_camera[3])

        self.image_maker = vtk.vtkWindowToImageFilter()
        self.image_maker.SetInput(self.renderer_window)

        self.recorder = vtk.vtkOggTheoraWriter()
        self.recorder.SetQuality(2)
        self.recorder.SetRate(self.opts.frames_per_second)
        self.recorder.SetFileName(os.path.splitext(self.opts.io_filename)[0]+'.avi')
        self.recorder.SetInputConnection(self.image_maker.GetOutputPort())

        self.writer = vtk.vtkPNGWriter()
        self.writer.SetInputConnection(self.image_maker.GetOutputPort())

        # Create a vtkLight, and set the light parameters.
        light = vtk.vtkLight()
        light.SetFocalPoint(0, 0, 0)
        light.SetPosition(0, 0, 500)
        # light.SetLightTypeToHeadlight()
        self.renderer.AddLight(light)

        hlight = vtk.vtkLight()
        hlight.SetFocalPoint(0, 0, 0)
        # hlight.SetPosition(0, 0, 500)
        hlight.SetLightTypeToHeadlight()
        self.renderer.AddLight(hlight)
        self.renderer.SetBackground(*self.config.get('background_color', [.0,.0,.0]))

    def setup_charts(self):
        # Warning! numpy support offer a view on numpy array
        # the numpy array must not be garbage collected!
        nxtime = self.solv_data[:, 0]
        nxiters = self.solv_data[:, 1]
        nprecs = self.solv_data[:, 2]
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
        self.iter_plot = chart_iter.AddPlot(vtk.vtkChart.LINE)
        self.iter_plot.SetLabel('Solver iterations')
        self.iter_plot.GetXAxis().SetTitle('time')
        self.iter_plot.GetYAxis().SetTitle('iterations')

        self.prec_plot = chart_prec.AddPlot(vtk.vtkChart.LINE)
        self.prec_plot.SetLabel('Solver precisions')
        self.prec_plot.GetXAxis().SetTitle('time')
        self.prec_plot.GetYAxis().SetTitle('precisions')

        add_compatiblity_methods(self.iter_plot)
        add_compatiblity_methods(self.prec_plot)

        self.iter_plot.SetInputData(table, 'time', 'iterations')
        self.prec_plot.SetInputData(table, 'time', 'precisions')
        self.iter_plot.SetWidth(5.0)
        self.prec_plot.SetWidth(5.0)
        self.iter_plot.SetColor(0, 255, 0, 255)
        self.prec_plot.SetColor(0, 255, 0, 255)

        tview_iter.GetInteractor().AddObserver('RightButtonReleaseEvent',
                                               self.input_observer.iter_plot_observer)

        tview_prec.GetInteractor().AddObserver('RightButtonReleaseEvent',
                                               self.input_observer.prec_plot_observer)

        # screen_size = self.renderer_window.GetScreenSize()
        self.renderer_window.SetSize(*self.config['window_size'])
        self.renderer_window.SetWindowName('vview: ' + self.opts.io_filename)
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
        self.tview_iter = tview_iter
        self.tview_prec = tview_prec

    def setup_sliders(self, times):
        if len(times) > 0:
            slider_repres = vtk.vtkSliderRepresentation2D()

            if self.min_time is None:
                self.min_time = times[0]
            if self.max_time is None:
                self.max_time = times[len(times) - 1]

            slider_repres.SetMinimumValue(self.min_time)
            slider_repres.SetMaximumValue(self.max_time)
            slider_repres.SetValue(self.min_time)
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
            
            background_color = self.config.get('background_color', [.0,.0,.0])
            reverse_background_color =numpy.ones(3) - background_color

            if (numpy.linalg.norm(background_color-reverse_background_color) < 0.2):
                reverse_background_color = numpy.ones(3)

            slider_repres.GetSliderProperty().SetColor(*reverse_background_color)
            slider_repres.GetTitleProperty().SetColor(*reverse_background_color);
            slider_repres.GetLabelProperty().SetColor(*reverse_background_color);
            slider_repres.GetTubeProperty().SetColor(*reverse_background_color);
            slider_repres.GetCapProperty().SetColor(*reverse_background_color);

            slider_widget = vtk.vtkSliderWidget()
            self.slider_widget = slider_widget
            slider_widget.SetInteractor(self.interactor_renderer)
            slider_widget.SetRepresentation(slider_repres)
            slider_widget.KeyPressActivationOff()
            slider_widget.SetAnimationModeToAnimate()
            slider_widget.SetEnabled(True)

            self.input_observer = InputObserver(self, times, slider_repres)
            slider_widget.AddObserver("InteractionEvent", self.input_observer.time)
        else:
            self.input_observer = InputObserver(self)
        self.interactor_renderer.AddObserver('KeyPressEvent', self.input_observer.key)

        self.interactor_renderer.AddObserver(
            'TimerEvent', self.input_observer.recorder_observer)

        if self.io.contact_forces_data().shape[0] > 0:
            self.slwsc, self.slrepsc = self.make_slider(
                'Scale',
                self.make_scale_observer([self.cone_glyph, self.cylinder_glyph,
                                          self.sphere_glypha, self.sphere_glyphb,
                                          self.arrow_glyph]),
                self.interactor_renderer,
                self.opts.cf_scale_factor, self.opts.cf_scale_factor -
                self.opts.cf_scale_factor / 2,
                self.opts.cf_scale_factor + self.opts.cf_scale_factor / 2,
                0.03, 0.03, 0.03, 0.7)

        if len(times) > 0:
            self.xslwsc, self.xslrepsc = self.make_slider(
                'Time scale',
                self.make_time_scale_observer(slider_repres,
                                              self.input_observer),
                self.interactor_renderer,
                self.opts.time_scale_factor, self.opts.time_scale_factor -
                self.opts.time_scale_factor / 2,
                self.opts.time_scale_factor + self.opts.time_scale_factor / 2,
                0.1, 0.9, 0.3, 0.9)

    def setup_axes(self):
        # display coordinates axes
        self.axes = vtk.vtkAxesActor()
        self.axes.SetTotalLength(1.0, 1.0, 1.0)
        self.widget = vtk.vtkOrientationMarkerWidget()
        # self.widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 )
        self.widget.SetOrientationMarker(self.axes)
        self.widget.SetInteractor(self.interactor_renderer)
        # self.widget.SetViewport( 0.0, 0.0, 40.0, 40.0 );
        self.widget.SetEnabled(True)
        self.widget.InteractiveOn()

    def export(self):
        big_data_writer = vtk.vtkXMLMultiBlockDataWriter()
        add_compatiblity_methods(big_data_writer)
        big_data_writer.SetInputConnection(self.big_data_source.GetOutputPort())

        times = list(set(self.dpos_data[:, 0]))
        times.sort()
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

            self.cf_prov._time = times[index]

            # fix: should be called by contact_source?
            self.cf_prov.xmethod()

            id_t = numpy.where(self.pos_data[:, 0] == times[index])

            if numpy.shape(self.spos_data)[0] > 0:
                self.set_position_v(self.spos_data[:, 1], self.spos_data[:, 2],
                              self.spos_data[:, 3],
                              self.spos_data[:, 4], self.spos_data[:, 5],
                              self.spos_data[:, 6],
                              self.spos_data[:, 7], self.spos_data[:, 8])

            self.set_position_v(
                self.pos_data[id_t, 1], self.pos_data[id_t, 2], self.pos_data[id_t, 3],
                self.pos_data[id_t, 4], self.pos_data[id_t, 5], self.pos_data[id_t, 6],
                self.pos_data[id_t, 7], self.pos_data[id_t, 8])

            id_tv = numpy.where(self.velo_data[:, 0] == times[index])[0]

            self.set_velocity_v(
                self.velo_data[id_tv, 1],
                self.velo_data[id_tv, 2],
                self.velo_data[id_tv, 3],
                self.velo_data[id_tv, 4],
                self.velo_data[id_tv, 5],
                self.velo_data[id_tv, 6],
                self.velo_data[id_tv, 7])

            self.set_translation_v(
                self.pos_data[id_t, 1],
                self.pos_data[id_t, 2],
                self.pos_data[id_t, 3],
                self.pos_data[id_t, 4],
            )

            big_data_writer.SetFileName('{0}-{1}.{2}'.format(
                os.path.splitext(os.path.basename(self.opts.io_filename))[0],
                index, big_data_writer.GetDefaultFileExtension()))
            big_data_writer.SetTimeStep(times[index])
            self.big_data_source.Update()

            if self.opts.ascii_mode:
                big_data_writer.SetDataModeToAscii()
            big_data_writer.Write()

    def initialize_vtk(self):
        self.big_data_source = vtk.vtkMultiBlockDataGroupFilter()
        add_compatiblity_methods(self.big_data_source)

        self.cf_prov = None
        if not self.opts.cf_disable:
            self.cf_prov = CFprov(self.cf_data, self.dom_data)
            for mu in self.cf_prov._mu_coefs:
                self.init_contact_pos(mu)

        times = self.dpos_data[:, 0]

        if (len(times) == 0):
            print('No dynamic data found!  Empty simulation.')

        if self.cf_prov is not None:
            self.cf_prov._time = min(times[:])
            self.cf_prov.xmethod()

        if self.cf_prov is not None:
            transform = vtk.vtkTransform()
            transform.Translate(-0.5, 0., 0.)
            for mu in self.cf_prov._mu_coefs:
                self.init_cf_sources(mu, transform)

        self.readers = dict()
        self.datasets = dict()
        self.mappers = dict()
        self.dynamic_actors = dict()
        self.static_actors = dict()
        self.vtk_reader = {'vtp': vtk.vtkXMLPolyDataReader,
                           'stl': vtk.vtkSTLReader}
        self.unfrozen_mappers = dict()

        self.pos_data = self.dpos_data[:]
        self.spos_data = self.spos_data[:]
        self.build_set_functions()

        self.renderer = vtk.vtkRenderer()
        self.renderer_window = vtk.vtkRenderWindow()
        self.interactor_renderer = vtk.vtkRenderWindowInteractor()

        self.init_shapes()
        self.init_instances()

        if self.cf_prov is not None:
            for mu in self.cf_prov._mu_coefs:
                self.renderer.AddActor(self.gactor[mu])
                self.renderer.AddActor(self.cactor[mu])

                self.renderer.AddActor(self.clactor[mu])
                self.renderer.AddActor(self.sactora[mu])
                self.renderer.AddActor(self.sactorb[mu])

        self.setup_initial_position()
        self.renderer.ResetCamera()

    def initialize_gui(self):

        self.setup_vtk_renderer()
        times = self.dpos_data[:, 0]
        self.setup_sliders(times)
        self.setup_charts()
        self.setup_axes()

        self.gui_initialized = True

    def run(self):
        self.initialize_vtk()
        self.initialize_gui()
        self.interactor_renderer.Start()

##
## Program starts
##

if __name__=='__main__':
    ## Persistent configuration
    config = VViewConfig()

    # Load it immediately
    config.load_configuration()

    # Parse command-line
    opts = VViewOptions()
    opts.parse()

# Heavier imports after command line parsing
import vtk
from vtk.util import numpy_support
from math import pi
import bisect
from numpy.linalg import norm
import numpy
import random

from siconos.io.mechanics_hdf5 import MechanicsHdf5
from siconos.io.mechanics_hdf5 import tmpfile as io_tmpfile
from siconos.io.mechanics_hdf5 import occ_topo_list, occ_load_file,\
    topods_shape_reader, brep_reader

if __name__=='__main__':
    ## Options and config already loaded above
    with MechanicsHdf5(io_filename=opts.io_filename, mode='r') as io:
        vview = VView(io, opts, config)
        vview.run()

    # Update configuration and save it
    config['window_size'] = vview.renderer_window.GetSize()
    config.save_configuration(force=False)
