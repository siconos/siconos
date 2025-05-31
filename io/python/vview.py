#!/usr/bin/env @Python_EXECUTABLE@
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
import traceback
import vtk
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.numpy_interface import dataset_adapter as dsa

from vtkmodules.vtkRenderingCore import vtkTextActor
import h5py

# Exports from this module
__all__ = ['VView', 'VViewOptions', 'VExportOptions', 'VViewConfig']

if hasattr(math, 'inf'):
    infinity = math.inf
else:
    infinity = float('inf')

def print_io_vview( *args, **kwargs):
    print('[io.vview]', *args, **kwargs)


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
                print_io_vview('Loaded configuration from ', self.filename)
                for k in self:
                    print_io_vview('  ', k,': ', self[k])
                self.should_save_config = True
            except:
                self.should_save_config = False
                print_io_vview("Warning: Error loading configuration `{}'".format(self.filename))

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
        if hasattr(vtk.vtkPolyDataMapper(), 'ImmediateModeRenderingOff'):
            self.imr = False
        else:
            # vtk 8
            self.imr = True
        self.depth_peeling = True
        self.maximum_number_of_peels = 100
        self.occlusion_ratio = 0.1
        self.global_filter = False
        self.initial_camera = [None] * 5
        self.visible_mode = 'all'
        self.export = False
        self.gen_para_script = False
        self.with_edges = False
        self.with_random_color = True
        self.with_charts= 0
        self.depth_2d=0.1
        self.verbose=0




    ## Print usage information
    def usage(self, long=False):
        print(__doc__); print()
        print('Usage: {0} [OPTION]... <HDF5>'
              .format(os.path.split(sys.argv[0])[1]))
        print()
        if not long:
            print("""[--help] [--tmin=<float value>] [--tmax=<float value>]
            [--cf-scale=<float value>] [--no-cf] [--imr] [--global-filter]
            [--no-depth-peeling] [--maximum-number-of-peels=<int value>]
            [--occlusion-ratio=<float value>]
            [--normalcone-ratio = <float value>]
            [--advance=<'fps' or float value>] [--fps=float value]
            [--camera=x,y,z] [--lookat=x,y,z] [--up=x,y,z] [--clipping=near,far] [--ortho=scale]
            [--with-charts=<int value>]
            [--visible=all,avatars,contactors] [--with-edges] [--verbose]
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
     --imr
       immediate-mode-rendering, use less memory at the price of
       slower rendering
     --global-filter (default : off)
       With export mode, concatenates all blocks in a big vtkPolyData.
       This option is for when the number of objects is huge.
       With vview, the display is done with only one vtk
       actor. Note that global-filter use a vtkCompositeDataGeometryFilter
       which is slow.
     --no-depth-peeling (default : on)
       do not use vtk depth peeling
     --maximum-number-of-peels= value
       maximum number of peels when depth peeling is on
     --occlusion-ratio= value
       occlusion-ratio when depth peeling is on
     --normalcone-ratio = value  (default : 1.0 )
       introduce a ratio between the representation of the contact
       forces arrows the normal cone and the contact points. useful
       when the contact forces are small with respect to the
       characteristic dimesion
     --advance= value or 'fps'
       automatically advance time during recording (default : don't
       advance)
     --fps= value
       frames per second of generated video (default 25)
     --camera=x,y,z
       initial position of the camera (default=above looking down)
       current camera position can be obtained with the key c
     --lookat=x,y,z
       initial direction to look (default=center of bounding box)
     --up=x,y,z
       initial up direction of the camera (default=y-axis)
     --clipping=x,y,z
       initial clippring range of the camera
     --ortho=scale
       start in ortho mode with given parallel scale
       (default=perspective)
     --with-charts=value
       display convergence charts
     --visible=all
       all: view all contactors and avatars
       avatars: view only avatar if an avatar is defined (for each
       object) contactors: ignore avatars, view only contactors where
       avatars are contactors with collision_group=-1
     --with-edges
       add edges in the rendering (experimental for primitives)
     --with-fixed-color
       use fixed color defined in the config file
     --depth-2d=<value>
       specify a depth for 2D objects
     --verbose=<int verbose_level>
    """)

    def parse(self):
        ## Parse command line
        try:
            opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                           ['help', 'version',
                                            'dat', 'tmin=', 'tmax=',
                                            'no-cf', 'imr', 'global-filter',
                                            'no-depth-peeling',
                                            'maximum-number-of-peels=',
                                            'occlusion-ratio=',
                                            'cf-scale=', 'normalcone-ratio=',
                                            'advance=', 'fps=',
                                            'camera=', 'lookat=', 'up=', 'clipping=', 'ortho=', 'visible=',
                                            'with-edges', 'with-fixed-color', 'with-charts=', 'depth-2d=', 'verbose='])
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
                print('{0} @SICONOS_VERSION@'.format(
                    os.path.split(sys.argv[0])[1]))
                exit(0)

            elif o == '--tmin':
                self.min_time = float(a)

            elif o == '--tmax':
                self.max_time = float(a)

            elif o == '--cf-scale':
                self.cf_scale_factor = float(a)

            elif o == '--no-cf':
                self.cf_disable = True

            elif o == '--imr':
                self.imr = True

            elif o == '--no-depth-peeling':
                self.depth_peeling=False

            elif o == '--maximum-number-of-peels':
                self.maximum_number_of_peels = int(a)

            elif o == '--occlusion-ratio':
                self.occlusion_ratio = float(a)

            elif o == '--global-filter':
                self.global_filter = True

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

            elif o == '--clipping':
                self.initial_camera[4] = map(float, a.split(','))

            elif o == '--ortho':
                self.initial_camera[3] = float(a)

            elif o == '--with-charts=':
                self.with_charts = int(a)

            elif o == '--depth-2d':
                self.depth_2d = float(a)

            elif o == '--visible':
                self.visible_mode = a

            elif o == '--with-edges':
                self.with_edges = True

            elif o == '--with-fixed-color':
                self.with_random_color = False

            elif o == '--verbose':
                self.verbose = int(a)


        if self.verbose >=1:
            self.display()


        if self.frames_per_second == 0:
            self.frames_per_second = 25

        if len(args) > 0:
            self.io_filename = args[0]

        else:
            self.usage()
            exit(1)

    def display(self):

        display_str = \
        """[io.VViewOptions] Display vview options:
        min_time : {0}
        min_time : {1}
        cf_disable : {2}
        cf_scale_factor : {3}
        normalcone_ratio {4}
        time_scale_factor {5}
        advance_by_time {6}
        frames_per_second {7}
        imr {8}
        depth_peeling {9}
        maximum_number_of_peels {10}
        occlusion_ratio {11}
        global_filter {12}
        initial_camera {13}
        visible_mode {14}
        export {15}
        gen_para_script {16}
        with_edges {17}
        with_random_color {18}
        with_charts {19}
        depth_2d {20}
        verbose {21}""".format(self.min_time,
                   self.max_time,
                   self.cf_disable,
                   self.cf_scale_factor,
                   self.normalcone_ratio,
                   self.time_scale_factor,
                   self.advance_by_time,
                   self.frames_per_second,
                   self.imr,
                   self.depth_peeling,
                   self.maximum_number_of_peels,
                   self.occlusion_ratio,
                   self.global_filter,
                   self.initial_camera,
                   self.visible_mode,
                   self.export,
                   self.gen_para_script,
                   self.with_edges,
                   self.with_random_color,
                   self.with_charts,
                   self.depth_2d,
                   self.verbose
                   )





        print(display_str)
        # self.cf_scale_factor = 1
        # self.normalcone_ratio = 1
        # self.time_scale_factor = 1
        # self.advance_by_time = None
        # self.frames_per_second = 25
        # self.cf_disable = False
        # if hasattr(vtk.vtkPolyDataMapper(), 'ImmediateModeRenderingOff'):
        #     self.imr = False
        # else:
        #     # vtk 8
        #     self.imr = True
        # self.depth_peeling = True
        # self.maximum_number_of_peels = 100
        # self.occlusion_ratio = 0.1
        # self.global_filter = False
        # self.initial_camera = [None] * 5
        # self.visible_mode = 'all'
        # self.export = False
        # self.gen_para_script = False
        # self.with_edges = False
        # self.with_random_color = True
        # self.with_charts= 0
        # self.depth=0.1
        # self.verbose=0




class VExportOptions(VViewOptions):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.export = True
        self.ascii_mode = False
        self.start_step = 0
        self.end_step = None
        self.stride = 1
        self.nprocs = 1
        self.gen_para_script = False

    def usage(self, long=False):
        print(__doc__); print()
        print('Usage:  {0} [--help] [--version] [--ascii] <HDF5>'
              .format(os.path.split(sys.argv[0])[1]))
        if long:
            print()
            print("""Options:
            --help               display this message
            --version            display version information
            --global-filter      one vtp file/time step
            --start-step=n       integer, set the first simulation time step
                                 number (default: 0)
            --end-step=n         integer, set the last simulation time step
                                 number (default: None)
            --stride=n           integer, set export time step/simulation time step
                                 (default: 1)
            --ascii              export file in ascii format
            --gen-para-script=n generation of a gnu parallel command for n processus
            """)

    def parse(self):
        ## Parse command line
        try:
            opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                           ['help', 'version', 'ascii',
                                            'start-step=', 'end-step=',
                                            'stride=', 'global-filter',
                                            'gen-para-script=',
                                            'depth-2d=', 'verbose='])
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
                print('{0} @SICONOS_VERSION@'.format(
                    os.path.split(sys.argv[0])[1]))
                exit(0)
            if o == '--global-filter':
                self.global_filter = True
            if o == '--start-step':
                self.start_step = int(a)
            if o == '--end-step':
                self.end_step = int(a)
            if o == '--stride':
                self.stride = int(a)
            if o == '--gen-para-script':
                self.gen_para_script = True
                self.nprocs = int(a)
            if o in ('--ascii'):
                self.ascii_mode = True
            if o in ('--depth-2d'):
                self.depth_2d = float(a)
            if o in ('--verbose'):
                self.verbose = int(a)

        if self.verbose >=1:
            self.display()

        if len(args) > 0:
            self.io_filename = args[0]

        else:
            self.usage()
            exit(1)

class VRawDataExportOptions(VViewOptions):
    def __init__(self, io_filename = None):
        super(self.__class__, self).__init__()
        self.export = True
        self._export_position = True
        self._export_velocity = True
        self._export_cf = False
        self._export_velocity_in_absolute_frame = False
        self.start_step = 0
        self.end_step = None
        self.stride = 1

        self.io_filename = io_filename
    def usage(self, long=False):
        print(__doc__); print()
        print('Usage:  {0} [--help]  <HDF5>'
              .format(os.path.split(sys.argv[0])[1]))
        if long:
            print()
            print("""Options:
            --help               display this message
            --version            display version information
            --start-step=n       integer, set the first simulation time step
                                 number (default: 0)
            --end-step=n         integer, set the last simulation time step
                                 number (default: None)
            --stride=n           integer, set export time step/simulation time step
                                 (default: 1)
            --no-export-position do not export position
            --no-export-velocity do not export position
            --export-cf          do export of contact friction data
            --export-velocity-in-absolute-frame          do export of contact friction data

            """)

    def parse(self):
        ## Parse command line
        try:
            opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                           ['help', 'version', 'ascii',
                                            'start-step=', 'end-step=',
                                            'stride=',
                                            'no-export-position',
                                            'no-export-velocity',
                                            'export-cf',
                                            'export-velocity-in-absolute-frame'])
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
                print('{0} @SICONOS_VERSION@'.format(
                    os.path.split(sys.argv[0])[1]))
                exit(0)
            if o == '--start-step':
                self.start_step = int(a)
            if o == '--end-step':
                self.end_step = int(a)
            if o == '--stride':
                self.stride = int(a)
            if o == '--no-export-position':
                self._export_position = False
            if o == '--no-export-velocity':
                self._export_velocity = False
            if o == '--export-cf':
                self._export_cf = True
            if o == '--export-velocity-in-absolute-frame':
                self._export_velocity_in_absolute_frame = True

        if self.io_filename is  None:
            if len(args) > 0 :
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
        self.vview.print_verbose('InputObserver -  update at time', self._time)

        self.vview.io_reader.SetTime(self._time)

        if self._times is None:
            self.vview.renderer_window.Render()
            return

        if not self.vview.opts.cf_disable:
            self.vview.io_reader.Update()
            for mu in self.vview.io_reader._mu_coefs:
                self.vview.contact_posa[mu].Update()
                self.vview.contact_posb[mu].Update()
                self.vview.contact_pos_force[mu].Update()
                self.vview.contact_pos_norm[mu].Update()

        self.vview.set_dynamic_actors_visibility(self.vview.io_reader._time)
        self.vview.set_static_actors_visibility(self.vview.io_reader._time)

        pos_data = self.vview.io_reader.pos_data

        self.vview.set_position(pos_data)

        self._slider_repres.SetValue(self.vview.io_reader._time)

        self._current_id.SetNumberOfValues(1)
        self._current_id.SetValue(0, self.vview.io_reader._index)
        if self.vview.opts.with_charts:
            self.vview.iter_plot.SetSelection(self._current_id)
            self.vview.prec_plot.SetSelection(self._current_id)

        self.vview.renderer_window.Render()

    def set_opacity(self):
        for instance, actors in self.vview.dynamic_actors.items():
            for actor,_,_ in actors:
                actor.GetProperty().SetOpacity(self._opacity)

    def set_opacity_static(self):
        for instance, actors in self.vview.static_actors.items():
            for actor,_,_ in actors:
                actor.GetProperty().SetOpacity(self._opacity_static)

    def set_opacity_contact(self):
        for mu in self.vview.io_reader._mu_coefs:
            self.vview.cactor[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.gactor[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.clactor[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.sactora[mu].GetProperty().SetOpacity(self._opacity_contact)
            self.vview.sactorb[mu].GetProperty().SetOpacity(self._opacity_contact)


    def key(self, obj, event):
        key = obj.GetKeySym()

        self.vview.print_verbose('InputObserver -  key', key)
        key_recognized = True
        if key == 'r':
            self.vview.reload()
            self._slider_repres.SetMinimumValue(self.vview.min_time)
            self._slider_repres.SetMaximumValue(self.vview.max_time)
            self.update()

        elif key == 'p':
            self._image_counter += 1
            self.vview.image_maker.Update()
            self.vview.writer.SetFileName(
                'vview-{0}.png'.format(self._image_counter))
            self.vview.writer.Write()

        elif key == 'Up':
            self._time_step = self._time_step * 2.
            self._time += self._time_step

        elif key == 'Down':
            self._time_step = self._time_step / 2.
            self._time -= self._time_step

        elif key == 'Left':
            self._time -= self._time_step

        elif key == 'Right':
            self._time += self._time_step

        elif key == 't':
            print('Decrease the opacity of bodies')
            self._opacity -= .1
            self.set_opacity()

        elif key == 'T':
            print('Increase the opacity of bodies')
            self._opacity += .1
            self.set_opacity()

        elif key == 'y':
            print('Decrease the opacity of static bodies')
            self._opacity_static -= .1
            self.set_opacity_static()

        elif key == 'Y':
            print('Increase the opacity of static bodies')
            self._opacity_static += .1
            self.set_opacity_static()

        elif key == 'u':
            print('Decrease the opacity of contact elements')
            self._opacity_contact -= .1
            self.set_opacity_contact()

        elif key == 'U':
            print('Increase the opacity of contact elements')
            self._opacity_contact += .1
            self.set_opacity_contact()
        elif key == 'c':
            print('camera position:', self._renderer.GetActiveCamera().GetPosition())
            print('camera focal point', self._renderer.GetActiveCamera().GetFocalPoint())
            print('camera clipping range', self._renderer.GetActiveCamera().GetClippingRange())
            print('camera up vector', self._renderer.GetActiveCamera().GetViewUp())
            if self._renderer.GetActiveCamera().GetParallelProjection() != 0:
                print('camera parallel scale', self._renderer.GetActiveCamera().GetParallelScale())

        elif key == 'o':
            self._renderer.GetActiveCamera().SetParallelProjection(
                1 - self._renderer.GetActiveCamera().GetParallelProjection())

        elif key == 'v':
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

        elif key == 'C':
            self._renderer.ResetCameraClippingRange()

        elif key == 's':
            self.toggle_recording(True)

        elif key == 'e':
            # Note 'e' has the effect to also "end" the program due to
            # default behaviour of vtkInteractorStyle, see class
            # documentation.
            self.toggle_recording(False)
        # else:
        #     self.vview.print_verbose('InputObserver -  key not recognized')
        #     key_recognized=False

        # if(key_recognized):
        self.update()

    def time(self, obj, event):
        slider_repres = obj.GetRepresentation()
        self._time = slider_repres.GetValue()
        self.update()

    def toggle_recording(self, recording):
        self.vview.print_verbose('InputObserver -  toggle recording')
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


class CellConnector(vtk.vtkProgrammableFilter):
    """
    Add Arrays to Cells
    """
    def __init__(self, instance, data_names, data_sizes):
        vtk.vtkProgrammableFilter.__init__(self)
        self._instance = instance
        self._data_names = data_names
        self._data_sizes = data_sizes
        self.SetExecuteMethod(self.method)
        self._datas = [numpy.zeros(s) for s in data_sizes]
        self._vtk_datas = [None]*len(data_sizes)
        self._index = list(enumerate(data_names))
        for i, data_name in self._index:
            self._vtk_datas[i] = vtk.vtkFloatArray()
            self._vtk_datas[i].SetName(data_name)
            self._vtk_datas[i].SetNumberOfComponents(data_sizes[i])

    def method(self):
        input = self.GetInput()
        output = self.GetOutput()
        output.ShallowCopy(input)
        ncells = output.GetNumberOfCells()

        for i, data_name in self._index:
            self._vtk_datas[i].SetNumberOfTuples(ncells)

            if output.GetCellData().GetArray(data_name) is None:
                output.GetCellData().AddArray(self._vtk_datas[i])

            data = self._datas[i]

            data_t = data[0:self._data_sizes[i]]

            for c in range(ncells):
                output.GetCellData().GetArray(data_name).SetTuple(c, data_t)



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


# attempt for a vtk reader
# only half the way, the reading part is ok but the output is only used
# in vview and export from python members
class IOReader(VTKPythonAlgorithmBase):

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
                                        nInputPorts=0,
                                        nOutputPorts=1,
                                        outputType='vtkPolyData')
        self._io = None
        self._with_contact_forces = False
        self.cf_data = None
        self.time = 0
        self.timestep = 0
        self.points = vtk.vtkPoints()

    def RequestInformation(self, request, inInfo, outInfo):

        info = outInfo.GetInformationObject(0)


        info.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(),
                 self._times,
                 len(self._times))

        info.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(),
                 [self._times[0], self._times[-1]], 2)

        return 1

    def RequestData(self, request, inInfo, outInfo):

        info = outInfo.GetInformationObject(0)
        output = vtk.vtkPolyData.GetData(outInfo)
        output.SetPoints(self.points)

        # The time step requested
        t = info.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())

        id_t = max(0, numpy.searchsorted(self._times, t, side='right') - 1)
        if id_t < len(self._indices)-1:
            self._id_t_m = list(range(self._indices[id_t],
                                      self._indices[id_t+1]))
        else:
            self._id_t_m = [self._indices[id_t]]

        self._time = self._times[id_t]
        self._index = id_t
        self.pos_data = self._idpos_data[self._id_t_m, :]
        self.velo_data = self._ivelo_data[self._id_t_m, :]

        static_id_t = max(0, numpy.searchsorted(self._static_times, t, side='right') - 1)

        if len(self._static_indices) >0:
            if static_id_t < len(self._static_indices)-1:
                self._static_id_t_m = list(range(self._static_indices[static_id_t],
                                                 self._static_indices[static_id_t+1]))
            else:
                self._static_id_t_m = [self._static_indices[static_id_t]]
        else:
            self._static_id_t_m = None
        if self._static_id_t_m:
            self.pos_static_data = self._ispos_data[self._static_id_t_m, :]
        else:
            self.pos_static_data = None

        vtk_pos_data = dsa.numpyTovtkDataArray(self.pos_data)
        vtk_pos_data.SetName('pos_data')

        vtk_velo_data = dsa.numpyTovtkDataArray(self.velo_data)
        vtk_velo_data.SetName('velo_data')

        vtk_points_data = dsa.numpyTovtkDataArray(self.pos_data[:, 2:5])

        self.points.SetData(vtk_points_data)

        output.GetPointData().AddArray(vtk_velo_data)
        try:
            if self._with_contact_forces:
                ncfindices = len(self._cf_indices)
                id_t_cf = min(numpy.searchsorted(self._cf_times, t,
                                                 side='right'),
                              ncfindices-1)
                # Check the duration between t and last impact.
                # If it is superior to current time step, we consider there
                # is no contact (rebound).
                # The current time step is the max between slider timestep
                # and simulation timestep
                ctimestep = max(self.timestep, self._avg_timestep)
                if (id_t_cf > 0 and abs(t-self._cf_times[id_t_cf-1])
                    <= ctimestep):
                    if id_t_cf < ncfindices-1:
                        self._id_t_m_cf = list(range(self._cf_indices[id_t_cf-1],
                                                     self._cf_indices[id_t_cf]))
                        self.cf_data = self._icf_data[self._id_t_m_cf, :]

                    else:
                        self.cf_data = self._icf_data[self._cf_indices[
                            id_t_cf]:, :]

                    self._cf_time = self._cf_times[id_t_cf]

                    vtk_cf_data = dsa.numpyTovtkDataArray(self.cf_data)
                    vtk_cf_data.SetName('cf_data')

                    output.GetFieldData().AddArray(vtk_cf_data)
                else:
                    # there is no contact forces at this time
                    self.cf_data = None
                    vtk_cf_data = dsa.numpyTovtkDataArray(numpy.array([]))
                    vtk_cf_data.SetName('cf_data')
                    output.GetFieldData().AddArray(vtk_cf_data)

            if self.cf_data is not None:
                self.contact = True

                data = self.cf_data

                for mu in self._mu_coefs:

                    imu = numpy.where(
                        abs(data[:, 1] - mu) < 1e-15)[0]

                    #dom_imu = None
                    #dom_imu = numpy.where(
                    #    self._dom_data[:,-1] == data[id_f[imu],-1]
                    #)[0]

                    if len(imu) > 0:
                        self.cpa_at_time[mu] = data[
                            imu, 2:5]
                        self.cpb_at_time[mu] = data[
                            imu, 5:8]
                        self.cn_at_time[mu] = - data[
                            imu, 8:11]
                        self.cf_at_time[mu] = data[
                            imu, 11:14]
                        if data[imu, :].shape[1] > 26:
                            self.ids_at_time[mu] = data[
                                imu, 23:26].astype(int)
                        else:
                            self.ids_at_time[mu] = None

            else:
                try:
                    for mu in self._mu_coefs:
                        self.cpa_at_time[mu] = [[nan, nan, nan]]
                        self.cpb_at_time[mu] = [[nan, nan, nan]]
                        self.cn_at_time[mu] =  [[nan, nan, nan]]
                        self.cf_at_time[mu] =  [[nan, nan, nan]]
                        self.ids_at_time[mu] = None
                except:
                    pass

            if self._with_contact_forces:
                for mu in self._mu_coefs:
                    self.cpa[mu] = numpy_support.numpy_to_vtk(
                        self.cpa_at_time[mu])
                    self.cpa[mu].SetName('contact_positions_a')

                    self.cpb[mu] = numpy_support.numpy_to_vtk(
                        self.cpb_at_time[mu])
                    self.cpb[mu].SetName('contact_positions_b')

                    self.cn[mu] = numpy_support.numpy_to_vtk(
                        self.cn_at_time[mu])
                    self.cn[mu].SetName('contact_normals')

                    self.cf[mu] = numpy_support.numpy_to_vtk(
                        self.cf_at_time[mu])
                    self.cf[mu].SetName('contact_forces')

                    # field info for vview (should go in point data)
                    self._contact_field[mu].AddArray(self.cpa[mu])
                    self._contact_field[mu].AddArray(self.cpb[mu])
                    self._contact_field[mu].AddArray(self.cn[mu])
                    self._contact_field[mu].AddArray(self.cf[mu])

                    # contact points
                    self._points[mu].SetData(self.cpa[mu])
                    self._output[mu].GetPointData().AddArray(self.cpb[mu])
                    self._output[mu].GetPointData().AddArray(self.cn[mu])
                    self._output[mu].GetPointData().AddArray(self.cf[mu])

                    if self.ids_at_time[mu] is not None:
                        self.ids[mu] = numpy_support.numpy_to_vtk(
                            self.ids_at_time[mu])
                        self.ids[mu].SetName('ids')
                        self._contact_field[mu].AddArray(self.ids[mu])
                        self._output[mu].GetPointData().AddArray(self.ids[mu])

                        dsa_ids = numpy.unique(self.ids_at_time[mu][:, 1])
                        dsb_ids = numpy.unique(self.ids_at_time[mu][:, 2])
                        _i, _i, dsa_pos_ids = numpy.intersect1d(
                            self.pos_data[:, 1],
                            dsa_ids, return_indices=True)
                        _i, _i, dsb_pos_ids = numpy.intersect1d(
                            self.pos_data[:, 1],
                            dsb_ids, return_indices=True)

                        # objects a & b translations
                        obj_pos_a = self.pos_data[dsa_pos_ids, 2:5]
                        obj_pos_b = self.pos_data[dsb_pos_ids, 2:5]

                        self._all_objs_pos[mu] = numpy.vstack((obj_pos_a,
                                                               obj_pos_b))

                        self._all_objs_pos_vtk[mu] = numpy_support.numpy_to_vtk(
                            self._all_objs_pos[mu])
                        self._objs_points[mu].SetData(self._all_objs_pos_vtk[mu])

                        self._objs_output[mu].GetPointData().AddArray(self.cn[mu])
                        self._objs_output[mu].GetPointData().AddArray(self.cf[mu])


                        #if dom_imu is not None:
                        #    self.dom_at_time[mu] = self._dom_data[
                        #        dom_imu, 1]
                        #    self.dom[mu] = numpy_support.numpy_to_vtk(
                        #        self.dom_at_time[mu])
                        #    self.dom[mu].SetName('domains')
                        #    self._contact_field[mu].AddArray(self.dom[mu])

        except Exception:
            traceback.print_exc()

        return 1

    def SetIO(self, io):
        self._io = io

        self._ispos_data = self._io.static_data()
        self._idpos_data = self._io.dynamic_data()
        try:
            self._idom_data = self._io.domains_data()
        except ValueError:
            self._idom_data = None

        self._icf_data = self._io.contact_forces_data()
        self._isolv_data = self._io.solver_data()
        self._ivelo_data = self._io.velocities_data()

        self._spos_data = self._ispos_data[:, :]

        # all times as hdf5 slice
        self._raw_times = self._idpos_data[:, 0]

        # build times steps
        self._times, self._indices = numpy.unique(self._raw_times,
                                                  return_index=True)
        # all times as hdf5 slice for static objects
        self._static_raw_times = self._ispos_data[:, 0]

        # build times steps for static objects
        self._static_times, self._static_indices = numpy.unique(self._static_raw_times,
                                                                return_index=True)
        dcf = self._times[1:]-self._times[:-1]
        self._avg_timestep = numpy.mean(dcf)
        self._min_timestep = numpy.min(dcf)

        # self._times.sort()
        # self._indices = ?
        # we assume times must be sorted
#        assert all(self._times[i] <= self._times[i+1]
#                   for i in range(len(self._times)-1))

#        if self._with_contact_forces:
#            assert all(self._cf_times[i] <= self._cf_times[i+1]
#                       for i in range(len(self._cf_times)-1))

        self.Modified()
        return 1

    def SetTime(self, time):
        self.GetOutputInformation(0).Set(
            vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(),
            time)
        self.timestep = abs(self.time-time)
        self.time = time

        # with a True pipeline: self.Modified()
        # but as the consumers (VView class, export function) are
        # not (yet) vtk filters, the Update is needed here
        self.Update()

        # contact forces provider
    def ContactForcesOn(self):

        self._cf_raw_times = self._icf_data[:, 0]
        self._cf_times, self._cf_indices = numpy.unique(self._cf_raw_times,
                                                        return_index=True)
        self._mu_coefs = numpy.unique(self._icf_data[:, 1],
                                      return_index=False)

        self.cpa_at_time = dict()
        self.cpa = dict()

        self.cpb_at_time = dict()
        self.cpb = dict()

        self.cf_at_time = dict()
        self.cf = dict()

        self.cn_at_time = dict()
        self.cn = dict()

        self.ids_at_time = dict()
        self.ids = dict()

        self.dom_at_time = [dict(), None][self._idom_data is None]
        self.dom = dict()

        self._all_objs_pos = dict()
        self._all_objs_pos_vtk = dict()

        self._points = dict()
        self._contact_field = dict()
        self._output = dict()
        self._objs_points = dict()
        self._objs_output = dict()

        for mu in self._mu_coefs:
            # the contact points
            self._points[mu] = vtk.vtkPoints()
            self._contact_field[mu] = vtk.vtkPointData()
            self._output[mu] = vtk.vtkPolyData()
            self._output[mu].SetPoints(self._points[mu])
            self._output[mu].SetFieldData(self._contact_field[mu])

            # the objects translations
            self._objs_points[mu] = vtk.vtkPoints()
            self._objs_output[mu] = vtk.vtkPolyData()
            self._objs_output[mu].SetPoints(self._objs_points[mu])

        self._with_contact_forces = True
        self.Update()

    def ContactForcesOff(self):
        self._with_contact_forces = False
        self.Update()

    def ExportOn(self):
        self._export = True

    def ExportOff(self):
        self._export = False


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
        self.mass = dict()
        self.inertia = dict()

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
        self.cell_connectors = dict()
        self.times_of_birth = dict()
        self.times_of_death = dict()
        self.min_time = self.opts.min_time
        self.max_time = self.opts.max_time

        self.transforms = dict()
        self.transformers = dict()

        self.offsets = dict()

        self.io_reader = IOReader()

        self.io_reader.SetIO(io=self.io)

        if self.opts.cf_disable:
            self.io_reader.ContactForcesOff()
        else:
            self.io_reader.ContactForcesOn()

        if self.opts.export:
            self.io_reader.ExportOn()
        else:
            self.io_reader.ExportOff()

    def print_verbose(self, *args, **kwargs):
        if self.opts.verbose:
            print('[io.vview]', *args, **kwargs)

    def print_verbose_level(self, level, *args, **kwargs):
        if level <= self.opts.verbose:
            print('[io.vview]', *args, **kwargs)

    def reload(self):

        if self.opts.cf_disable:
            self.io_reader.ContactForcesOff()

        else:
            self.io_reader._time = min(times[:])
            for mu in self.io_reader._mu_coefs:
                self.contact_posa[mu].SetInputData(self.io_reader._output[mu])
                self.contact_posa[mu].Update()
                self.contact_posb[mu].SetInputData(self.io_reader._output[mu])
                self.contact_posb[mu].Update()
                self.contact_pos_force[mu].Update()
                self.contact_pos_norm[mu].Update()

        self.min_time = self.io_reader._times[0]
        self.max_time = self.io_reader._times[-1]
        self.set_dynamic_actors_visibility(self.time0)

    def init_contact_pos(self, mu):
        self.print_verbose_level(2,'contact positions', mu)
        self.contact_posa[mu] = vtk.vtkDataObjectToDataSetFilter()
        self.contact_posb[mu] = vtk.vtkDataObjectToDataSetFilter()

        add_compatiblity_methods(self.contact_posa[mu])
        add_compatiblity_methods(self.contact_posb[mu])

        self.contact_pos_force[mu] = vtk.vtkFieldDataToAttributeDataFilter()
        self.contact_pos_norm[mu] = vtk.vtkFieldDataToAttributeDataFilter()

        self.contact_posa[mu].SetDataSetTypeToPolyData()
        self.contact_posa[mu].SetPointComponent(0, "contact_positions_a", 0)
        self.contact_posa[mu].SetPointComponent(1, "contact_positions_a", 1)
        self.contact_posa[mu].SetPointComponent(2, "contact_positions_a", 2)

        self.contact_posb[mu].SetDataSetTypeToPolyData()
        self.contact_posb[mu].SetPointComponent(0, "contact_positions_b", 0)
        self.contact_posb[mu].SetPointComponent(1, "contact_positions_b", 1)
        self.contact_posb[mu].SetPointComponent(2, "contact_positions_b", 2)

        self.contact_pos_force[mu].SetInputConnection(
            self.contact_posa[mu].GetOutputPort())
        self.contact_pos_force[mu].SetInputFieldToDataObjectField()
        self.contact_pos_force[mu].SetOutputAttributeDataToPointData()
        self.contact_pos_force[mu].SetVectorComponent(0, "contact_forces", 0)
        self.contact_pos_force[mu].SetVectorComponent(1, "contact_forces", 1)
        self.contact_pos_force[mu].SetVectorComponent(2, "contact_forces", 2)

        self.contact_pos_norm[mu].SetInputConnection(
            self.contact_posa[mu].GetOutputPort())
        self.contact_pos_norm[mu].SetInputFieldToDataObjectField()
        self.contact_pos_norm[mu].SetOutputAttributeDataToPointData()
        self.contact_pos_norm[mu].SetVectorComponent(0, "contact_normals", 0)
        self.contact_pos_norm[mu].SetVectorComponent(1, "contact_normals", 1)
        self.contact_pos_norm[mu].SetVectorComponent(2, "contact_normals", 2)

#       if self.cf_prov.dom_at_time is not None:
#            self.contact_pos_norm[mu].SetScalarComponent(0, "domains", 0)

    def init_cf_sources(self, mu, transform):
        self.print_verbose_level(2,'contact sources', mu)
        self.cf_collector.AddInputData(self.io_reader._output[mu])
        self.cf_collector.AddInputData(self.io_reader._objs_output[mu])
        self.contact_posa[mu].SetInputData(self.io_reader._output[mu])
        self.contact_posa[mu].Update()
        self.contact_posb[mu].SetInputData(self.io_reader._output[mu])
        self.contact_posb[mu].Update()

        self.contact_pos_force[mu].Update()
        self.contact_pos_norm[mu].Update()

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

        self.cone_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contact_normals')
        self.cone_glyph[mu].OrientOn()

        # Don't allow scalar to affect size of glyph
        self.cone_glyph[mu].SetScaleModeToDataScalingOff()

        # Allow scalar to affect color of glyph
        self.cone_glyph[mu].SetColorModeToColorByScalar()

        self.cmapper[mu] = vtk.vtkPolyDataMapper()
        if not self.opts.imr:
            self.cmapper[mu].ImmediateModeRenderingOff()
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
        if self.io_reader.dom_at_time is not None:
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

        self.arrow_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contact_forces')
        self.arrow_glyph[mu].SetInputArrayToProcess(3, 0, 0, 0, 'contact_forces')
        self.arrow_glyph[mu].OrientOn()

        self.gmapper[mu] = vtk.vtkPolyDataMapper()
        if not self.opts.imr:
            self.gmapper[mu].ImmediateModeRenderingOff()
        self.gmapper[mu].SetInputConnection(self.arrow_glyph[mu].GetOutputPort())
        self.gmapper[mu].SetScalarModeToUsePointFieldData()
        self.gmapper[mu].SetColorModeToMapScalars()
        self.gmapper[mu].ScalarVisibilityOn()
        self.gmapper[mu].SelectColorArray('contact_forces')
        # gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contact_forces').GetRange())

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

        self.cylinder_glyph[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contact_normals')
        self.cylinder_glyph[mu].OrientOn()
        self.cylinder_glyph[mu]._scale_fact = self.opts.normalcone_ratio
        self.cylinder_glyph[mu].SetScaleFactor(
            self.cylinder_glyph[mu]._scale_fact * self.opts.cf_scale_factor)

        self.clmapper[mu] = vtk.vtkPolyDataMapper()
        if not self.opts.imr:
            self.clmapper[mu].ImmediateModeRenderingOff()
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

        # self.sphere_glyphb[mu].SetInputArrayToProcess(1, 0, 0, 0, 'contact_normals')
        # self.sphere_glyph.OrientOn()

        self.smappera[mu] = vtk.vtkPolyDataMapper()
        if not self.opts.imr:
            self.smappera[mu].ImmediateModeRenderingOff()
        self.smappera[mu].SetInputConnection(self.sphere_glypha[mu].GetOutputPort())
        self.smapperb[mu] = vtk.vtkPolyDataMapper()
        if not self.opts.imr:
            self.smapperb[mu].ImmediateModeRenderingOff()
        self.smapperb[mu].SetInputConnection(self.sphere_glyphb[mu].GetOutputPort())

        # self.cmapper.SetScalarModeToUsePointFieldData()
        # self.cmapper.SetColorModeToMapScalars()
        # self.cmapper.ScalarVisibilityOn()
        # self.cmapper.SelectColorArray('contact_normals')
        # self.gmapper.SetScalarRange(contact_pos_force.GetOutput().GetPointData().GetArray('contact_forces').GetRange())

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
        self.print_verbose_level(2,'init_shape', shape_name)
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
                ## fix compatibility with h5py version: to be removed in the future
                if (h5py.version.version_tuple.major >=3 ):
                    tmpf[0].write((self.io.shapes()[shape_name][:][0]).decode('utf-8'))
                else:
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
                        contents=(self.io.shapes()[shape_name][:][0]).decode("utf-8")) as tmpf:
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
            if not self.opts.imr:
                mapper.ImmediateModeRenderingOff()
            mapper.SetInputConnection(delaunay.GetOutputPort())
            add_compatiblity_methods(mapper)
            self.mappers[shape_name] = None
            self.mappers[shape_name] = (x for x in [mapper])

        elif shape_type == 'convex':
            # a convex shape
            points = vtk.vtkPoints()
            convex = vtk.vtkConvexPointSet()
            data = self.io.shapes()[shape_name][:]
            if self.io.dimension() == 3:
                convex.GetPointIds().SetNumberOfIds(data.shape[0])
                for id_, vertice in enumerate(data):
                    points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
                    convex.GetPointIds().SetId(id_, id_)
            elif self.io.dimension() == 2:
                number_of_vertices = data.shape[0]
                convex.GetPointIds().SetNumberOfIds(data.shape[0]*2)
                for id_, vertice in enumerate(data):
                    points.InsertNextPoint(vertice[0], vertice[1], - self.opts.depth_2d/2.0)
                    convex.GetPointIds().SetId(id_, id_)
                    points.InsertNextPoint(vertice[0], vertice[1], + self.opts.depth_2d/2.0)
                    convex.GetPointIds().SetId(id_+number_of_vertices, id_+number_of_vertices)

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

            elif primitive == 'Disk':
                cylinder = vtk.vtkCylinderSource()
                cylinder.SetResolution(100)
                cylinder.SetRadius(attrs[0])
                cylinder.SetHeight(self.opts.depth_2d)
                cylinder.Update()

                line1 = vtk.vtkLineSource()
                line1.SetPoint1(attrs[0]/2, self.opts.depth_2d, 0)
                line1.SetPoint2(attrs[0], self.opts.depth_2d, 0)

                tube1 = vtk.vtkTubeFilter()
                tube1.SetInputConnection(line1.GetOutputPort())
                tube1.SetRadius(attrs[0]/20)
                tube1.SetNumberOfSides(10)
                tube1.Update()

                data = vtk.vtkMultiBlockDataSet()
                data.SetNumberOfBlocks(2)
                data.SetBlock(0, tube1.GetOutput())
                data.SetBlock(1, cylinder.GetOutput())

                source = vtk.vtkMultiBlockDataGroupFilter()
                add_compatiblity_methods(source)
                source.AddInputData(data)

            elif primitive == 'Circle':
                circle = vtk.vtkDiskSource()
                circle.SetOuterRadius(attrs[0]*1.005)
                circle.SetInnerRadius(attrs[0])
                circle.SetRadialResolution(90)
                circle.SetCircumferentialResolution(90)

                source = circle

            elif primitive == 'Box2d':
                source = vtk.vtkCubeSource()
                source.SetXLength(attrs[0])
                source.SetYLength(attrs[1])
                source.SetZLength(self.opts.depth_2d)

            elif primitive == 'Segment':
                line = vtk.vtkLineSource()
                (x1, y1, x2, y2) = attrs

                line.SetPoint1(x1, y1, 0)
                line.SetPoint2(x2, y2, 0)

                source = vtk.vtkTubeFilter()
                source.SetInputConnection(line.GetOutputPort())
                source.SetRadius(.1)
                source.SetNumberOfSides(30)
                source.Update()

            elif primitive == 'Line':
                line = vtk.vtkLineSource()
                (a, b, c) = attrs
                a2pb2 = a*a + b*b
                assert (a2pb2 > 0)

                if (b != 0):
                    x0 = 0.
                    y0 = -c / b
                else:
                    x0 =-c / a
                    y0 = 0.

                x1 = x0 - b
                y1 = y0 + a

                x2 = x0 + b
                y2 = y0 - a

                line.SetPoint1(x1, y1, 0)
                line.SetPoint2(x2, y2, 0)

                source = vtk.vtkTubeFilter()
                source.SetInputConnection(line.GetOutputPort())
                source.SetRadius(.1)
                source.SetNumberOfSides(30)
                source.Update()

            self.readers[shape_name] = source
            mapper = vtk.vtkCompositePolyDataMapper()
            if not self.opts.imr:
                mapper.ImmediateModeRenderingOff()
            mapper.SetInputConnection(source.GetOutputPort())
            self.mappers[shape_name] = (x for x in [mapper])
            if self.opts.with_edges:
                mapper_edge = vtk.vtkCompositePolyDataMapper()
                if not self.opts.imr:
                    mapper_edge.ImmediateModeRenderingOff()
                    mapper_edge.SetInputConnection(source.GetOutputPort())
                self.mappers_edges[shape_name] = (y for y in [mapper_edge])



    def init_shapes(self):
        self.print_verbose_level(1,'init_shapes')
        for shape_name in self.io.shapes():
            self.init_shape(shape_name)

        for shape_name in self.mappers.keys():
            if shape_name not in self.unfrozen_mappers:
                self.unfrozen_mappers[shape_name] = next(self.mappers[shape_name])
        if self.opts.with_edges:
            for shape_name in self.mappers_edges.keys():
                if shape_name not in self.unfrozen_mappers_edges:
                    self.unfrozen_mappers_edges[shape_name] = next(self.mappers_edges[shape_name])


    def init_contactor(self, contactor_instance_name, instance, instid):
        self.print_verbose_level(2,'init_contactor',  contactor_instance_name)
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

        if not (self.opts.global_filter or self.opts.export):
            actor = vtk.vtkActor()
            if self.opts.with_edges:
                actor_edge = vtk.vtkActor()

            if instance.attrs.get('mass', 0) > 0:
                # objects that may move
                self.dynamic_actors[instid].append((actor, contact_shape_indx,
                                                    collision_group))
                actor.GetProperty().SetOpacity(
                    self.config.get('dynamic_opacity', 0.7))

                actor.GetProperty().SetColor(
                    instance.attrs.get('color',
                                       self.config.get('dynamic_bodies_color', [0.3,0.3,0.3])))
                if self.opts.with_edges:
                    self.dynamic_actors[instid].append((actor_edge, contact_shape_indx,
                                                    collision_group))
                    actor_edge.GetProperty().SetOpacity(
                        self.config.get('dynamic_opacity', 1.0))
                    actor_edge.GetProperty().SetRepresentationToWireframe()

            else:
                # objects that are not supposed to move
                self.static_actors[instid].append((actor, contact_shape_indx,
                                                   collision_group))
                actor.GetProperty().SetOpacity(
                    self.config.get('static_opacity', 1.0))
                actor.GetProperty().SetColor(
                    instance.attrs.get('color',
                                       self.config.get('static_bodies_color', [0.5,0.5,0.5])))

            if self.opts.with_random_color :
                actor.GetProperty().SetColor(random_color())
                if self.opts.with_edges:
                    actor_edge.GetProperty().SetColor(random_color())

            actor.SetMapper(self.unfrozen_mappers[contact_shape_indx])
            if self.opts.with_edges:
                actor_edge.SetMapper(self.unfrozen_mappers_edges[contact_shape_indx])

            if not (self.opts.global_filter or self.opts.export):
                self.renderer.AddActor(actor)
                if self.opts.with_edges:
                    self.renderer.AddActor(actor_edge)

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
            if not (self.opts.global_filter or self.opts.export):
                actor.SetUserTransform(scale_transform)
                if self.opts.with_edges:
                    actor_edge.SetUserTransform(scale_transform)
        else:
            transformer.SetTransform(transform)
            if not (self.opts.global_filter or self.opts.export):
                actor.SetUserTransform(transform)
                if self.opts.with_edges:
                    actor_edge.SetUserTransform(transform)

        self.transformers[contact_shape_indx] = transformer

        self.transforms[instid].append(transform)

        if 'center_of_mass' in instance.attrs:
            center_of_mass = instance.\
                             attrs['center_of_mass'].astype(float)
        else:
            center_of_mass = [0., 0., 0.]

        offset_orientation= contactor.attrs['orientation'].astype(float)

        # for disk, we change the offset since cylinder source are directed along the y axis by default
        # since the disk shapemis invariant with respect to the rotation w.r.t to z-axis
        # we propose to erase it.
        try:
            if self.io.shapes()[contact_shape_name].attrs['primitive'] == 'Disk':
                offset_orientation = [math.cos(pi/4.0), math.sin(pi/4.0), 0., 0.]
        except:
            pass
        self.offsets[instid].append(
            (numpy.subtract(contactor.attrs['translation'].astype(float),
                            center_of_mass),
             offset_orientation))

        cc = CellConnector(
            instid,
            data_names=['instance', 'translation',
                        'velocity(linear)','velocity(angular)',
                        'kinetic_energy', 'displacement'],
            data_sizes=[1, 3, 3, 3, 1, 3])

        if self.cell_connectors.get(instid, None) == None :
            self.cell_connectors[instid] = [ cc ]
        else:
            self.cell_connectors[instid].append(cc)


        cc.SetInputConnection(
            transformer.GetOutputPort())

        self.objects_collector.AddInputConnection(
            cc.GetOutputPort())
        cc.Update()


    def init_instance(self, instance_name):
        self.print_verbose_level(2,'init_instance',  instance_name)
        instance = self.io.instances()[instance_name]
        instid = int(instance.attrs['id'])
        self.transforms[instid] = []
        self.offsets[instid] = []
        if 'time_of_birth' in instance.attrs:
            self.times_of_birth[instid] = instance.attrs['time_of_birth']
        if 'time_of_death' in instance.attrs:
            self.times_of_death[instid] = instance.attrs['time_of_death']
        if 'mass' in instance.attrs:
            # a dynamic instance
            self.mass[instid] = instance.attrs['mass']
            if 'inertia' in instance.attrs:
                inertia = instance.attrs['inertia']
                if self.io.dimension() ==3 :
                    if len(inertia.shape) > 1 and inertia.shape[0] == inertia.shape[1] == 3:
                        self.inertia[instid] = inertia
                    else:
                        self.inertia[instid] = numpy.zeros((3, 3))
                        self.inertia[instid][0, 0] = inertia[0]
                        self.inertia[instid][1, 1] = inertia[1]
                        self.inertia[instid][2, 2] = inertia[2]
                elif self.io.dimension() ==2 :
                     self.inertia[instid] = inertia
            else:
                if self.io.dimension() ==3 :
                    self.inertia[instid] = numpy.eye(3)
                elif self.io.dimension() ==2 :
                    self.inertia[instid] = 1.0

        else:
            pass

        if instid >= 0:
            self.dynamic_actors[instid] = list()
        else:
            self.static_actors[instid] = list()

        for contactor_instance_name in instance:
            self.init_contactor(contactor_instance_name, instance,  instid)

    def init_instances(self):
        self.print_verbose_level(1,'init_instances')
        for instance_name in self.io.instances():
            self.init_instance(instance_name)

    # this sets the position for all transforms associated to an instance
    def set_position_i(self, instance, q0, q1, q2, q3, q4, q5, q6):
        # all objects are set to a nan position at startup,
        # so they are invisibles


        if (numpy.any(numpy.isnan([q0, q1, q2, q3, q4, q5, q6]))
            or numpy.any(numpy.isinf([q0, q1, q2, q3, q4, q5, q6]))):
            print('Bad position for object number', int(instance),' :',  q0, q1, q2, q3, q4, q5, q6)
        else:
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
        self.print_verbose_level(2,'set_position')
        self.set_position_v(data[:, 1],
                            data[:, 2],
                            data[:, 3],
                            data[:, 4],
                            data[:, 5],
                            data[:, 6],
                            data[:, 7],
                            data[:, 8])

    def build_set_functions(self, cc=None):
        if cc is None: cc = self.cell_connectors

        # the numpy vectorization is ok on column vectors for each args
        self.set_position_v = numpy.vectorize(self.set_position_i)

        # here the numpy vectorization is used with a column vector and a
        # scalar for the time arg
        self.set_visibility_v = numpy.vectorize(self.set_dynamic_instance_visibility)
        if self.static_actors:
            self.set_visibility_static_v = numpy.vectorize(self.set_static_instance_visibility)
        else:
            self.set_visibility_static_v = None

        def set_velocity(instance, v0, v1, v2, v3, v4, v5):
            if instance in cc:
                for ccc in cc[instance]:
                    ccc._datas[2][:] = [v0, v1, v2]
                    ccc._datas[3][:] = [v3, v4, v5]
                    if self.io.dimension() == 3 :
                        ccc._datas[4][:] = \
                            0.5*(self.mass[instance]*(v0*v0+v1*v1+v2*v2) +
                                 numpy.dot([v3, v4, v5],
                                           numpy.dot(self.inertia[instance],
                                                     [v3, v4, v5])))
                    elif self.io.dimension() == 2 :
                        kinetic_energy= 0.5*(self.mass[instance]*(v0*v0+v1*v1+v2*v2) +
                                             self.inertia[instance]*(v5*v5))
                        #print('velo', instance, v0, v1, v2, v3, v4, v5)
                        #print('mass', self.mass[instance])
                        #print('inertia', self.inertia[instance])
                        #print('kinetic_energy', kinetic_energy)
                        ccc._datas[4][:] = kinetic_energy



        self.set_velocity_v = numpy.vectorize(set_velocity)

        def set_translation(instance, x0, x1, x2, x0_0, x1_0, x2_0 ):
            if instance in cc:
                for ccc in cc[instance]:
                    ccc._datas[1][:] = [x0, x1, x2]
                    ccc._datas[5][:] = [x0-x0_0, x1-x1_0, x2-x2_0]

        self.set_translation_v = numpy.vectorize(set_translation)

        def set_instance(instance):
            if instance in cc:
                for ccc in cc[instance]:
                    ccc._datas[0][:] = [instance]

        self.set_instance_v = numpy.vectorize(set_instance)

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

    # set visibility for all actors associated to a static instance
    def set_static_instance_visibility(self, instance, time):
        tob = self.times_of_birth.get(instance, -1)
        tod = self.times_of_death.get(instance, infinity)
        has_avatar = False
        if self.opts.visible_mode=='avatars' or self.opts.visible_mode=='contactors':
            for actor, index, group in self.static_actors[instance]:
                if group==-1:
                    has_avatar = True
                    break
        if (tob <= time and tod >= time):
            for actor, index, group in self.static_actors[instance]:
                if not has_avatar or visible_mode=='all':
                    actor.VisibilityOn()
                elif visible_mode=='avatars' and group==-1 and has_avatar:
                    actor.VisibilityOn()
                elif visible_mode=='contactors' and group!=-1 and has_avatar:
                    actor.VisibilityOn()
                else:
                    actor.VisibilityOff()
        else:
            for actor, index, group in self.static_actors[instance]:
                actor.VisibilityOff()

    def set_dynamic_actors_visibility(self, time):
        self.set_visibility_v(list(self.dynamic_actors.keys()), time)

    def set_static_actors_visibility(self, time):
        if self.set_visibility_static_v:
            self.set_visibility_static_v(list(self.static_actors.keys()), time)

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

        if self.opts.export:
            # For time_of_birth specifications with export mode:
            # a 0 scale for objects whose existence is deferred.
            # The correct transform will be set in set_position when
            # the objects appears in pos_data.
            # One have to disable vtkMath generic warnings in order to avoid
            # plenty of 'Unable to factor linear system'
            vtk.vtkMath.GlobalWarningDisplayOff()
            for instance_name in self.io.instances():
                instance = self.io.instances()[instance_name]
                instid = int(instance.attrs['id'])
                for transform in self.transforms[instid]:
                    transform.Scale(0, 0, 0)
        self.time0 = None
        if len(self.io_reader._times) > 0:
            # Positions at first time step
            self.time0 = self.io_reader._times[0]
            self.io_reader.SetTime(self.time0)
            #self.pos_t0 = dsa.WrapDataObject(self.io_reader.GetOutputDataObject(0).GetFieldData().GetArray('pos_data'))

            self.pos_t0 = [self.io_reader.pos_data]

        else:
            # this is for the case simulation has not been ran and
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

        if numpy.shape(self.io_reader._spos_data)[0] > 0:
            self.set_position(self.io_reader._spos_data)
        #    print('self.io_reader._spos_data', self.io_reader._spos_data)
            # static objects are always visible
            #for instance, actors in self.static_actors.items():
            #    for actor,_,_ in actors:
            #         actor.VisibilityOn()

        self.set_position(*self.pos_t0)

        self.set_static_actors_visibility(self.time0)
        self.set_dynamic_actors_visibility(self.time0)

    def setup_vtk_renderer(self):
        self.print_verbose_level(1,'setup_tvk_renderer')
        self.renderer_window.AddRenderer(self.renderer)
        self.interactor_renderer.SetRenderWindow(self.renderer_window)
        self.interactor_renderer.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

        # http://www.itk.org/Wiki/VTK/Depth_Peeling

        if self.opts.depth_peeling:
            # Use a render window with alpha bits (as initial value is 0 (false) ):
            self.renderer_window.SetAlphaBitPlanes(1)

            # Force to not pick a framebuffer with a multisample buffer ( as initial
            # value is 8):
            self.renderer_window.SetMultiSamples(0)

            # Choose to use depth peeling (if supported) (initial value is 0
            # (false) )
            self.renderer.SetUseDepthPeeling(1)

            # Set depth peeling parameters.
            self.renderer.SetMaximumNumberOfPeels(
                self.opts.maximum_number_of_peels)

            # Set the occlusion ratio (initial value is 0.0, exact image)
            self.renderer.SetOcclusionRatio(self.opts.occlusion_ratio)

        # Set the initial camera position and orientation if specified
        if self.opts.initial_camera[0] is not None:
            self.renderer.GetActiveCamera().SetPosition(*self.opts.initial_camera[0])
        if self.opts.initial_camera[1] is not None:
            self.renderer.GetActiveCamera().SetFocalPoint(*self.opts.initial_camera[1])
        if self.opts.initial_camera[2] is not None:
            self.renderer.GetActiveCamera().SetViewUp(*self.opts.initial_camera[2])
        if self.opts.initial_camera[4] is not None:
            self.renderer.GetActiveCamera().SetClippingRange(*self.opts.initial_camera[4])
        else:
            self.renderer.ResetCameraClippingRange()
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

        self.renderer_window.SetSize(*self.config['window_size'])
        self.renderer_window.SetWindowName('vview: ' + self.opts.io_filename)

    def setup_charts(self):
        self.print_verbose_level(1,'setup_charts')
        # Warning! numpy support offer a view on numpy array
        # the numpy array must not be garbage collected!
        nxtime = self.io_reader._isolv_data[:, 0]
        nxiters = self.io_reader._isolv_data[:, 1]
        nprecs = self.io_reader._isolv_data[:, 2]
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
        self.print_verbose_level(1,'setup_sliders')
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

        self.input_observer.set_opacity()

        self.interactor_renderer.AddObserver('KeyPressEvent', self.input_observer.key)

        self.interactor_renderer.AddObserver(
            'TimerEvent', self.input_observer.recorder_observer)

        if self.io.contact_forces_data().shape[0] > 0:
            self.slwsc, self.slrepsc = self.make_slider(
                'CF scale',
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
        self.print_verbose_level(1,'setup_axes')
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

    def setup_siconos_text(self):
        self.print_verbose_level(1,'setup_siconos_text')
        # Create a text actor.
        txt = vtkTextActor()
        txt.SetInput('Siconos')
        txtprop = txt.GetTextProperty()
        txtprop.SetFontFamilyToArial()
        txtprop.BoldOn()
        txtprop.SetFontSize(88)
        txtprop.ShadowOn()
        txtprop.SetShadowOffset(4, 4)
        background_color = self.config.get('background_color', [.0,.0,.0])
        reverse_background_color =numpy.ones(3) - background_color

        txtprop.SetColor(*reverse_background_color)
        txt.SetDisplayPosition(20, 30)

        # Assign actor to the renderer.
        self.renderer.AddActor(txt)
        #self._renderer.SetBackground(colors.GetColor3d('DarkGreen'))



    # this should be extracted from the VView class
    def export(self):
        self.print_verbose('export start')
        times = self.io_reader._times[
                self.opts.start_step:self.opts.end_step:self.opts.stride]
        ntime = len(times)

        if self.opts.gen_para_script:
            # just the generation of a parallel command
            options_str = ''
            if self.opts.ascii_mode:
                options_str += '--ascii '
            if self.opts.global_filter:
                options_str += '--global-filter '
            if self.opts.depth_2d != 0.1:
                options_str += ' --depth-2d='+str(self.opts.depth_2d)

            ntimes_proc = int(ntime / self.opts.nprocs)
            s = ''
            for i in range(self.opts.nprocs):
                s += '{0}/{1} '.format(ntimes_proc*i, ntimes_proc*(i+1))
            print('#!/bin/sh')
            print('parallel --verbose', sys.argv[0], self.opts.io_filename,
                  options_str, '--start-step={//} --end-step={/} :::', s)

        else:
            # export
            big_data_writer = vtk.vtkXMLMultiBlockDataWriter()
            add_compatiblity_methods(big_data_writer)

            big_data_writer.SetInputConnection(self.big_data_collector.GetOutputPort())

            if self.opts.ascii_mode:
                big_data_writer.SetDataModeToAscii()
            k = self.opts.start_step

            packet = int(ntime/100)+1

            initial_time =True

            for time in times:
                k = k + self.opts.stride
                if (k % packet == 0):
                    sys.stdout.write('.')
                    sys.stdout.flush()

                self.io_reader.SetTime(time)

                spos_data = self.io_reader.pos_static_data
                #print('spos_data at time 1' , time, spos_data)
                if spos_data.size > 0:
                    self.set_position_v(spos_data[:, 1], spos_data[:, 2],
                                        spos_data[:, 3],
                                        spos_data[:, 4], spos_data[:, 5],
                                        spos_data[:, 6],
                                        spos_data[:, 7], spos_data[:, 8])


                pos_data = self.io_reader.pos_data
                velo_data = self.io_reader.velo_data

                if initial_time:
                    initial_pos_data = self.io_reader.pos_data
                    initial_velo_data = self.io_reader.velo_data
                    initial_time = False


                self.set_position_v(
                    pos_data[:, 1], pos_data[:, 2], pos_data[:, 3],
                    pos_data[:, 4], pos_data[:, 5], pos_data[:, 6],
                    pos_data[:, 7], pos_data[:, 8])

                self.set_velocity_v(
                    velo_data[:, 1],
                    velo_data[:, 2],
                    velo_data[:, 3],
                    velo_data[:, 4],
                    velo_data[:, 5],
                    velo_data[:, 6],
                    velo_data[:, 7])

                self.set_translation_v(
                    pos_data[:, 1],
                    pos_data[:, 2],
                    pos_data[:, 3],
                    pos_data[:, 4],
                    initial_pos_data[:, 2],
                    initial_pos_data[:, 3],
                    initial_pos_data[:, 4]
                )

                self.set_instance_v(pos_data[:, 1])

                big_data_writer.SetFileName('{0}-{1}.{2}'.format(
                    os.path.splitext(os.path.basename(self.opts.io_filename))[0],
                    k, big_data_writer.GetDefaultFileExtension()))
                if self.opts.global_filter:
                    self.big_data_geometry_filter.Update()
                else:
                    self.big_data_collector.Update()

                big_data_writer.Write()

            big_data_writer.Write()

    def export_raw_data(self):

        times = self.io_reader._times[
                self.opts.start_step:self.opts.end_step:self.opts.stride]
        ntime = len(times)
        export_2d = False
        if self.io.dimension() ==2 :
            export_2d=True
            print('We export raw data for 2D object')
        # export
        k = self.opts.start_step
        packet = int(ntime/100)+1

        # ######## position output ########

        # nvalue = ndyna*7+1

        # position_output = numpy.empty((ntime,nvalue))
        # #print('position_output shape', numpy.shape(position_output))
        # position_output[:,0] = times[:]
        position_output = {}
        velocity_output = {}
        velocity_absolute_output = {}

        for time in times:
            k = k + self.opts.stride
            if (k % packet == 0):
                sys.stdout.write('.')
            self.io_reader.SetTime(time)

            pos_data = self.io_reader.pos_data
            velo_data = self.io_reader.velo_data

            ndyna=pos_data.shape[0]


            for i in range(ndyna):
                bdy_id = int(pos_data[i,1])

                ######## position output ########
                if self.opts._export_position :
                    nvalue=pos_data.shape[1]
                    position_output_bdy = position_output.get(bdy_id)
                    if position_output_bdy is None:
                        position_output[bdy_id] = []
                    position_output_body =  position_output[bdy_id]
                    position_output_body.append([])
                    position_output_body[-1].append(time)
                    if export_2d:
                        data_2d =  [pos_data[i,2],pos_data[i,3],numpy.acos(pos_data[i,5]/2.0)]
                        position_output_body[-1].extend(data_2d)
                        #position_output_body[-1].extend(pos_data[i,2:nvalue])
                    else:
                        position_output_body[-1].extend(pos_data[i,2:nvalue])
                    position_output_body[-1].append(bdy_id)


                ######## velocity output ########
                if self.opts._export_velocity :
                    nvalue=velo_data.shape[1]
                    velocity_output_bdy = velocity_output.get(bdy_id)
                    if velocity_output_bdy is None:
                        velocity_output[bdy_id] = []
                    velocity_output_body =  velocity_output[bdy_id]
                    velocity_output_body.append([])
                    velocity_output_body[-1].append(time)
                    velocity_output_body[-1].extend(velo_data[i,2:nvalue])
                    velocity_output_body[-1].append(bdy_id)

                ######## velocity in absolute frame output ########
                if self.opts._export_velocity_in_absolute_frame :
                    nvalue=velo_data.shape[1]
                    [q1,q2,q3,q4] = pos_data[i,5:9]
                    q = Quaternion((q1, q2, q3, q4))
                    velo = q.rotate(velo_data[i,5:8])
                    velocity_absolute_output_bdy = velocity_absolute_output.get(bdy_id)
                    if velocity_absolute_output_bdy is None:
                        velocity_absolute_output[bdy_id] = []
                    velocity_absolute_output_body =  velocity_absolute_output[bdy_id]
                    velocity_absolute_output_body.append([])
                    velocity_absolute_output_body[-1].append(time)
                    velocity_absolute_output_body[-1].extend(velo_data[i,2:5])
                    velocity_absolute_output_body[-1].extend(velo[:])
                    velocity_absolute_output_body[-1].append(bdy_id)

        for bdy_id in position_output.keys():
            output = numpy.array(position_output[bdy_id])
            filename_output = '{0}-position-body_{1}.dat'.format(
                os.path.splitext(os.path.basename(self.opts.io_filename))[0],
            bdy_id)
            numpy.savetxt(filename_output, output)

        for bdy_id in velocity_output.keys():
            output = numpy.array(velocity_output[bdy_id])
            filename_output = '{0}-velocity-body_{1}.dat'.format(
                os.path.splitext(os.path.basename(self.opts.io_filename))[0],
            bdy_id)
            numpy.savetxt(filename_output, output)

        for bdy_id in velocity_absolute_output.keys():
            output = numpy.array(velocity_absolute_output[bdy_id])
            filename_output = '{0}-velocity-absolute-body_{1}.dat'.format(
                os.path.splitext(os.path.basename(self.opts.io_filename))[0],
            bdy_id)
            numpy.savetxt(filename_output, output)

        cf_output = {}
        for time in times:
            #print('time', time)
            k = k + self.opts.stride
            if (k % packet == 0):
                sys.stdout.write('.')

            self.io_reader.SetTime(time)

            cf_data = self.io_reader.cf_data
            #print('cf_data', cf_data)

            if cf_data is not None and self.opts._export_cf :
                ncontact=cf_data.shape[0]

                for i in range(ncontact):
                    contact_id = int(cf_data[i,23])
                    #print('contact_id', contact_id)
                    ######## contact output ########
                    nvalue=cf_data.shape[1]
                    cf_output_contact = cf_output.get(contact_id)
                    if cf_output_contact is None:
                        cf_output[contact_id] = []
                    cf_output_contact =  cf_output[contact_id]
                    cf_output_contact.append([])
                    cf_output_contact[-1].append(time)
                    cf_output_contact[-1].extend(cf_data[i,2:nvalue])
                    cf_output_contact[-1].append(contact_id)

        for contact_id in cf_output.keys():
            output = numpy.array(cf_output[contact_id])
            filename_output = '{0}-cf-contact_{1}.dat'.format(
                os.path.splitext(os.path.basename(self.opts.io_filename))[0],
            contact_id)
            numpy.savetxt(filename_output, output)


        sys.stdout.write('\n')


    def initialize_vtk(self):
        self.print_verbose('initialize_vtk')
        if not self.opts.gen_para_script:

            self.objects_collector = vtk.vtkMultiBlockDataGroupFilter()
            add_compatiblity_methods(self.objects_collector)
            self.cf_collector = vtk.vtkMultiBlockDataGroupFilter()
            add_compatiblity_methods(self.cf_collector)

            self.big_data_collector = vtk.vtkMultiBlockDataGroupFilter()
            add_compatiblity_methods(self.big_data_collector)

            self.big_data_collector.AddInputConnection(self.cf_collector.GetOutputPort())

            if self.opts.global_filter:

                self.big_data_geometry_filter = vtk.vtkCompositeDataGeometryFilter()
                add_compatiblity_methods(self.big_data_geometry_filter)
                self.big_data_geometry_filter.SetInputConnection(self.objects_collector.GetOutputPort())

                self.big_data_collector.AddInputConnection(self.big_data_geometry_filter.GetOutputPort())
            else:
                self.big_data_collector.AddInputConnection(self.objects_collector.GetOutputPort())


            if self.opts.global_filter and not self.opts.export:
                self.big_data_mapper = vtk.vtkCompositePolyDataMapper()
                add_compatiblity_methods(self.big_data_mapper)
                self.big_data_mapper.SetInputConnection(self.big_data_collector.GetOutputPort())

                if not self.opts.imr:
                    self.big_data_mapper.ImmediateModeRenderingOff()

                self.big_actor = vtk.vtkActor()
                self.big_actor.SetMapper(self.big_data_mapper)

            times = self.io_reader._times

            if (len(times) == 0):
                print('No dynamic data found!  Empty simulation.')

            self.readers = dict()
            self.datasets = dict()
            self.mappers = dict()
            self.mappers_edges = dict()
            self.dynamic_actors = dict()
            self.static_actors = dict()
            self.vtk_reader = {'vtp': vtk.vtkXMLPolyDataReader,
                               'stl': vtk.vtkSTLReader}
            self.unfrozen_mappers = dict()
            self.unfrozen_mappers_edges = dict()

            self.build_set_functions()

            self.renderer = vtk.vtkRenderer()
            self.renderer_window = vtk.vtkRenderWindow()
            self.interactor_renderer = vtk.vtkRenderWindowInteractor()

            self.init_shapes()
            self.init_instances()

            if self.opts.cf_disable:
                self.io_reader.ContactForcesOff()
                self.setup_initial_position()
            else:
                self.io_reader.ContactForcesOn()
                for mu in self.io_reader._mu_coefs:
                    self.init_contact_pos(mu)

                self.setup_initial_position()
                transform = vtk.vtkTransform()
                transform.Translate(-0.5, 0., 0.)
                for mu in self.io_reader._mu_coefs:
                    self.init_cf_sources(mu, transform)

            if not self.opts.export:
                if not self.opts.cf_disable and not self.opts.global_filter:
                    for mu in self.io_reader._mu_coefs:
                        self.renderer.AddActor(self.gactor[mu])
                        self.renderer.AddActor(self.cactor[mu])

                        self.renderer.AddActor(self.clactor[mu])
                        self.renderer.AddActor(self.sactora[mu])
                        self.renderer.AddActor(self.sactorb[mu])

                if self.opts.global_filter:
                    self.renderer.AddActor(self.big_actor)

                self.renderer.ResetCamera()

    def initialize_gui(self):
        self.print_verbose('initialize_gui')
        self.setup_vtk_renderer()
        self.setup_sliders(self.io_reader._times)
        if self.opts.with_charts:
            self.setup_charts()
        self.setup_axes()
        self.setup_siconos_text()

        self.gui_initialized = True

    def run(self):
        self.print_verbose('vview start')
        self.initialize_vtk()
        self.initialize_gui()
        self.interactor_renderer.Start()
        self.print_verbose('vview end')
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

nan = numpy.nan
if __name__=='__main__':
    ## Options and config already loaded above
    with MechanicsHdf5(io_filename=opts.io_filename, mode='r') as io:
        vview = VView(io, opts, config)
        vview.run()

    # Update configuration and save it
    config['window_size'] = vview.renderer_window.GetSize()
    config.save_configuration(force=False)
