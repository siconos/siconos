#!/usr/bin/env @PYTHON_EXECUTABLE@

import OCC

from OCC import VERSION

from OCC.gp import gp_Ax1, gp_Pnt, gp_Dir, gp_Trsf, gp_Quaternion,gp_Vec, gp_XYZ
from OCC.TopLoc import TopLoc_Location

#from OCC.Display.SimpleGui import get_backend

from OCC.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCC.Interface import Interface_Static_SetCVal

from OCC.BRep import BRep_Builder
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.TopoDS import TopoDS_Compound

from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity

import OCC.Graphic3d as Graphic3d
from OCC.Quantity import Quantity_NOC_DARKVIOLET, Quantity_NOC_BLUE1, Quantity_NOC_GREEN, Quantity_NOC_RED, Quantity_NOC_ORANGE, Quantity_NOC_SALMON, Quantity_NOC_YELLOW

import vtk
from vtk.util import numpy_support

import sys

#from Siconos.Mechanics import IO
import siconos.io.mechanics_io as IO

from siconos.io.SimpleGui import get_backend,init_display

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
        #assert(rv[0] == 0)
        return [rv[1], rv[2], rv[3]]

    def axisAngle(self):
        r = [0,0,0]
        a = self._data.GetRotationAngleAndAxis(r)
        return r, a






def memoize(f):
    """ Memoization decorator for a function taking one or more arguments. """
    class memodict(dict):
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            ret = self[key] = f(*key)
            return ret

    return memodict().__getitem__


def make_slider(minv, maxv, vstep):
    from PyQt4 import QtCore, QtGui, QtOpenGL
    class SlidersGroup(QtGui.QGroupBox):

        valueChanged = QtCore.pyqtSignal(int)

        def __init__(self, orientation, title, parent=None):
            super(SlidersGroup, self).__init__(title, parent)

            self.slider = QtGui.QSlider(orientation)
            self.slider.setFocusPolicy(QtCore.Qt.StrongFocus)
            self.slider.setTickPosition(QtGui.QSlider.TicksBothSides)
            self.slider.setTickInterval(10)
            self.slider.setSingleStep(1)

            self.slider.valueChanged.connect(self.setValue)
            self.slider.valueChanged.connect(self.valueChanged)

            if orientation == QtCore.Qt.Horizontal:
                direction = QtGui.QBoxLayout.TopToBottom
            else:
                direction = QtGui.QBoxLayout.LeftToRight

            slidersLayout = QtGui.QBoxLayout(direction)
            slidersLayout.addWidget(self.slider)
            self.setLayout(slidersLayout)

        def setValue(self, value):
            vstep(value)
            self.slider.setValue(value)

        def setMinimum(self, value):
            self.slider.setMinimum(value)

        def setMaximum(self, value):
            self.slider.setMaximum(value)


    class SliderWindow(QtGui.QWidget):
        def __init__(self):
            super(SliderWindow, self).__init__()

            self.horizontalSliders = SlidersGroup(QtCore.Qt.Horizontal,
                    "Steps")

            self.stackedWidget = QtGui.QStackedWidget()
            self.stackedWidget.addWidget(self.horizontalSliders)

            self.createControls("Controls")

            self.valueSpinBox.valueChanged.connect(self.horizontalSliders.setValue)
            self.horizontalSliders.slider.valueChanged.connect(self.valueSpinBox.setValue)

            layout = QtGui.QHBoxLayout()
            layout.addWidget(self.controlsGroup)
            layout.addWidget(self.stackedWidget)
            self.setLayout(layout)

            self.minimumSpinBox.setValue(minv)
            self.maximumSpinBox.setValue(maxv)
            self.valueSpinBox.setValue(minv)

            self.setWindowTitle("Step")

        def createControls(self, title):
            self.controlsGroup = QtGui.QGroupBox(title)

            minimumLabel = QtGui.QLabel("Minimum step:")
            maximumLabel = QtGui.QLabel("Maximum step:")
            valueLabel = QtGui.QLabel("Current step:")

            self.minimumSpinBox = QtGui.QSpinBox()
            self.minimumSpinBox.setRange(minv, maxv)
            self.minimumSpinBox.setSingleStep(1)

            self.maximumSpinBox = QtGui.QSpinBox()
            self.maximumSpinBox.setRange(minv, maxv)
            self.maximumSpinBox.setSingleStep(1)

            self.valueSpinBox = QtGui.QSpinBox()
            self.valueSpinBox.setRange(minv, maxv)
            self.valueSpinBox.setSingleStep(1)

            self.minimumSpinBox.valueChanged.connect(self.horizontalSliders.setMinimum)
            self.maximumSpinBox.valueChanged.connect(self.horizontalSliders.setMaximum)

            controlsLayout = QtGui.QGridLayout()
            controlsLayout.addWidget(minimumLabel, 0, 0)
            controlsLayout.addWidget(maximumLabel, 1, 0)
            controlsLayout.addWidget(valueLabel, 2, 0)
            controlsLayout.addWidget(self.minimumSpinBox, 0, 1)
            controlsLayout.addWidget(self.maximumSpinBox, 1, 1)
            controlsLayout.addWidget(self.valueSpinBox, 2, 1)

            self.controlsGroup.setLayout(controlsLayout)

    return SliderWindow()



with IO.Hdf5('siconos-mechanisms.hdf5', 'r') as io:

    display, start_display, add_menu, add_function_to_menu, win, app = init_display()

    dpos_data = io.dynamic_data()[:]
    nbobjs = len(filter(lambda x: io.instances()[x].attrs['id'] >= 0, io.instances()))

    nbsteps = dpos_data.shape[0] / nbobjs

    assert nbsteps * nbobjs == dpos_data.shape[0]

    current_color = 0
    @memoize
    def make_shape(shape_name):
        global current_color

        # cf CADMBTB_API, but cannot get the same color order
        colors = list(reversed([Quantity_NOC_DARKVIOLET, Quantity_NOC_BLUE1, Quantity_NOC_GREEN, Quantity_NOC_RED, Quantity_NOC_ORANGE, Quantity_NOC_SALMON, Quantity_NOC_YELLOW]))

        with IO.tmpfile(contents=io.shapes()[shape_name][:][0]) as tmpfile:

            step_reader = STEPControl_Reader()

            status = step_reader.ReadFile(tmpfile[1])

            if status == IFSelect_RetDone:  # check status
                failsonly = False
                step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
                step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)

                ok = step_reader.TransferRoot(1)
                nbs = step_reader.NbShapes()

                l=[]
                for i in range(1, nbs+1):
                    ais_shape = display.DisplayShape(step_reader.Shape(i), update=True, transparency=.55)
                    ais_shape.GetObject().SetColor(colors[current_color % 6])
                    current_color += 1
                    ais_shape.GetObject().SetMaterial(Graphic3d.Graphic3d_NOM_PLASTIC)
                    l.append(ais_shape)

                return l


    obj_by_id = dict()
    for instance in io.instances():
        obj_by_id[io.instances()[instance].attrs['id']] = instance

    def get_offset(instance_name, shape_name):
        return (io.instances()[instance_name][shape_name].attrs['translation'],
                io.instances()[instance_name][shape_name].attrs['orientation'])

    def shape_names(obj):
        return [io.instances()[obj][shape].attrs['name']
                for shape in io.instances()[obj]]

    @memoize
    def avatars(obj):
        l = [make_shape(io.instances()[obj][shape].attrs['name'])
             for shape in io.instances()[obj]]
        # flatten
        return [item for sublist in l for item in sublist]

    @memoize
    def write_step(stshape):
        # initialize the STEP exporter
        step_writer = STEPControl_Writer()

        # missing person => load failure
        #Interface_Static_SetCVal("write.step.schema", "AP203")

        step_str, shape = stshape

        step_writer.Transfer(shape, STEPControl_AsIs)
        status = step_writer.Write('siconos-mechanisms-{0}.stp'.format(step_str))

    def vstep(step_str):

        step = int(step_str)

        positions = dpos_data[nbobjs*step:nbobjs*step+nbobjs, 2:]

        builder = BRep_Builder()
        comp = TopoDS_Compound()
        builder.MakeCompound(comp)

        for _id in range(positions.shape[0]):

            q0, q1, q2, q3, q4, q5, q6 = [float(x) for x in positions[_id,:]]

            obj = obj_by_id[_id+1]

            q = Quaternion((q3, q4, q5, q6))

            for shape_name, avatar in zip(io.instances()[obj], avatars(obj)):
                offset = get_offset(obj, shape_name)
                p = q.rotate(offset[0])
                r = q*Quaternion(offset[1])

                tr = gp_Trsf()
                qocc = gp_Quaternion(r[1],
                                     r[2],
                                     r[3],
                                     r[0])
                tr.SetRotation(qocc)
                xyz = gp_XYZ(q0 + p[0], q1 + p[1], q2 + p[2])
                vec = gp_Vec(xyz)
                tr.SetTranslationPart(vec)
                loc = TopLoc_Location(tr)

                display.Context.SetLocation(avatar, loc)

                moved_shape = BRepBuilderAPI_Transform(avatar.GetObject().Shape(), tr, True).Shape()

                builder.Add(comp, moved_shape)

            display.Context.UpdateCurrentViewer()

        write_step((step_str, comp))


#    add_menu('run')
#    add_function_to_menu('run', run)

    from PyQt4 import QtGui, QtCore

    sl = make_slider(0, nbsteps, vstep)
    dw = QtGui.QDockWidget()
    win.addDockWidget(QtCore.Qt.DockWidgetArea(QtCore.Qt.TopDockWidgetArea), dw)
    dw.setWidget(sl)

    for instance in io.instances():
        avatars(instance)

    start_display()
