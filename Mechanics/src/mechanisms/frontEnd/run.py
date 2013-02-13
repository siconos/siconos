#!/usr/bin/env python
import mbtb
import cadmbtb
import numpy
import array
import os
import sys

SiconosMechanisms_BUILD=os.environ.get("SiconosMechanisms_BUILD")
SiconosMechanisms_DIR=os.environ.get("SiconosMechanisms_DIR")
my_PI=3.14159265

if SiconosMechanisms_BUILD is None:
   print "Can not find SiconosMechanisms_BUILD from the environement, need it."
   sys.exit()


if SiconosMechanisms_DIR is None:
   print "Can not find SiconosMechanisms_DIR from the environement, need it."
   sys.exit()

execfile(SiconosMechanisms_DIR+"frontEnd/mbtbDefaultOptions.py")
print "run.py: ../mbtbDefaultOptions.py loaded"
 
try:
   execfile("mbtbLocalOptions.py")
   print "run.py: mbtbLocalOptions.py loaded"
except :
   print "run.py, info: mbtbLocalOptions.py not defined"

execfile("bodydef.py")
print "run.py: bodydef.py loaded"

if with3D:
    from OCC.BRepPrimAPI import *
    from OCC.gp import *
    from OCC.TopLoc import *
    from OCC.AIS import *
    from OCC.Display.SimpleGui import *
    display, start_display, add_menu, add_function_to_menu = init_display()
    cadmbtb.CADMBTB_setGraphicContext(display.Context)
    v = display.GetView().GetObject()
    cadmbtb.CADMBTB_setGraphicView(v)
    display.View_Iso()


mbtb.MBTB_init(NBBODIES,NBJOINTS,NBCONTACTS)

if NBCONTACTS>0:
   mbtb.MBTB_ContactSetDParam(2,0,0,ContactArtefactLength)
   mbtb.MBTB_ContactSetDParam(3,0,0,ArtefactThershold)
   mbtb.MBTB_ContactSetDParam(4,0,0,NominalForce)

for idContact in range(NBCONTACTS):
   cadmbtb.CADMBTB_setContactDParam(0,idContact,0,contactTrans1[idContact]) # visu
   cadmbtb.CADMBTB_setContactDParam(0,idContact,1,contactTrans2[idContact]) # visu


#build contacts (must be done before the dynamical system building
for idContact in range(NBCONTACTS):
    mbtb.MBTB_ContactLoadCADFile(idContact,afileContact1[idContact],afileContact2[idContact],contactDraw1[idContact],contactDraw2[idContact])
    mbtb.MBTB_ContactBuild(idContact,contactName[idContact],contactBody1[idContact],contactBody2[idContact],contactType3D[idContact],contactmu[idContact],contacten[idContact],et)
    mbtb.MBTB_ContactSetDParam(1,idContact,0,contactOffset[idContact]) # offset
    mbtb.MBTB_ContactSetIParam(0,idContact,0,contactOffsetP1[idContact]) # application point P1 or Trans P1.
    mbtb.MBTB_ContactSetIParam(1,idContact,0,contactNormalFromFace1[idContact]) # normal computed from Face1.

for idBody in range(NBBODIES):
   cadmbtb.CADMBTB_setShapeDParam(0,idBody,bodyTrans[idBody]) # visu

#build dynamical systems
for idBody in range(NBBODIES):
    mbtb.MBTB_BodyLoadCADFile(idBody,afile[idBody],bodyDraw[idBody])
    mbtb.MBTB_BodyBuild(idBody, body[idBody], m[idBody],
                        initPos[idBody], initCenterMass[idBody],
                        inertialMatrix[idBody],plugin, fctf[idBody],plugin,fctm[idBody])
    mbtb.MBTB_BodySetVelocity(idBody,initVel[idBody])

#build joint
for idJoint in range(NBJOINTS):
    mbtb.MBTB_JointBuild(idJoint,jointName[idJoint],jointType[idJoint],jointBody1[idJoint],jointBody2[idJoint],jointPos[idJoint])

try:
   for idJoint in range(NBJOINTS):
      mbtb.MBTB_setJointPoints(idJoint,jointPoints[2*idJoint],jointPoints[2*idJoint+1])
except NameError:
   print "run.py, info: joint points not defined."

for idArtefact in range(NBARTEFACTS):
   cadmbtb.CADMBTB_loadArtefactCADFile(Artefactfile[idArtefact],ArtefactTrans[idArtefact])
   

mbtb.MBTB_SetDParam(0,TSdeactivateYPosThreshold)
mbtb.MBTB_SetDParam(1,TSdeactivateYVelThreshold)
mbtb.MBTB_SetDParam(2,TSactivateYPosThreshold)
mbtb.MBTB_SetDParam(3,TSactivateYVelThreshold)
mbtb.MBTB_SetDParam(4,TSProjectionMaxIteration)
mbtb.MBTB_SetDParam(5,TSConstraintTol)
mbtb.MBTB_SetDParam(6,TSConstraintTolUnilateral)
mbtb.MBTB_SetDParam(7,TSLevelOfProjection)





mbtb.MBTB_initSimu(stepSize,withProj)
mbtb.MBTB_setGraphicFreq(freqUpdate)
mbtb.MBTB_setOutputFreq(freqOutput)
mbtb.MBTB_setSolverIOption(2,withReduced);
mbtb.MBTB_setSolverIOption(0,solverIt);
mbtb.MBTB_setSolverDOption(0,solverTol);
if gotoPos:
    mbtb.MBTB_moveBodyToPosWithSpeed(0,aPos0,aVel0)
    mbtb.MBTB_moveBodyToPosWithSpeed(1,aPos1,aVel1)
    mbtb.MBTB_moveBodyToPosWithSpeed(2,aPos2,aVel2)

if with3D:
   display.View_Iso()
   display.FitAll()
   

                         # ord('W'): self._display.SetModeWireFrame,
                         # ord('S'): set_shade_mode,
                         # ord('A'): self._display.EnableAntiAliasing,
                         # ord('B'): self._display.DisableAntiAliasing,
                         # ord('Q'): self._display.SetModeQuickHLR,
                         # ord('E'): self._display.SetModeExactHLR,
                         # ord('F'): self._display.FitAll(
ais_boxshp=None

def STEP(event=None):
    mbtb.MBTB_step()
    
    
def STEP_1000(event=None):
   ii=0
   while ii < 1000:
      mbtb.MBTB_step()
      cadmbtb.CADMBTB_DumpGraphic()
      ii=ii+1

    
def RUN10(event=None):
    mbtb.MBTB_run(10)
    
def RUN(event=None):
    mbtb.MBTB_run(stepNumber)

def RUN100(event=None):
    mbtb.MBTB_run(100)

def RUN400(event=None):
    mbtb.MBTB_run(400)


def RUN1000(event=None):
    mbtb.MBTB_run(1000)

def RUN7000(event=None):
    mbtb.MBTB_run(7000)

    
def RUN10000(event=None):
    mbtb.MBTB_run(10000)

def RUN50000(event=None):
    mbtb.MBTB_run(50000)

def RUN300000(event=None):
    mbtb.MBTB_run(300000)

def QUIT(event=None):
    exit();


def GRAPHIC_ON(event=None):
    cadmbtb.CADMBTB_setGraphicContext(display.Context)
    v = display.GetView().GetObject()
    cadmbtb.CADMBTB_setGraphicView(v)


def GRAPHIC_OFF(event=None):
    cadmbtb.CADMBTB_disableGraphic()

def DISABLE_PROJ(event=None):
    mbtb.MBTB_doProj(0)

def ENABLE_PROJ(event=None):
    mbtb.MBTB_doProj(1)

def DISABLE_ONLYPROJ(event=None):
    mbtb.MBTB_doOnlyProj(0)

def ENABLE_ONLYPROJ(event=None):
    mbtb.MBTB_doOnlyProj(1)

def GRAPHIC_DUMP(event=None):
   cadmbtb.CADMBTB_DumpGraphic()

def ENABLE_VERBOSE_MBTB_PRINT_DIST(event=None):
   mbtb.MBTB_print_dist(1)
   cadmbtb.CADMBTB_print_dist(1)
   
def DISABLE_VERBOSE_MBTB_PRINT_DIST(event=None):
   mbtb.MBTB_print_dist(0)
   cadmbtb.CADMBTB_print_dist(0)
 
def ENABLE_VERBOSE_MBTB_DISPLAY_STEP_BODIES(event=None):
   mbtb.MBTB_displayStep_bodies(1)
   
def DISABLE_VERBOSE_MBTB_DISPLAY_STEP_BODIES(event=None):
   mbtb.MBTB_displayStep_bodies(0)

def ENABLE_VERBOSE_MBTB_DISPLAY_STEP_JOINTS(event=None):
   mbtb.MBTB_displayStep_joints(1)
   
def DISABLE_VERBOSE_MBTB_DISPLAY_STEP_JOINTS(event=None):
   mbtb.MBTB_displayStep_joints(0)

def ENABLE_VERBOSE_MBTB_DISPLAY_STEP_CONTACTS(event=None):
   mbtb.MBTB_displayStep_contacts(1)
   
def DISABLE_VERBOSE_MBTB_DISPLAY_STEP_CONTACTS(event=None):
   mbtb.MBTB_displayStep_contacts(0)

def VIEW_TOP(event=None):
   display.View_Top()
def VIEW_BOTTOM(event=None):
   display.View_Bottom()
def VIEW_LEFT(event=None):
   display.View_Left()
def VIEW_RIGHT(event=None):
   display.View_Right()
def VIEW_REAR(event=None):
   display.View_Rear()
def VIEW_FRONT(event=None):
   display.View_Front()
def VIEW_ISO(event=None):
   display.View_Iso()
def VIEW_FITALL(event=None):
   display.FitAll()
    
def VIEW_TOP(event=None):
   display.View_Top()



if with3D and __name__ == '__main__':
   cadmbtb.CADMBTB_setIParam(0,dumpGraphic)
   mbtb.MBTB_ContactSetIParam(2,0,0,drawMode)
   add_menu('MBTB')
   add_menu('OPTION')
   add_menu('VERBOSE')
   add_menu('DISPLAY')
   add_menu('VIEW')
   add_function_to_menu('MBTB', STEP)
   add_function_to_menu('MBTB', RUN)
   add_function_to_menu('MBTB', RUN10)
   add_function_to_menu('MBTB', RUN100)
   add_function_to_menu('MBTB', RUN400)
   add_function_to_menu('MBTB', RUN1000)
   add_function_to_menu('MBTB', RUN7000)
   add_function_to_menu('MBTB', RUN10000)
   add_function_to_menu('MBTB', RUN50000)
   add_function_to_menu('MBTB', RUN300000)
   add_function_to_menu('MBTB', QUIT)
   add_function_to_menu('DISPLAY',GRAPHIC_OFF)
   add_function_to_menu('DISPLAY',GRAPHIC_ON)
   add_function_to_menu('DISPLAY',GRAPHIC_DUMP)
   add_function_to_menu('VERBOSE',ENABLE_VERBOSE_MBTB_PRINT_DIST)
   add_function_to_menu('VERBOSE',DISABLE_VERBOSE_MBTB_PRINT_DIST)
   add_function_to_menu('VERBOSE',ENABLE_VERBOSE_MBTB_DISPLAY_STEP_BODIES)
   add_function_to_menu('VERBOSE',DISABLE_VERBOSE_MBTB_DISPLAY_STEP_BODIES)
   add_function_to_menu('VERBOSE',ENABLE_VERBOSE_MBTB_DISPLAY_STEP_JOINTS)
   add_function_to_menu('VERBOSE',DISABLE_VERBOSE_MBTB_DISPLAY_STEP_JOINTS)
   add_function_to_menu('VERBOSE',ENABLE_VERBOSE_MBTB_DISPLAY_STEP_CONTACTS)
   add_function_to_menu('VERBOSE',DISABLE_VERBOSE_MBTB_DISPLAY_STEP_CONTACTS)
   add_function_to_menu('OPTION',DISABLE_PROJ)
   add_function_to_menu('OPTION',ENABLE_PROJ)
   add_function_to_menu('OPTION',DISABLE_ONLYPROJ)
   add_function_to_menu('OPTION',ENABLE_ONLYPROJ)
   add_function_to_menu('VIEW',VIEW_TOP)
   add_function_to_menu('VIEW',VIEW_BOTTOM)
   add_function_to_menu('VIEW',VIEW_RIGHT)
   add_function_to_menu('VIEW',VIEW_LEFT)
   add_function_to_menu('VIEW',VIEW_FRONT)
   add_function_to_menu('VIEW',VIEW_REAR)
   add_function_to_menu('VIEW',VIEW_ISO)
   add_function_to_menu('VIEW',VIEW_FITALL)
   start_display()
else:
   mbtb.MBTB_run(stepNumber)
