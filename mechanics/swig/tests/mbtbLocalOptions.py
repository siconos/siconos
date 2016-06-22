"""! \addtogroup mbtbLocalOption

This documentation is build from the python example Tests/SliderCrank, it contains an example of local options.

This file is optional, if it doesn't exit, the default voptions will be use.

@{
"""
import array
import numpy
#bodies view Param
## To define if a body is drew.
bodyDraw=array.array('I',[1,1,1,0])
## Ti define the transparency level of the body.
bodyTrans=array.array('d',[
        2.5,
        0.7,
        0.5])

# contacts view parameters
ContactArtefactLength=0.1
ArtefactThershold=1e-12
## To define if the fisrt surface of the contact is drew.
contactDraw1=array.array('I',[
        1,
        1,
        1])
## To define if the second surface of the contact is drew.
contactDraw2=array.array('I',[
        0,
        0,
        1])



## To define the transparency level of the first contact surface.
contactTrans1=array.array('d',[
        2.,
        2.,
        2.])
## To define the transparency level of the second contact surface.
contactTrans2=array.array('d',[
        0.7,
        0.7,
        2.])

#3D parameters
## It must be set to 1 to run in a 3D view.
with3D = 0
## 3D viewer update frequency and output frequency.
freqOutput=10
freqUpdate=10

apple=0


#Simulation parameters
## Simulation parameters time step size.
stepSize=1e-4
## Simulation parameters number of steps. Useful if with3D=0.
stepNumber=20

TSTheta=0.5
TSGamma=0.5

TSNewtonTolerance=1e-10
TSNewtonMaxIteration=15


TSdeactivateYPosThreshold=1e-5
TSdeactivateYVelThreshold=0.0
TSactivateYPosThreshold=1.e-7
TSactivateYVelThreshold=100


TSProjectionMaxIteration=100
TSConstraintTol=1e-8
TSConstraintTolUnilateral=1e-7
TSLevelOfProjection=0

#solver parameters
## To activate the projection algorithm.
withProj=2
## Solver option
withReduced=2

## Solver option
solverTol=1e-10
## Solver option
solverIt=1000

gotoPos=0

## The number of artefacts
NBARTEFACTS=1
## CAD file of the artefacts
Artefactfile=[
        './CAD/artefact2.step']
## transparency of the artefacts
ArtefactTrans=array.array('d',[
        0.8])
##! @}
