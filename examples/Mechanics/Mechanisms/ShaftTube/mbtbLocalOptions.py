"""! \addtogroup mbtbLocalOption

This file is optional, if it doesn't exit, the default voptions will be use.

@{
"""
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
ContactArtefactLength=20.0
ArtefactThershold=1e-12
## To define if the fisrt surface of the contact is drew.
contactDraw1=array.array('I',[
        1,
        1,
        1])
## To define if the second surface of the contact is drew.
contactDraw2=array.array('I',[
        1,
        1,
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
with3D=1
## 3D viewer update frequency and output frequency.
freqOutput=20
freqUpdate=20

apple=0


#Simulation parameters
## Simulation parameters time step size.
stepSize=1e-4
## Simulation parameters number of steps. Useful if with3D=0.
stepNumber=2000

TSdeactivateYPosThreshold=1e-1
TSdeactivateYVelThreshold=0.0
TSactivateYPosThreshold=1.e-1
TSactivateYVelThreshold=100


TSProjectionMaxIteration=100
TSConstraintTol=1e-1
TSConstraintTolUnilateral=1e-1
TSLevelOfProjection=0

#solver parameters
## To activate the projection algorithm.
withProj=1
## Solver option
withReduced=2

## Solver option
solverTol=1e-1
## Solver option
solverIt=100000

gotoPos=0

## The number of artefacts
NBARTEFACTS=1
## CAD file of the artefacts
Artefactfile=[
        './CAD/shaft.stp']
## transparency of the artefacts
ArtefactTrans=array.array('d',[
        0.8])
##! @}
