# Example 3 - Shaft and tube assembly for contact detection.
import math
# REQUIRED number of bodies
NBBODIES=1
## Identifier of the word, an object attached to the referential frame.
BATIE=-1
## Identifier of the first body.
PART1=0

## REQUIRED the name of the bodies
body=numpy.array(['part1'])
## REQUIRED the initial position of the bodies
initPos=numpy.array([(0,0,0, 1,0,0*math.pi/180.0,0*math.pi/180.0)]) #0*math.pi/180.0
## REQUIRED the initial velocyties of the bodies
initVel=numpy.array([(0,0,0,0,0,0)])

## REQUIRED Initial position of the center of mass.
# It is the position of the center of mass just after loading, ie with the position (0,0,0,1,0,0,0).
initCenterMass=numpy.array([(0.0,0.00,0)])
## REQUIRED mass of the bodies.
m=array.array('d',[7.6e-2])
## REQUIRED inertial matrices.
inertialMatrix=numpy.array([((1.91e+1,0,0),(0,1.91e+1,0),(0,0,2.68e+1))])

## REQUIRED the CAD files.
afile=['./CAD/tube.stp']

## REQUIRED the library for the pluged forces.
if apple :
    plugin='TubePlugin.dylib'
else :
    plugin='TubePlugin'    

## REQUIRED the external forces.
fctf=numpy.array(['externalForceG'])
## REQUIRED the external momentums.
fctm=numpy.array(['externalMomentumY'])
#Description###################################################################################################
NBJOINTS=0 # number of joints


#CONTACTS DESCRIPTION
## REQUIRED the number of contacts
NBCONTACTS=3
## REQUIRED the names of each contact.
contactName=[
         'contact_boby0_1','contact_boby1_2','contact_boby2_3']
## REQUIRED the CAD file attached to the first body involved in the contact.
afileContact1=[
       './CAD/1.stp','./CAD/tube_bottom_surface.stp','./CAD/tube_top_side_surface.stp'
       ]
## REQUIRED the CAD file attached to the second body involved in the contact.
afileContact2=[
       './CAD/shaft_surface.stp','./CAD/shaft_surface.stp','./CAD/shaft_top_surface.stp'
        ]
## REQUIRED the identifier of the first body involved in the contact.
contactBody1=array.array('I',[
             PART1,PART1,PART1])
## REQUIRED the identifier of the second body involved in the contact.
contactBody2=array.array('i',[
        BATIE,BATIE,BATIE
        ])
## the artificial offset sub to the geometrical computation.
contactOffset=array.array('d',[
        0.2,0.2,0.2
       ])
## defining if the offset is applied to the first surface. Useful to place the contact point.
contactOffsetP1=array.array('I',[
        0,0,0
        ])
## defining if the normal is computed from the first surface.
contactNormalFromFace1=array.array('I',[
        0,0,0
        ])
## for each contact, 1 if 3D friction, 0 if perfect unilateral constraints.
contactType3D=array.array('I',[
        1,1,1
       ])
## for each contact, the mu parameter, useful in the case of friction.
contactmu=array.array('d',[
        0.001, 0.001,0.001
        ])
## the restitution coeficient of each contact.
contacten=array.array('d',[
        0.0,0.0,0.0
        ])
##! @}
