
"""! \addtogroup SliderCrankexample

This documentation is build from the python example Tests/SliderCrank, it contains an example of a simple multibody systems.

This example is the body system corresponding to the example developed in the paper :
MODELING AND ANALYSIS OF MULTIBODY SYSTEMS WITH TRANSLATIONAL CLEARANCE JOINTS BASED ON THE NSDS 
Paulo Flores, Remco Leine, Christoph Glocker.
python ../run.py and result.gp lead to the Fig 11 (d) of this paper.

This example is based on cad files located in Multibody/Tests/CAO/SliderCrank.

The option WITH_CLEARANCE_ON_RODE can be set to 1 to add clearance between the rode 1 and the rode 2.
@{
"""



## A user parameter of this example.
#The option WITH_CLEARANCE_ON_RODE can be set to 1 to add clearance between the rode 1 and the rode 2.
WITH_CLEARANCE_ON_RODE=1
## A user parameter of this example.
l1=0.153
## A user parameter of this example.
l2=0.306
## A user parameter of this example.
a=0.05
## A user parameter of this example.
b=0.025
## A user parameter of this example.
c=0.001
## A user parameter of this example.
w10=-150
## A user parameter of this example.
w20=75
## A user parameter of this example.
w30=0

## REQUIRE BODIES DESCRIPTION
# REQUIRED number of bodies
NBBODIES=3
## Identifier of the word, an object attached to the referential frame.
BATIE=-1
## Identifier of the first body.
PART1=0
## Identifier of the second body.
PART2=1
## Identifier of the third body.
PISTON=2

## REQUIRED the name of the bodies
body=numpy.array(['part1','part2','slider'])
## REQUIRED the initial position of the bodies
initPos=numpy.array([(0,0,0, 0,0,1,0),
                     (l1,0.0,0.0, 0,0,1,0),
                     (l1+l2-a,0,0, 0,1,0,0)])
## REQUIRED the initial velocyties of the bodies
initVel=numpy.array([(0,0,-0.5*w10*l1,0,w10,0),
                     (0,0,-0.5*w10*l1,0,w20,0),
                     (0,0,0,0,w30,0)])

## REQUIRED Initial position of the center of mass.
# It is the position of the center of mass just after loading, ie with the position (0,0,0,1,0,0,0).
initCenterMass=numpy.array([(0.5*l1,0.00,0.00),
                            (0.5*l2,0.00,0.00),
                            (a,0,0)])
## REQUIRED mass of the bodies.
m=array.array('d',[0.038,0.038,0.076])
## REQUIRED inertial matrices.
inertialMatrix=numpy.array([((1,0,0),(0,7.4e-5,0),(0,0,1)),
                            ((1,0,0),(0,5.9e-4,0),(0,0,1)),
                            ((1,0,0),(0,2.7e-6,0),(0,0,1))])
inertialMatrix=numpy.array([((1,0,0),(0,7.4e-5,0),(0,0,1)),
                            ((1,0,0),(0,5.9e-4,0),(0,0,1)),
                            ((1,0,0),(0,2.7e-6,0),(0,0,1))])
## REQUIRED the CAD files.
afile=numpy.array(['./CAD/body1.step',
                   './CAD/body2.step',
                   './CAD/Slider.step'])

## REQUIRED the library for the plugged forces.
#if apple :
#    plugin="SliderCrankPlugin.dylib"
#else :
plugin="SliderCrankPlugin.so"

## REQUIRED the external forces.
fctfext=numpy.array(['externalForcesB1','externalForcesB2','externalForcesS'])
## REQUIRED the external momentums.
#fctmext=numpy.array(['','',''])

## REQUIRED the internal forces.
#fctfint=numpy.array(['internalForcesB1','',''])
## REQUIRED the internal momentums.
#fctmint=numpy.array(['internalMomentsB1','',''])

## REQUIRED the internal forces.
#fctfintjacq=numpy.array(['internalForcesB1_Jacq','',''])
## REQUIRED the internal momentums.
#fctmintjacq=numpy.array(['internalMomentsB1_Jacq','',''])

## REQUIRED the internal forces.
#fctfintjacv=numpy.array(['','',''])
## REQUIRED the internal momentums.
#fctmintjacv=numpy.array(['','',''])

## REQUIRED Boundary condition
# we impose a boundary condition for the foirst body given by the plugin
boundaryCondition=numpy.array(['prescribedvelocityB1','',''])
# we prescribe the 4th component of the velocity of the first body
# i.e we prescribe a angular velocity around the y-axis.
boundaryConditionIndex=numpy.array([numpy.array([4]), numpy.array([]),numpy.array([])]) # we impose the velocity


#JOINTS DESCRIPTION
## REQUIRED the number of joints
NBJOINTS=3
if WITH_CLEARANCE_ON_RODE:
    NBJOINTS=2
else:
    NBJOINTS=3


## REQUIRED the name of the joints
jointName=numpy.array(['Part1_0',
                       'Part2_Piston',
                       'Part1_2'])
## REQUIRED the joint type.
jointType=array.array('I',[mbtb.PIVOT_0,
                           mbtb.PIVOT_1,
                           mbtb.PIVOT_1])
## REQUIRED the index of the first body involved in the joint.
jointBody1=array.array('I',[PART1,
                            PART2,
                            PART1])
## REQUIRED the index of the second body involved in the joint.
jointBody2=array.array('I',[0,
                            PISTON,
                            PART2])
## REQUIRED the joint position.
jointPos=numpy.array([(0,1,0,  0.00,00,00),
                      (0,1,0,  l2,00,00),
                      (0,1,0,  l1,00,00)])


#CONTACTS DESCRIPTION
## REQUIRED the number of contacts
NBCONTACTS=3
if WITH_CLEARANCE_ON_RODE:
    NBCONTACTS=3
else:
    NBCONTACTS=2

## REQUIRED the names of each contact.
contactName=numpy.array([
        'contact_h',
        'contact_b',
        'contact_boby1_2'])
## REQUIRED the CAD file attached to the first body involved in the contact.
afileContact1=numpy.array([
       './CAD/contact_b_cyl.step',
       './CAD/contact_h_cyl.step',
       './CAD/RingBody1.stp'
        ])
## REQUIRED the CAD file attached to the second body involved in the contact.
afileContact2=numpy.array([
       './CAD/chamber.step',
       './CAD/chamber.step',
       './CAD/AxisBody2.stp'
        ])
## REQUIRED the identifier of the first body involved in the contact.
contactBody1=array.array('I',[
        PISTON,
        PISTON,
        PART1])
## REQUIRED the identifier of the second body involved in the contact.
contactBody2=array.array('i',[
        BATIE,
        BATIE,
        PART2])
## the artificial offset sub to the geometrical computation.
contactOffset=array.array('d',[
        0.024,
        0.024,
        0.006])
## defining if the offset is applied to the first surface. Useful to place the contact point.
contactOffsetP1=array.array('I',[
        0,
        0,
        0])
## defining if the normal is computed from the first surface.
contactNormalFromFace1=array.array('I',[
        0,
        0,
        0])
## for each contact, 1 if 3D friction, 0 if perfect unilateral constraints.
contactType3D=array.array('I',[
        1,
        1,
        1])
## for each contact, the mu parameter, useful in the case of friction.
contactmu=array.array('d',[
        0.01,
        0.01,
        0.01])
## the restitution coeficient of each contact.
contacten=array.array('d',[
        0.4,
        0.4,
        0.0])
##! @}
