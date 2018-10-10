"""Simulation of a slider-crank system.
Check examples manual for details.

"""

import numpy as np
import array
import siconos.mechanics.mechanisms.mbtb as mbtb

WITH_CLEARANCE_ON_RODE = 1
"""if true, add clearance between rodes 1 and 2."""

# -- Slider-crank system parameters --

l1 = 0.153
"""Crank length"""
l2 = 0.306
"""Connecting rod length"""
a = 0.05
"""half-length of the slider"""
b = 0.025
"""half-height of the slider"""
c = 0.001
"""clearance between slider and guide"""
w10 = -150.
"""initial angular speed for the crank"""
w20 = 75.
"""initial angular speed for the connecting rod"""
w30 = 0.
"""initial angular speed for the slider"""

NBBODIES = 3
"""Number of bodies """

BATIE = -1
"""Identifier of the world, an object attached to the referential frame."""

PART1 = 0
"""First body id"""

PART2 = 1
"""2nd body id """

PISTON = 2
"""Third body id"""

body = np.array(['part1', 'part2', 'slider'])
""" Bodies names"""

initPos = np.array([(0, 0, 0, 0, 0, 1, 0),
                    (l1, 0.0, 0.0, 0, 0, 1, 0),
                    (l1 + l2 - a, 0, 0, 0, 1, 0, 0)])
"""initial positions of the bodies """

initVel = np.array([(0, 0, -0.5 * w10 * l1, 0, w10, 0),
                    (0, 0, -0.5 * w10 * l1, 0, w20, 0),
                    (0, 0, 0, 0, w30, 0)])
"""initial velocities of the bodies """

initCenterMass = np.array([(0.5 * l1, 0.0, 0.0),
                           (0.5 * l2, 0.0, 0.0),
                           (a, 0., 0.)])
"""initial positions of the centers of mass of the bodies
positions, ie with the position (0,0,0,1,0,0,0)"""

m = array.array('d', [0.038, 0.038, 0.076])
"""masses of the bodies"""

inertialMatrix = np.array([((1, 0, 0), (0, 7.4e-5, 0), (0, 0, 1)),
                           ((1, 0, 0), (0, 5.9e-4, 0), (0, 0, 1)),
                           ((1, 0, 0), (0, 2.7e-6, 0), (0, 0, 1))])
"""inertia matrix"""

afile = ['./CAD/body1.step',
         './CAD/body2.step',
         './CAD/Slider.step']
"""CAD files """

plugin = "SliderCrankPlugin.so"
"""plugins library name"""

fctfext = np.array(['externalForcesB1', 'externalForcesB2', 'externalForcesS'])
"""external forces"""

# REQUIRED the external momentums.
#fctmext = np.array(['','',''])

# REQUIRED the internal forces.
#fctfint = np.array(['internalForcesB1','',''])
# REQUIRED the internal momentums.
#fctmint = np.array(['internalMomentsB1','',''])

# REQUIRED the internal forces.
#fctfintjacq = np.array(['internalForcesB1_Jacq','',''])
# REQUIRED the internal momentums.
#fctmintjacq = np.array(['internalMomentsB1_Jacq','',''])

# REQUIRED the internal forces.
#fctfintjacv = np.array(['','',''])
# REQUIRED the internal momentums.
#fctmintjacv = np.array(['','',''])

boundaryCondition = np.array(['prescribedvelocityB1', '', ''])
""""a boundary condition (given by the plugin)
 is enforced for the first body"""

boundaryConditionIndex = np.array([np.array([4]),
                                   np.array([]),
                                   np.array([])])
""""we prescribe the 4th component of the velocity of the first body
i.e we prescribe a angular velocity around the y-axis."""

# --- JOINTS DESCRIPTION ---
NBJOINTS = 3
"""number of joints"""

if WITH_CLEARANCE_ON_RODE:
    NBJOINTS = 2
else:
    NBJOINTS = 3


jointName = np.array(['Part1_0',
                      'Part2_Piston',
                      'Part1_2'])
"""joints' names"""

jointType = array.array('I', [mbtb.PIVOT_0,
                              mbtb.PIVOT_1,
                              mbtb.PIVOT_1])
"""joints' types"""

jointBody1 = array.array('I', [PART1,
                               PART2,
                               PART1])
"""index of the first body involved in the joint."""

jointBody2 = array.array('I', [0,
                               PISTON,
                               PART2])
"""index of the second body involved in the joint."""

jointPos = np.array([(0, 1, 0, 0.0, 0, 0.),
                     (0, 1, 0, l2, 0., 0.),
                     (0, 1, 0, l1, 0., 0)])
"""joints' positions"""

#--- CONTACTS DESCRIPTION ---

NBCONTACTS = 3
"""number of contacts"""
if WITH_CLEARANCE_ON_RODE:
    NBCONTACTS = 3
else:
    NBCONTACTS = 2

contactName = np.array([
    'contact_h',
    'contact_b',
    'contact_boby1_2'])
"""contacts' names"""

afileContact1 = [
    './CAD/contact_b_cyl.step',
    './CAD/contact_h_cyl.step',
    './CAD/RingBody1.stp']
"""CAD files attached to the first body involved in the contact"""

afileContact2 = [
    './CAD/chamber.step',
    './CAD/chamber.step',
    './CAD/AxisBody2.stp']
"""CAD files attached to the second body involved in the contact"""

contactBody1 = array.array('I', [
    PISTON,
    PISTON,
    PART1])
"""identifier of the first body involved in the contact."""

contactBody2 = array.array('i', [
    BATIE,
    BATIE,
    PART2])
"""identifier of the second body involved in the contact"""

contactOffset = array.array('d', [
    0.024,
    0.024,
    0.006])
"""the artificial offset sub to the geometrical computation."""

contactOffsetP1 = array.array('I', [
    0,
    0,
    0])
"""defining if the offset is applied to the first surface.
Useful to place the contact point."""

contactNormalFromFace1 = array.array('I', [
    0,
    0,
    0])
"""defining if the normal is computed from the first surface."""

contactType3D = array.array('I', [
    1,
    1,
    1])
"""for each contact, 1 if 3D friction, 0 if perfect unilateral constraints."""

contactmu = array.array('d', [
    0.01,
    0.01,
    0.01])
""""friction coeff, only useful in the case of friction."""

contacten = array.array('d', [
    0.4,
    0.4,
    0.0])
"""restitution coeff, for each contact"""
