// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file BodiesDeclaration.hpp
 *  \brief Declaration of all bodies classes and all classes used to describe collision between
 * bodies (relations)
 */

// Thif file must only contain forward declarations !  No implementation, no definitions !

#ifndef BODIES_DECL

namespace siconos::collision::native::bodies {

// Bodies
class CircularDS;
class SphereLDS;
class Disk;
class Circle;
class SphereNEDS;
class ExternalBody;

// Relations

class CircularR;
class SphereNEDSSphereNEDSR;
class SphereNEDSPlanR;
class DiskMovingPlanR;
class DiskPlanR;
class SphereLDSPlanR;
class SphereLDSSphereLDSR;
class CircleCircleR;
class DiskDiskR;
}  // namespace siconos::collision::native::bodies

#endif
