/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
