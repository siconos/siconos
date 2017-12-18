/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*! \file SiconosCollisionManager.hpp
\brief A mechanics world is a Siconos InteractionManager that supports
  static contactors and dynamic contactors attached to special
  Dynamical Systems (BodyDS, derived from NewtonEulerDS) found in the
  NonSmoothDynamicalSystem.
*/

#ifndef SiconosCollisionManager_h
#define SiconosCollisionManager_h

#include <InteractionManager.hpp>
#include <SiconosContactor.hpp>
#include <SiconosCollisionQueryResult.hpp>

class SiconosCollisionManager : public InteractionManager
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosCollisionManager);

public:
  SiconosCollisionManager() : InteractionManager() {}
  virtual ~SiconosCollisionManager() {}

  /** An opaque handle can be used to refer to a specific static
   * contactor set previously added to the collision manager. */
  typedef void* StaticContactorSetID;

  /** Remove a body from the collision detector. This must be done
   *  after removing a body from the NonSmoothDynamicalSystem
   *  otherwise contact will occur with a non-graph body which results
   *  in failure. */
  virtual void removeBody(const SP::BodyDS& body) {}

  /** Perform an intersection test on all shapes in the contactors and
   * return a vector of all results, ordered by distance from start.
   \param start The starting point of the line segment in inertial
                frame (world) coordinates.
   \param end The ending point of the line segment in inertial frame
              (world) coordinates.
   \param closestOnly If true, indicates only interested in first
                      result closest to half-space boundary, max size
                      of returned vector = 1.
   \param sorted If true, results are sorted by distance.
   \return A vector of SiconosCollisionQueryResult that contain
           information about the query results.
  */
  virtual std::vector<SP::SiconosCollisionQueryResult>
  lineIntersectionQuery(const SiconosVector& start,
                        const SiconosVector& end,
                        bool closestOnly=false,
                        bool sorted=true)
    { return std::vector<SP::SiconosCollisionQueryResult>(); }

  /** Find all shapes that are within a sphere defined by a point and
   * a radius and return them in an ordered list based on distance to
   * the center.
   \param center The center of the sphere in inertial frame (world) coordinates.
   \param radius The radius of the sphere.
   \param closestOnly If true, indicates only interested in first
                      result closest to half-space boundary, max size
                      of returned vector = 1.
   \param sorted If true, results are sorted by distance.
   \return A vector of SiconosCollisionQueryResult that contain
           information about the query results.
  */
  virtual std::vector<SP::SiconosCollisionQueryResult>
  inSphereQuery(const SiconosVector& center,
                double radius,
                bool closestOnly=false,
                bool sorted=true)
    { return std::vector<SP::SiconosCollisionQueryResult>(); }

  /** Find all shapes that are within a box defined by a center point
   * and a dimensions (3-vector), and return them in an ordered list
   * based on distance to the center.
   \param center The center of the box in inertial frame (world)
   coordinates.
   \param dimensions The dimensions of the box (3-vector).
   \param closestOnly If true, indicates only interested in first
                      result closest to half-space boundary, max size
                      of returned vector = 1.
   \param sorted If true, results are sorted by distance.
   \return A vector of SiconosCollisionQueryResult that contain
           information about the query results.
  */
  virtual std::vector<SP::SiconosCollisionQueryResult>
  inBoxQuery(const SiconosVector& center,
             const SiconosVector& dimensions,
             bool closestOnly=false,
             bool sorted=true)
    { return std::vector<SP::SiconosCollisionQueryResult>(); }

  /** Find all shapes that are inside a half-space, defined by a point
   * and a normal direction.
   \param point The point defining the boundary of the half-space.
   \param normal The normal pointing away from the surface of the half-space.
   \param closestOnly If true, indicates only interested in first
                      result closest to half-space boundary, max size
                      of returned vector = 1.
   \param sorted If true, results are sorted by distance.
   \return A vector of SiconosCollisionQueryResult that contain
           information about the query results.
  */
  virtual std::vector<SP::SiconosCollisionQueryResult>
  inHalfSpaceQuery(const SiconosVector& point,
                   const SiconosVector& normal,
                   bool closestOnly=false,
                   bool sorted=true)
    { return std::vector<SP::SiconosCollisionQueryResult>(); }

public:
  virtual StaticContactorSetID insertStaticContactorSet(
    SP::SiconosContactorSet cs, SP::SiconosVector position = SP::SiconosVector()) = 0;

  virtual bool removeStaticContactorSet(StaticContactorSetID id) = 0;
};

#endif /* SiconosCollisionManager.hpp */
