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

/*! \file SpaceFilter.hpp
 *  \brief Spatial filtering of interactions for 2D/3D objects
 */

/** Basic broad phase contact detection between 2D/3D mechanical systems
 *
 *  algorithm description:
 *   Optimized Spatial Hashing for Collision Detection of Deformable Objects
 *   M. Teschner, B. Heidelberger, M. Mueller, D. Pomeranets, M. Gross
 *   Proceedings of VMV'03
 *   Munich, Germany
 *   pp. 47-54
 *   November 19-21, 2003
 */

#ifndef SpaceFilter_hpp
#define SpaceFilter_hpp

#include <memory>
#include <unordered_set>

#include "BodiesAndRelationsDeclaration.hpp"  // For CircularDS, SphereNEDS, DiskPlanR ...
#include "DynamicalSystemVisitor.hpp"
#include "FMatrix.hpp"  // For FMatrix
#include "InteractionManager.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"

namespace siconos::collision::native {

namespace internal {
class Hashed;
}  // namespace internal

using DiskPlanRDeclared = std::array<double, 6>;

class SpaceFilter : public siconos::simulation::InteractionManager {
 protected:
  ACCEPT_SERIALIZATION(SpaceFilter);

  using SpaceHash = std::unordered_multiset<std::shared_ptr<internal::Hashed>,
                                            boost::hash<std::shared_ptr<internal::Hashed>>>;

  /* relations pool */
  using CircleCircleRDeclared = std::pair<double, double>;
  using DiskDiskRDeclared = std::pair<double, double>;
  using CircleCircleRDeclaredPool =
      std::map<CircleCircleRDeclared,
               std::shared_ptr<siconos::collision::native::bodies::CircularR>>;
  using DiskDiskRDeclaredPool =
      std::map<DiskDiskRDeclared,
               std::shared_ptr<siconos::collision::native::bodies::CircularR>>;
  using DiskPlanRDeclaredPool =
      std::map<DiskPlanRDeclared,
               std::shared_ptr<siconos::collision::native::bodies::DiskPlanR>>;

  /** the bounding box factor is multiplicated by the largest object
      dimension */
  unsigned int _bboxfactor{0};

  /** the cell size */
  unsigned int _cellsize{0};

  /** plans */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _plans{nullptr};

  /** moving plans */
  std::shared_ptr<siconos::collision::native::FMatrix> _moving_plans{nullptr};

  /* the hash table */
  std::shared_ptr<SpaceHash> _hash_table{nullptr};

  /* relations pool */
  std::shared_ptr<DiskDiskRDeclaredPool> diskdisk_relations{nullptr};
  std::shared_ptr<DiskPlanRDeclaredPool> diskplan_relations{nullptr};
  std::shared_ptr<CircleCircleRDeclaredPool> circlecircle_relations{nullptr};

  void _PlanCircularFilter(std::shared_ptr<siconos::simulation::Simulation>, double A,
                           double B, double C, double xCenter, double yCenter, double width,
                           std::shared_ptr<siconos::collision::native::bodies::CircularDS> ds);

  void _MovingPlanCircularFilter(
      std::shared_ptr<siconos::simulation::Simulation>, unsigned int i,
      std::shared_ptr<siconos::collision::native::bodies::CircularDS> ds, double time);

  void _PlanSphereLDSFilter(std::shared_ptr<siconos::simulation::Simulation>, double A,
                            double B, double C, double D,
                            std::shared_ptr<siconos::collision::native::bodies::SphereLDS> ds);

  void _PlanSphereNEDSFilter(
      std::shared_ptr<siconos::simulation::Simulation>, double A, double B, double C, double D,
      std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> ds);

  // Visitors for each kind of proximity detection
  struct CircularFilterVisitor_;
  struct SphereLDSFilterVisitor_;
  struct SphereNEDSFilterVisitor_;

  /* The body hasher, a visitor

     Its visit functions are supposed to insert a specific ds into
     the SpaceFilter hash_table.
   */
  struct BodyHashVisitor_;

  /* the proximity detection */
  struct FindInteractionsVisitor_;

  /* to compare relation */
  struct _IsSameDiskPlanR;
  struct _IsSameDiskMovingPlanR;
  struct _IsSameSpherePlanR;

  /* to compute distance */
  struct DiskDistanceVisitor_;

  // friend struct SpaceFilter::_CircularFilter;
  // friend struct SpaceFilter::_SphereLDSFilter;
  // friend struct SpaceFilter::_SphereNEDSFilter;
  // friend struct SpaceFilter::_BodyHash;
  // friend struct SpaceFilter::_FindInteractions;
  // friend struct SpaceFilter::_IsSameDiskPlanR;
  // friend struct SpaceFilter::_IsSameDiskMovingPlanR;
  // friend struct SpaceFilter::_IsSameSpherePlanR;
  // friend struct SpaceFilter::_DiskDistance;

  SpaceFilter() = delete;  // Do we need a default?
  SpaceFilter(const SpaceFilter&) = delete;
  SpaceFilter(SpaceFilter&&) = delete;
  SpaceFilter operator=(const SpaceFilter&) = delete;
  SpaceFilter operator=(SpaceFilter&&) = delete;

 public:
  SpaceFilter(unsigned int bboxfactor, unsigned int cellsize,
              std::shared_ptr<siconos::algebra::SiconosMatrix> plans,
              std::shared_ptr<siconos::collision::native::FMatrix> moving_plans);

  SpaceFilter(unsigned int bboxfactor, unsigned int cellsize,
              std::shared_ptr<siconos::algebra::SiconosMatrix> plans);

  // SpaceFilter::SpaceFilter()
  //     : SpaceFilter(0, 0, nullptr, nullptr)
  // {
  // }

  /** Destructor */
  virtual ~SpaceFilter() noexcept = default;

  /** 2D/3D objects insertion
   *
   */
  void insert(std::shared_ptr<siconos::collision::native::bodies::Disk>, int, int, int);

  void insert(std::shared_ptr<siconos::collision::native::bodies::Circle>, int, int, int);

  void insert(std::shared_ptr<siconos::collision::native::bodies::SphereLDS>, int, int, int);

  void insert(std::shared_ptr<siconos::collision::native::bodies::SphereNEDS>, int, int, int);

  /** general hashed object
   */
  void insert(std::shared_ptr<siconos::collision::native::internal::Hashed>);

  /** get parameters
   */

  inline unsigned int bboxfactor() { return _bboxfactor; };
  inline unsigned int cellsize() { return _cellsize; };

  void setBBoxfactor(unsigned int value) { _bboxfactor = value; }

  void setCellsize(unsigned int value) { _cellsize = value; }

  /**
     Just test the presence of neighbours.

     \param h hashed component of a body.
   */
  bool haveNeighbours(std::shared_ptr<siconos::collision::native::internal::Hashed> h);

  /**
     Give the minimal distance.

     \param h hashed component of a body.
   */
  double minDistance(std::shared_ptr<siconos::collision::native::internal::Hashed> h);

  /** Broadphase contact detection: add interactions in indexSet 0.
   *
   *  \param simulation the current simulation setup
   */
  void updateInteractions(
      std::shared_ptr<siconos::simulation::Simulation> simulation) override;

  void insertLine(double a, double b, double c);
};

bool operator==(DiskPlanRDeclared const& a, DiskPlanRDeclared const& b);

}  // namespace siconos::collision::native
#endif /* SpaceFilter_hpp */
