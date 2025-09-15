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

#include "MechanicsFwd.hpp"

#include <InteractionManager.hpp>
#include <SiconosFwd.hpp>
#include <SiconosSerialization.hpp>
#include <SiconosVisitor.hpp>

/* local forwards (see SpaceFilter_impl.hpp) */
DEFINE_SPTR(space_hash);
DEFINE_SPTR(DiskDiskRDeclaredPool);
DEFINE_SPTR(DiskPlanRDeclaredPool);
DEFINE_SPTR(CircleCircleRDeclaredPool);
DEFINE_SPTR(Hashed);

class SpaceFilter : public InteractionManager,
                    public std::enable_shared_from_this<SpaceFilter> {

protected:

  ACCEPT_SERIALIZATION(SpaceFilter);

  /** the bounding box factor is multiplicated by the largest object
      dimension */
  unsigned int _bboxfactor;

  /** the cell size */
  unsigned int _cellsize;

  /** plans */
  SP::SiconosMatrix _plans;

  /** moving plans */
  SP::FMatrix _moving_plans;

  /* the hash table */
  SP::space_hash _hash_table;

  /* relations pool */
  SP::DiskDiskRDeclaredPool diskdisk_relations;
  SP::DiskPlanRDeclaredPool diskplan_relations;
  SP::CircleCircleRDeclaredPool circlecircle_relations;

  void _PlanCircularFilter(SP::Simulation, double A, double B, double C,
                           double xCenter, double yCenter, double width,
                           SP::CircularDS ds);

  void _MovingPlanCircularFilter(SP::Simulation, unsigned int i,
                                 SP::CircularDS ds, double time);

  void _PlanSphereLDSFilter(SP::Simulation, double A, double B, double C,
                            double D, SP::SphereLDS ds);

  void _PlanSphereNEDSFilter(SP::Simulation, double A, double B, double C,
                             double D, SP::SphereNEDS ds);

  /* visitors defined as Inner class */
  /* note : cf Thinking in C++, vol2, the inner class idiom. */

  /* each kind of proximity detection */
  struct _CircularFilter;
  struct _SphereLDSFilter;
  struct _SphereNEDSFilter;

  /* the body hasher */
  struct _BodyHash;

  /* the proximity detection */
  struct _FindInteractions;

  /* to compare relation */
  struct _IsSameDiskPlanR;
  struct _IsSameDiskMovingPlanR;
  struct _IsSameSpherePlanR;

  /* to compute distance */
  struct _DiskDistance;

  friend struct SpaceFilter::_CircularFilter;
  friend struct SpaceFilter::_SphereLDSFilter;
  friend struct SpaceFilter::_SphereNEDSFilter;
  friend struct SpaceFilter::_BodyHash;
  friend struct SpaceFilter::_FindInteractions;
  friend struct SpaceFilter::_IsSameDiskPlanR;
  friend struct SpaceFilter::_IsSameDiskMovingPlanR;
  friend struct SpaceFilter::_IsSameSpherePlanR;
  friend struct SpaceFilter::_DiskDistance;

public:
  SpaceFilter(unsigned int bboxfactor, unsigned int cellsize,
              SP::SiconosMatrix plans, SP::FMatrix moving_plans);

  SpaceFilter(unsigned int bboxfactor, unsigned int cellsize,
              SP::SiconosMatrix plans);

  SpaceFilter();

  /** 2D/3D objects insertion
   *
   */
  void insert(SP::Disk, int, int, int);

  void insert(SP::Circle, int, int, int);

  void insert(SP::SphereLDS, int, int, int);

  void insert(SP::SphereNEDS, int, int, int);

  /** general hashed object
   */
  void insert(SP::Hashed);

  /** get parameters
   */

  inline unsigned int bboxfactor() { return _bboxfactor; };
  inline unsigned int cellsize() { return _cellsize; };

  void setBBoxfactor(unsigned int value) { _bboxfactor = value; }

  void setCellsize(unsigned int value) { _cellsize = value; }

  /** get the neighbours
   * */
  //  std::pair<space_hash::iterator, space_hash::iterator>
  //  neighbours(SP::Hashed h);

  /**
     Just test the presence of neighbours.
     
     \param h hashed component of a body.
   */
  bool haveNeighbours(SP::Hashed h);

  /**
     Give the minimal distance.
     
     \param h hashed component of a body.
   */
  double minDistance(SP::Hashed h);

  /** Broadphase contact detection: add interactions in indexSet 0.
   *
   *  \param simulation the current simulation setup
   */
  void updateInteractions(SP::Simulation simulation) override;

  void insertLine(double a, double b, double c);

  /** Destructor.
   */
  virtual ~SpaceFilter(){};

  ACCEPT_STD_VISITORS();
};

#endif /* SpaceFilter_hpp */
