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
#include <SiconosFwd.hpp>
#include <SiconosSerialization.hpp>
#include <SiconosVisitor.hpp>

/* local forwards (see SpaceFilter_impl.hpp) */
DEFINE_SPTR(space_hash);
DEFINE_SPTR(DiskDiskRDeclaredPool);
DEFINE_SPTR(DiskPlanRDeclaredPool);
DEFINE_SPTR(CircleCircleRDeclaredPool);
DEFINE_SPTR(Hashed);

class SpaceFilter : public std11::enable_shared_from_this<SpaceFilter>
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SpaceFilter);


  /** the bounding box factor is multiplicated by the largest object
      dimension */
  unsigned int _bboxfactor;

  /** the cell size */
  unsigned int _cellsize;

  /** the siconos model */
  SP::Model _model;

  /** nslaws */
  SP::NSLawMatrix _nslaws;

  /** plans */
  SP::SiconosMatrix _plans;

  /** moving plans */
  SP::FMatrix _moving_plans;

  /* kee track of one step ns integrator initialization */
  bool _osnsinit;

  /* the hash table */
  SP::space_hash _hash_table;

  /* relations pool */
  SP::DiskDiskRDeclaredPool  diskdisk_relations;
  SP::DiskPlanRDeclaredPool  diskplan_relations;
  SP::CircleCircleRDeclaredPool circlecircle_relations;

  void _PlanCircularFilter(double A, double B, double C,
                           double xCenter, double yCenter, double width,
                           SP::CircularDS ds);

  void _MovingPlanCircularFilter(unsigned int i,
                                 SP::CircularDS ds,
                                 double time);

  void _PlanSphereLDSFilter(double A, double B, double C, double D,
                            SP::SphereLDS ds);

  void _PlanSphereNEDSFilter(double A, double B, double C, double D,
                             SP::SphereNEDS ds);

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

  SpaceFilter(unsigned int bboxfactor,
              unsigned int cellsize,
              SP::Model model,
              SP::SiconosMatrix plans,
              SP::FMatrix moving_plans);

  SpaceFilter(unsigned int bboxfactor,
              unsigned int cellsize,
              SP::Model model,
              SP::SiconosMatrix plans);

  SpaceFilter(SP::Model model);

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

  /* Insert a NonSmoothLaw between contactors i & contactors j.
   * \param nslaw the non smooth law
   * \param class1 contactor id
   * \param class2 contactor id
   */
  void insert(SP::NonSmoothLaw nslaw, long unsigned int class1,
              long unsigned int class2);

  /** get parameters
   */

  inline unsigned int bboxfactor()
  {
    return _bboxfactor;
  };
  inline unsigned int cellsize()
  {
    return _cellsize;
  };

  /** Get the model.
      \return a Model object.
   */
  SP::Model model()
  {
    return _model;
  };

  /** Get non smooth laws.
      \return a non smooth laws matrix : SP::NSLawMatrix
   */
  SP::NSLawMatrix nslaws()
  {
    return _nslaws;
  };

  /** Get a non smooth law.
   * \param class1 collision group id of first contactor
   * \param class2 collision group id of second contactor
   * \return a pointer on a NonSmoothLaw object
   */
  SP::NonSmoothLaw nslaw(long unsigned int class1, long unsigned class2);

  /** get the neighbours
   * */
//  std::pair<space_hash::iterator, space_hash::iterator> neighbours(SP::Hashed h);


  /** Just test the presence of neighbours.
      \param h hashed component of a body.
   */
  bool haveNeighbours(SP::Hashed h);

  /** Give the minimal distance.
      \param h hashed component of a body.
   */
  double minDistance(SP::Hashed h);

  /** Broadphase contact detection: add interactions in indexSet 0.
   *  \param time the current time.
   */
  virtual void buildInteractions(double time);

  /** Destructor.
   */
  virtual ~SpaceFilter() {};

  /** visitor hook
   */
  ACCEPT_STD_VISITORS();
};

#endif /* SpaceFilter_hpp */
