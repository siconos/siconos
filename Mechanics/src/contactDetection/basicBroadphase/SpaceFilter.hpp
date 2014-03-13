/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

  /** interaction counter */
  unsigned int _interID;

  /** the siconos model */
  SP::Model _model;

  /** nslaws */
  SP::NSLawMatrix _nslaws;

  /** plans */
  SP::SiconosMatrix _plans;

  /** moving plans */
  SP::FMatrix _moving_plans;

  /* the hash table */
  SP::space_hash _hash_table;

  /* kee track of one step ns integrator initialization */
  bool _osnsinit;

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

  /* insert a NonSmoothLaw between contactors i & contactors j
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

  /** get the model
   */
  SP::Model model()
  {
    return _model;
  };


  /** get non smooth laws
      \return SP::NSLawMatrix
   */
  SP::NSLawMatrix nslaws()
  {
    return _nslaws;
  };


  /** get an interaction id
      \return unsigned int
   * */
  unsigned int newInteractionId()
  {
    return _interID++;
  };

  /** get the neighbours
   * */
//  std::pair<space_hash::iterator, space_hash::iterator> neighbours(SP::Hashed h);


  /** just test the presence of neighbours
   */
  bool haveNeighbours(SP::Hashed h);

  /** give the minimal distance
   */
  double minDistance(SP::Hashed h);


  /** insert a new interaction and link it to 1 or 2 ds.
      \param inter the new interaction
      \param ds1 a SP::DynamicalSystem
      \param ds2 a SP::DynamicalSystem (optional)
  */
  void link(SP::Interaction inter, SP::DynamicalSystem,
            SP::DynamicalSystem = SP::DynamicalSystem());

  /** broadphase contact detection: add interactions in indexSet 0
   *  \param the current time
   */
  virtual void buildInteractions(double);


  /** destructor
   */
  virtual ~SpaceFilter() {};


  /** visitor hook
   */
  ACCEPT_STD_VISITORS();
};

#endif /* SpaceFilter_hpp */
