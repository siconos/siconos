/* Siconos-Kernel, Copyright INRIA 2005-2010.
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


/** Very simple spatial filtering of interactions between 2D/3D
 *  systems => to see how may be plugged more elaborated collision
 *  detection softwares.
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

#include "SiconosKernel.hpp"
#include "Circle.hpp"
#include "Disk.hpp"
#include "DiskDiskR.hpp"
#include "CircleCircleR.hpp"
#include "DiskPlanR.hpp"
#include "DiskMovingPlanR.hpp"
#include "SphereLDS.hpp"
#include "SphereLDSSphereLDSR.hpp"
#include "SphereLDSPlanR.hpp"
#include "ExternalBody.hpp"

#include <tr1/unordered_set>
#include <boost/throw_exception.hpp>
#include <boost/functional/hash.hpp>

DEFINE_SPTR(Hashed);

typedef std::tr1::unordered_multiset<SP::Hashed, boost::hash<SP::Hashed> > space_hash;

typedef ublas::matrix<FTime, ublas::column_major, std::vector<FTime> > FMatrix;

TYPEDEF_SPTR(FMatrix);

class SpaceFilter : public boost::enable_shared_from_this<SpaceFilter>
{

protected:

  /** the bounding box factor is multiplicated by the largest object
      dimension */
  unsigned int _bboxfactor;

  /** the cell size */
  unsigned int _cellsize;

  /** interaction counter */
  unsigned int _interID;

  /** plans */
  SP::SiconosMatrix _plans;

  /** moving plans */
  SP::FMatrix _moving_plans;

  /** the whole NonSmoothDynamicalSystem */
  SP::NonSmoothDynamicalSystem _nsds;

  /** only one nslaw */
  SP::NonSmoothLaw _nslaw;

  /* the hash table */
  space_hash _hash_table;

  void _PlanCircularFilter(double A, double B, double C,
                           double xCenter, double yCenter, double width,
                           SP::CircularDS ds);

  void _MovingPlanCircularFilter(unsigned int i, SP::CircularDS ds, double time);

  void _PlanSphereLDSFilter(double A, double B, double C, double D,
                            SP::SphereLDS ds);

  /* visitors defined as Inner class */
  /* note : cf Thinking in C++, vol2, the inner class idiom. */

  /* each kind of proximity detection */
  struct _CircularFilter;
  struct _SphereLDSFilter;



  /* the body hasher */
  struct _BodyHash;

  /* the proximity detection */
  struct _FindInteractions;

  /* to compare relation */
  struct _IsSameDiskPlanR;
  struct _IsSameDiskMovingPlanR;
  struct _IsSameSphereLDSPlanR;

  friend class SpaceFilter::_CircularFilter;
  friend class SpaceFilter::_SphereLDSFilter;
  friend class SpaceFilter::_BodyHash;
  friend class SpaceFilter::_FindInteractions;
  friend class SpaceFilter::_IsSameDiskPlanR;
  friend class SpaceFilter::_IsSameDiskMovingPlanR;
  friend class SpaceFilter::_IsSameSphereLDSPlanR;

  SpaceFilter() {};

public:

  SpaceFilter(unsigned int bboxfactor, unsigned int cellsize,
              SP::NonSmoothDynamicalSystem nsds, SP::NonSmoothLaw nslaw, SP::SiconosMatrix plans, SP::FMatrix moving_plans) :
    _bboxfactor(bboxfactor), _cellsize(cellsize), _interID(0),
    _nsds(nsds), _nslaw(nslaw), _plans(plans), _moving_plans(moving_plans)
  {};


  /** 2D/3D objects insertion
   *
   */
  void insert(SP::Disk, int, int, int);

  void insert(SP::Circle, int, int, int);

  void insert(SP::SphereLDS, int, int, int);


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


  /** search potential interactions
   *
   */
  virtual void buildInteractions(double);

};

#endif /* SpaceFilter_hpp */
