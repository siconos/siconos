/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

#include "SiconosKernel.hpp"
#include "Circle.hpp"
#include "Disk.hpp"
#include "DiskDiskR.hpp"
#include "CircleCircleR.hpp"
#include "DiskPlanR.hpp"
#include "DiskMovingPlanR.hpp"
#include "SphereLDS.hpp"
#include "SphereLDSSphereLDSR.hpp"
#include "SphereNEDSSphereNEDSR.hpp"
#include "SphereLDSPlanR.hpp"
#include "SphereNEDS.hpp"
#include "SphereNEDSPlanR.hpp"
#include "ExternalBody.hpp"

#ifndef __GCCXML__
#include <tr1/unordered_set>
#else
/* gccxml fail to parse some gcc 4.3 builtins, so we provide a fake
 * unordered_multiset*/
namespace std
{
namespace tr1
{
template<class A, class B>
class unordered_multiset
{
public:
  typedef void iterator;
};
}
}
#endif

#include <boost/throw_exception.hpp>
#include <boost/functional/hash.hpp>

/** hash container
 */
class Hashed : public boost::enable_shared_from_this<Hashed>
{
public:
  SP::DynamicalSystem body;
  int i;
  int j;
  int k;
  Hashed(SP::DynamicalSystem body, int i, int j, int k = 0) :
    body(body), i(i), j(j), k(k) {};

  Hashed(int i, int j, int k = 0)
    : i(i), j(j), k(k) {};

  ~Hashed() {};

};

DEFINE_SPTR(Hashed);

typedef std::tr1::unordered_multiset < SP::Hashed,
        boost::hash<SP::Hashed> > space_hash;

typedef ublas::matrix < FTime, ublas::column_major,
        std::vector<FTime> > FMatrix;

TYPEDEF_SPTR(FMatrix);

class SpaceFilter : public boost::enable_shared_from_this<SpaceFilter>
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

  /** the whole NonSmoothDynamicalSystem */
  SP::NonSmoothDynamicalSystem _nsds;

  /** only one nslaw */
  SP::NonSmoothLaw _nslaw;

  /** plans */
  SP::SiconosMatrix _plans;

  /** moving plans */
  SP::FMatrix _moving_plans;

  /* the hash table */
  space_hash _hash_table;

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


  friend class SpaceFilter::_CircularFilter;
  friend class SpaceFilter::_SphereLDSFilter;
  friend class SpaceFilter::_SphereNEDSFilter;
  friend class SpaceFilter::_BodyHash;
  friend class SpaceFilter::_FindInteractions;
  friend class SpaceFilter::_IsSameDiskPlanR;
  friend class SpaceFilter::_IsSameDiskMovingPlanR;
  friend class SpaceFilter::_IsSameSpherePlanR;
  friend class SpaceFilter::_DiskDistance;

public:

  SpaceFilter(unsigned int bboxfactor,
              unsigned int cellsize,
              SP::NonSmoothDynamicalSystem nsds,
              SP::NonSmoothLaw nslaw,
              SP::SiconosMatrix plans,
              SP::FMatrix moving_plans) :
    _bboxfactor(bboxfactor), _cellsize(cellsize), _interID(0),
    _nsds(nsds), _nslaw(nslaw), _plans(plans), _moving_plans(moving_plans)
  {};

  SpaceFilter(unsigned int bboxfactor,
              unsigned int cellsize,
              SP::NonSmoothDynamicalSystem nsds,
              SP::NonSmoothLaw nslaw,
              SP::SiconosMatrix plans) :
    _bboxfactor(bboxfactor), _cellsize(cellsize), _interID(0),
    _nsds(nsds), _nslaw(nslaw), _plans(plans)
  {};

  SpaceFilter()
  {};

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

  inline unsigned int bboxfactor()
  {
    return _bboxfactor;
  };
  inline unsigned int cellsize()
  {
    return _cellsize;
  };

  /** get non smooth dynamical system
   */
  SP::NonSmoothDynamicalSystem nsds()
  {
    return _nsds;
  };


  /** get non smooth law
   */
  SP::NonSmoothLaw nslaw()
  {
    return _nslaw;
  };


  /** get an interaction id
   * */
  unsigned int newInteractionId()
  {
    return _interID++;
  };

  /** get the neighbours
   * */
  std::pair<space_hash::iterator, space_hash::iterator> neighbours(SP::Hashed h);


  /** just test the presence of neighbours
   */
  bool haveNeighbours(SP::Hashed h);

  /** give the minimal distance
   */
  double minDistance(SP::Hashed h);



  /** search potential interactions
   *
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
