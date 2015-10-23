/* Siconos-Kernel, Copyright INRIA 2005-2013.
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

/*! \file SpaceFilter_impl.hpp
 *  \brief implementation details for moving plans
 */

#ifndef SpaceFilter_impl_hpp
#define SpaceFilter_impl_hpp

#include <map>

#include <SpaceFilter.hpp>
#include "DiskMovingPlanR.hpp"
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/unordered_set.hpp>
#include <boost/throw_exception.hpp>
#include <boost/functional/hash.hpp>

class FMatrix  : public ublas::matrix < FTime, ublas::column_major,
                                         std::vector<FTime> >
{
  ACCEPT_SERIALIZATION(FMatrix);
};

class NSLawMatrix : public ublas::symmetric_matrix < SP::NonSmoothLaw >
{
  ACCEPT_SERIALIZATION(NSLawMatrix);
};

class Hashed : public std11::enable_shared_from_this<Hashed>
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(Hashed);

  Hashed() {};

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

class space_hash : public boost::unordered_multiset < SP::Hashed,
                                                       boost::hash<SP::Hashed> >
{
  ACCEPT_SERIALIZATION(space_hash);
};

/* relations pool */
typedef std::pair<double, double> CircleCircleRDeclared;
typedef std::pair<double, double> DiskDiskRDeclared;
typedef std11::array<double, 6> DiskPlanRDeclared;


class CircleCircleRDeclaredPool : public std::map<CircleCircleRDeclared, SP::CircularR>
{
  ACCEPT_SERIALIZATION(CircleCircleRDeclaredPool);
};


class DiskDiskRDeclaredPool : public std::map<DiskDiskRDeclared, SP::CircularR>
{
  ACCEPT_SERIALIZATION(DiskDiskRDeclaredPool);
};


class DiskPlanRDeclaredPool : public std::map<DiskPlanRDeclared, SP::DiskPlanR>
{
  ACCEPT_SERIALIZATION(DiskPlanRDeclaredPool);
};

#endif
