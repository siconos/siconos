/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*! \file SpaceFilter_impl.hpp
 *  \brief implementation details for moving plans
 */

#ifndef SpaceFilter_impl_hpp
#define SpaceFilter_impl_hpp

#include <map>

#include <NSLawMatrix.hpp>
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
