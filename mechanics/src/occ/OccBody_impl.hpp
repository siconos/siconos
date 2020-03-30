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
/** \file OccBody_impl.hpp
    \brief OccBody implementation details
 */

#ifndef OccBody_impl_hpp
#define OccBody_impl_hpp

#include "OccContactShape.hpp"

#include <vector>
#include <boost/tuple/tuple.hpp>

/* same things in BulletDS_impl.hpp */
/* long unsigned int with bullet stuff ? */
typedef std::array<double, 7> OffSet;

struct ContactShapes : public std::vector<boost::tuple<SP::OccContactShape, OffSet, int > >
{
private:
  ACCEPT_SERIALIZATION(ContactShapes);
};

struct TopoDS_Shapes : public std::vector<boost::tuple<SP::TopoDS_Shape, OffSet > >
{
private:
  ACCEPT_SERIALIZATION(TopoDS_Shapes);
};


#endif
