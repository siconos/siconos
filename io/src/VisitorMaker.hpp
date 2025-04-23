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

/*! \file VisitorMaker.hpp
  \brief Generation of visitors on base classes
*/

#ifndef VisitorMaker_hpp
#define VisitorMaker_hpp

#include <boost/mpl/fold.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/type_traits.hpp>
#include <iostream>

#include "SiconosConfig.h"
namespace siconos::modeling {
class DynamicalSystem;
class LagrangianDS;
class NewtonEulerDS;
class LagrangianR;
class Lagrangian2d2DR;
class Lagrangian2d3DR;
class NewtonEulerR;
class NewtonEuler1DR;
class NewtonEuler3DR;
class NewtonEuler5DR;

}  // namespace siconos::modeling

namespace siconos::collision {

class ContactR;
class Contact5DR;
class Contact2dR;
class Contact2d3DR;
class RigidBodyDS;
class RigidBody2dDS;
#ifdef SICONOS_HAS_BULLET
namespace bullet {
class BulletR;
}
#endif
namespace native::bodies {
class Disk;
class Circle;
}  // namespace native::bodies
}  // namespace siconos::collision

namespace siconos::joints {
class CouplerJointR;
class CylindricalJointR;
class FixedJointR;
class JointFrictionR;
class JointStopR;
class KneeJointR;
class NewtonEulerJointR;
class PivotJointR;
class PrismaticJointR;
}  // namespace siconos::joints

namespace siconos::internal {

struct TypeNotFound {};

/* final visitor from Action for class T with base type BaseType */
template <typename T, typename Action, typename BaseType>
struct Call : public Action {
  using type = Call<T, Action, BaseType>;

  using Action::visit;

  virtual void visit(const T& x) { (*this)(static_cast<const BaseType&>(x)); }
  // What ??
  // If BaseType is found (if T inherits from BaseType), then visit(T) is transformed
  // into visit(BaseType) ==> we can use a visitor defined on BaseType for T
};

/* if base type is not found then the visitor is an empty function */
template <typename T, typename Action>
struct Call<T, Action, TypeNotFound> : public Action {
  using type = Call<T, Action, TypeNotFound>;

  using Action::visit;

  virtual void visit(const T& x) { std::cout << "Not good vis \n"; }
};

/* search for a base type of T in Pred::Action::Base and creation of a visitor */
template <typename T, typename Pred>
class VisitMaker {
 private:
  using BaseType =
      typename boost::mpl::fold<typename Pred::Action::Base, TypeNotFound,
                                boost::mpl::if_<boost::is_base_of<boost::mpl::_2, T>,
                                                boost::mpl::_2, boost::mpl::_1>>::type;

 public:
  using Action = typename Call<T, typename Pred::Action, BaseType>::type;
  // What ??
  // Look in Pred::Action::Base for the first class which is a base of T.
  // If found, use Call on it --> generate a working visitor. If not generate a visitor which
  // does nothing.
};

/* build the global visitor for all specified classes */
template <typename T>
struct GlobalDSVisitor {
  using Make = typename VisitMaker<
      siconos::modeling::DynamicalSystem,
      VisitMaker<
          siconos::modeling::LagrangianDS,
          VisitMaker<
              siconos::modeling::NewtonEulerDS,
              VisitMaker<siconos::collision::native::bodies::Disk,
                         VisitMaker<siconos::collision::native::bodies::Circle,
                                    VisitMaker<siconos::collision::RigidBodyDS,
                                               VisitMaker<siconos::collision::RigidBody2dDS,
                                                          T>>>>>>>::Action;
};

#ifdef SICONOS_HAS_BULLET

template <typename T>
struct GlobalRelationVisitor {
  using Make = typename VisitMaker<
      siconos::modeling::LagrangianR,
      VisitMaker<
          siconos::modeling::Lagrangian2d2DR,
          VisitMaker<
              siconos::modeling::Lagrangian2d3DR,
              VisitMaker<
                  siconos::modeling::NewtonEulerR,
                  VisitMaker<
                      siconos::modeling::NewtonEuler1DR,
                      VisitMaker<
                          siconos::modeling::NewtonEuler3DR,
                          VisitMaker<
                              siconos::modeling::NewtonEuler5DR,
                              VisitMaker<
                                  siconos::collision::ContactR,
                                  VisitMaker<
                                      siconos::collision::Contact5DR,
                                      VisitMaker<
                                          siconos::collision::Contact2dR,
                                          VisitMaker<
                                              siconos::collision::Contact2d3DR,
                                              VisitMaker<
                                                  siconos::collision::bullet::BulletR,
                                                  VisitMaker<
                                                      siconos::joints::CouplerJointR,
                                                      VisitMaker<
                                                          siconos::joints::CylindricalJointR,
                                                          VisitMaker<
                                                              siconos::joints::FixedJointR,
                                                              VisitMaker<
                                                                  siconos::joints::
                                                                      JointFrictionR,
                                                                  VisitMaker<
                                                                      siconos::joints::
                                                                          JointStopR,
                                                                      VisitMaker<
                                                                          siconos::joints::
                                                                              KneeJointR,
                                                                          VisitMaker<
                                                                              siconos::joints::
                                                                                  NewtonEulerJointR,
                                                                              VisitMaker<
                                                                                  siconos::joints::
                                                                                      PivotJointR,
                                                                                  VisitMaker<
                                                                                      siconos::
                                                                                          joints::
                                                                                              PrismaticJointR,

                                                                                      T>>>>>>>>>>>>>>>>>>>>>::
      Action;
};

#else
template <typename T>
struct GlobalRelationVisitor {
  using Make = typename VisitMaker<
      siconos::modeling::LagrangianR,
      VisitMaker<
          siconos::modeling::Lagrangian2d2DR,
          VisitMaker<
              siconos::modeling::Lagrangian2d3DR,
              VisitMaker<
                  siconos::modeling::NewtonEulerR,
                  VisitMaker<
                      siconos::modeling::NewtonEuler1DR,
                      VisitMaker<
                          siconos::modeling::NewtonEuler3DR,
                          VisitMaker<
                              siconos::modeling::NewtonEuler5DR,
                              VisitMaker<
                                  siconos::collision::ContactR,
                                  VisitMaker<
                                      siconos::collision::Contact5DR,
                                      VisitMaker<
                                          siconos::collision::Contact2dR,
                                          VisitMaker<
                                              siconos::collision::Contact2d3DR,
                                              VisitMaker<
                                                  siconos::joints::PivotJointR,
                                                  VisitMaker<siconos::joints::KneeJointR,
                                                             VisitMaker<siconos::joints::
                                                                            PrismaticJointR,
                                                                        T>>>>>>>>>>>>>>::
      Action;
};

#endif

template <typename... Ts>
struct Classes {
  using Base = boost::mpl::vector<Ts...>;
};

template <typename C, typename T>
struct Filter {
  struct _T : public T, public C {
    using Action = _T;
  };

  using Make = typename boost::mpl::fold<typename C::Base, _T,
                                         VisitMaker<boost::mpl::_2, boost::mpl::_1>>::type;
};

template <typename C, typename T>
struct DSVisitor {
  using LocalFilter = typename Filter<C, T>::Make;

  using Make = typename GlobalDSVisitor<LocalFilter>::Make;
};

template <typename C, typename T>
struct RelationVisitor {
  using LocalFilter = typename Filter<C, T>::Make;

  using Make = typename GlobalRelationVisitor<LocalFilter>::Make;
};

}  // namespace siconos::internal

#endif
