/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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

#include "SiconosConfig.h"
// #include <iostream>

namespace siconos::modeling {
class DynamicalSystem;
class LagrangianDS;
class NewtonEulerDS;
class LagrangianR;
class Lagrangian2dVelocityR;
class Lagrangian2dVelocityAngularVelocityR;
class NewtonEulerR;
class NewtonEulerNormalVelocityR;
class NewtonEulerVelocityR;
class NewtonEulerVelocityAngularVelocityR;

}  // namespace siconos::modeling

namespace siconos::collision {

class ContactR;
class ContactVelocityAngularVelocityR;
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
class DiskDiskR;
class DiskPlanR;
class CircleCircleR;
}  // namespace native::bodies
}  // namespace siconos::collision

#ifdef SICONOS_HAS_OpenCASCADE
namespace siconos::mechanics::occ {
class OccBody;
}
#endif

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

  virtual void visit(const T& x) {}
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

/* The classes that override accept (e.g. OccBody) call visit(const Class&): each one must be
   listed in the global visitors, otherwise the default visit of the base visitor throws. */
#ifdef SICONOS_HAS_OpenCASCADE
template <typename T>
using OccVisitMaker = VisitMaker<siconos::mechanics::occ::OccBody, T>;
#else
template <typename T>
using OccVisitMaker = T;
#endif

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
                                                          OccVisitMaker<T>>>>>>>>::Action;
};

#ifdef SICONOS_HAS_BULLET

template <typename T>
struct GlobalRelationVisitor {
  using Make = typename VisitMaker<
      siconos::modeling::LagrangianR,
      VisitMaker<
          siconos::modeling::Lagrangian2dVelocityR,
          VisitMaker<
              siconos::modeling::Lagrangian2dVelocityAngularVelocityR,
              VisitMaker<
                  siconos::modeling::NewtonEulerR,
                  VisitMaker<
                      siconos::modeling::NewtonEulerNormalVelocityR,
                      VisitMaker<
                          siconos::modeling::NewtonEulerVelocityR,
                          VisitMaker<
                              siconos::modeling::NewtonEulerVelocityAngularVelocityR,
                              VisitMaker<
                                  siconos::collision::ContactR,
                                  VisitMaker<
                                      siconos::collision::ContactVelocityAngularVelocityR,
                                      VisitMaker<
                                          siconos::collision::Contact2dR,
                                          VisitMaker<
                                              siconos::collision::Contact2d3DR,
                                              VisitMaker<
                                                  siconos::collision::bullet::BulletR,
                                                  VisitMaker<
                                                      siconos::joints::NewtonEulerJointR,
                                                      VisitMaker<
                                                          siconos::joints::CouplerJointR,
                                                          VisitMaker<
                                                              siconos::joints::
                                                                  CylindricalJointR,
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
                                                                                      PivotJointR,
                                                                                  VisitMaker<
                                                                                      siconos::
                                                                                          joints::
                                                                                              PrismaticJointR,

                                                                                      VisitMaker<
                                                                                          siconos::collision::
                                                                                              native::
                                                                                                  bodies::
                                                                                                      DiskDiskR,
                                                                                          VisitMaker<
                                                                                              siconos::collision::
                                                                                                  native::
                                                                                                      bodies::
                                                                                                          DiskPlanR,
                                                                                              VisitMaker<
                                                                                                  siconos::collision::
                                                                                                      native::
                                                                                                          bodies::
                                                                                                              CircleCircleR,
                                                                                                  T>>>>>>>>>>>>>>>>>>>>>>>>::
      Action;
};

#else
template <typename T>
struct GlobalRelationVisitor {
  using Make = typename VisitMaker<
      siconos::modeling::LagrangianR,
      VisitMaker<
          siconos::modeling::Lagrangian2dVelocityR,
          VisitMaker<
              siconos::modeling::Lagrangian2dVelocityAngularVelocityR,
              VisitMaker<
                  siconos::modeling::NewtonEulerR,
                  VisitMaker<
                      siconos::modeling::NewtonEulerNormalVelocityR,
                      VisitMaker<
                          siconos::modeling::NewtonEulerVelocityR,
                          VisitMaker<
                              siconos::modeling::NewtonEulerVelocityAngularVelocityR,
                              VisitMaker<
                                  siconos::collision::ContactR,
                                  VisitMaker<
                                      siconos::collision::ContactVelocityAngularVelocityR,
                                      VisitMaker<
                                          siconos::collision::Contact2dR,
                                          VisitMaker<
                                              siconos::collision::Contact2d3DR,
                                              VisitMaker<
                                                  siconos::joints::PivotJointR,
                                                  VisitMaker<
                                                      siconos::joints::KneeJointR,
                                                      VisitMaker<
                                                          siconos::joints::PrismaticJointR,
                                                          VisitMaker<
                                                              siconos::collision::native::
                                                                  bodies::DiskDiskR,
                                                              VisitMaker<
                                                                  siconos::collision::native::
                                                                      bodies::DiskPlanR,
                                                                  VisitMaker<
                                                                      siconos::collision::
                                                                          native::bodies::
                                                                              CircleCircleR,
                                                                      T>>>>>>>>>>>>>>>>>::
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
