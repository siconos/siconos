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

#include <SiconosVisitor.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/type_traits.hpp>

#include "KneeJointR.hpp"
#include "Lagrangian2d2DR.hpp"
#include "Lagrangian2d3DR.hpp"
#include "NewtonEuler3DR.hpp"
#include "RigidBodyDS.hpp"

/* With visitors on base classes, matches of derived classes is possible
   in a templated visitor operator, example:

struct GetVelocity : public SiconosVisitor
{

std::shared_ptr<SiconosVector> result;

  template<typename T>
  void operator()(const T& ds)
  {
    result = ds.velocity();
  }
};


Visitor < Classes < LagrangianDS, NewtonEulerDS >,
          GetVelocity >::Make getVelocity;

auto q = std::make_shared<SiconosVector>(3);
auto v = std::make_shared<SiconosVector>(3);

  (*q)(0) = 0.;
  (*q)(1) = 1.;
  (*q)(2) = 1.;

  (*v)(0) = 0;
  (*v)(1) = 0;
  (*v)(2) = 10.;

auto ds = std::make_shared<Disk>(1, 1, q, v);

ds->accept(getVelocity)->display();


*/
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

namespace native::bodies {
class Disk;
class Circle;
}  // namespace native::bodies
}  // namespace siconos::collision

namespace siconos::joints {
class PivotJointR;
class KneeJointR;
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
};

/* if base type is not found then the visitor is an empty
   function */
template <typename T, typename Action>
struct Call<T, Action, TypeNotFound> : public Action {
  using type = Call<T, Action, TypeNotFound>;

  using Action::visit;

  virtual void visit(const T& x) {}
};

/* search for a base type of T in Pred::Action::Base and creation
   of a visitor */
template <typename T, typename Pred>
class VisitMaker {
 private:
  using BaseType =
      typename boost::mpl::fold<typename Pred::Action::Base, TypeNotFound,
                                boost::mpl::if_<boost::is_base_of<boost::mpl::_2, T>,
                                                boost::mpl::_2, boost::mpl::_1>>::type;

 public:
  using Action = typename Call<T, typename Pred::Action, BaseType>::type;
};

/* build the global visitor for all specified classes */

template <typename T>
struct GlobalVisitor {
  using Make = VisitMaker<
      siconos::modeling::DynamicalSystem,
      VisitMaker<
          siconos::modeling::LagrangianDS,
          VisitMaker<
              siconos::modeling::NewtonEulerDS,
              VisitMaker<
                  siconos::modeling::LagrangianR,
                  VisitMaker<
                      siconos::collision::native::bodies::Disk,
                      VisitMaker<
                          siconos::collision::native::bodies::Circle,
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
                                                                  siconos::collision::
                                                                      Contact2d3DR,
                                                                  VisitMaker<
                                                                      siconos::joints::
                                                                          PivotJointR,
                                                                      VisitMaker<
                                                                          siconos::joints::
                                                                              KneeJointR,
                                                                          VisitMaker<
                                                                              siconos::joints::
                                                                                  PrismaticJointR,
                                                                              VisitMaker<
                                                                                  siconos::
                                                                                      collision::
                                                                                          RigidBodyDS,
                                                                                  VisitMaker<
                                                                                      siconos::
                                                                                          collision::
                                                                                              RigidBody2dDS,
                                                                                      T>>>>>>>>>>>>>>>>>>>>>::
      Action;
};

/* hide mpl::vector */
struct empty {};

template <typename T1 = empty, typename T2 = empty, typename T3 = empty, typename T4 = empty,
          typename T5 = empty, typename T6 = empty, typename T7 = empty, typename T8 = empty,
          typename T9 = empty>
struct Classes {
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6, T7, T8, T9> Base;
};

template <typename T1>
struct Classes<T1, empty, empty, empty, empty, empty, empty, empty, empty> {
  typedef typename boost::mpl::vector<T1> Base;
};

template <typename T1, typename T2>
struct Classes<T1, T2, empty, empty, empty, empty, empty, empty, empty> {
  typedef typename boost::mpl::vector<T1, T2> Base;
};

template <typename T1, typename T2, typename T3>
struct Classes<T1, T2, T3, empty, empty, empty, empty, empty, empty> {
  typedef typename boost::mpl::vector<T1, T2, T3> Base;
};

template <typename T1, typename T2, typename T3, typename T4>
struct Classes<T1, T2, T3, T4, empty, empty, empty, empty, empty> {
  typedef typename boost::mpl::vector<T1, T2, T3, T4> Base;
};

template <typename T1, typename T2, typename T3, typename T4, typename T5>
struct Classes<T1, T2, T3, T4, T5, empty, empty, empty, empty> {
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5> Base;
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct Classes<T1, T2, T3, T4, T5, T6, empty, empty, empty> {
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6> Base;
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
          typename T7>
struct Classes<T1, T2, T3, T4, T5, T6, T7, empty, empty> {
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6, T7> Base;
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
          typename T7, typename T8>
struct Classes<T1, T2, T3, T4, T5, T6, T7, T8, empty> {
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6, T7, T8> Base;
};

/* build final visitor */
template <typename C, typename T>
struct Filter {
  struct _T : public T, public C {
    typedef _T Action;
  };

  typedef typename boost::mpl::fold<typename C::Base, _T,
                                    VisitMaker<boost::mpl::_2, boost::mpl::_1>>::type Make;
};

template <typename C, typename T>
struct Visitor {
  typedef typename Filter<C, T>::Make LocalFilter;

  typedef typename GlobalVisitor<LocalFilter>::Make Make;
};

}  // namespace siconos::internal

#endif
