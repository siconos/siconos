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

/*! \file VisitorMaker.hpp
  \brief Generation of visitors on base classes
*/

#ifndef VisitorMaker_hpp
#define VisitorMaker_hpp

#include <SiconosVisitor.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/fold.hpp>
/* With visitors on base classes, matches of derived classes is possible
   in a templated visitor operator, example:

struct GetVelocity : public SiconosVisitor
{

  SP::SiconosVector result;

  template<typename T>
  void operator()(const T& ds)
  {
    result = ds.velocity();
  }
};


Visitor < Classes < LagrangianDS, NewtonEulerDS >,
          GetVelocity >::Make getVelocity;

SP::SiconosVector q(new SiconosVector(3));
SP::SiconosVector v(new SiconosVector(3));

  (*q)(0) = 0.;
  (*q)(1) = 1.;
  (*q)(2) = 1.;

  (*v)(0) = 0;
  (*v)(1) = 0;
  (*v)(2) = 10.;

SP::DynamicalSystem ds(new Disk(1, 1, q, v));

ds->accept(getVelocity)->display();


*/

namespace Experimental {

  template<typename T, typename Action, typename IsModifier>
struct Call : public Action
{
  typedef Call<T, Action, IsModifier> type;
  typedef typename Action::arguments_type arguments_type;
  typedef typename boost::mpl::if_<IsModifier, T, const T>::type visitable_type;

  Call() : Action() {};
  Call(arguments_type& args) : Action(args) {};

  using Action::visit;

  virtual void visit(visitable_type& x)
  {
    (*this)(x);
  }

};

  template<typename T, typename Action, typename IsModifier>
struct NoCall : public Action
{
  typedef NoCall<T, Action, IsModifier> type;

  typedef typename Action::arguments_type arguments_type;
  typedef typename boost::mpl::if_<IsModifier, T, const T>::type visitable_type;

  NoCall() : Action() {};
  NoCall(arguments_type& args) : Action(args) {};

  using Action::visit;

  virtual void visit(visitable_type& x)
  {
  }

};


  template<typename T, typename Pred, typename IsModifier>
class VisitMaker
{
private:
  typedef typename
  boost::mpl::fold<
  typename Pred::Action::Base,
  boost::mpl::false_,
  boost::mpl::if_<boost::is_base_of<boost::mpl::_2, T>,
                  boost::mpl::true_,
                  boost::mpl::_1> >::type Condition;

public:
  typedef typename
  boost::mpl::eval_if<Condition,
                      Call<T, typename Pred::Action, IsModifier>,
                      NoCall<T, typename Pred::Action, IsModifier> >::type Action;

};


#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#define REGISTER(X) VisitMaker<X,
#define REGISTER_STRUCT(X)
#define REGISTER_BASE(X, Y) REGISTER(X)
#define REGISTER_BASE_EXTERN(X, Y)

  template<typename T, typename IsModifier>
  struct GlobalVisitor
  {
    typedef typename
    VISITOR_CLASSES()
    T

#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#define REGISTER(X) , IsModifier >
#define REGISTER_STRUCT(X)
#define REGISTER_BASE(X, Y) REGISTER(X)
#define REGISTER_BASE_EXTERN(X, Y)
    VISITOR_CLASSES()
      ::Action Make;
  };


/* hide mpl::vector */
struct empty{};

template<typename T1 = empty, typename T2 = empty, typename T3 = empty,
         typename T4 = empty, typename T5 = empty, typename T6 = empty,
         typename T7 = empty, typename T8 = empty, typename T9 = empty>
struct Classes
{
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6, T7, T8, T9> Base;
};

template<typename T1>
struct Classes<T1, empty, empty, empty, empty, empty, empty, empty, empty>
{
  typedef typename boost::mpl::vector<T1> Base;
};

template<typename T1, typename T2>
struct Classes<T1, T2, empty, empty, empty, empty, empty, empty, empty>
{
  typedef typename boost::mpl::vector<T1, T2> Base;
};

template<typename T1, typename T2, typename T3>
struct Classes<T1, T2, T3, empty, empty, empty, empty, empty, empty>
{
  typedef typename boost::mpl::vector<T1, T2, T3> Base;
  };

template<typename T1, typename T2, typename T3, typename T4>
struct Classes<T1, T2, T3, T4, empty, empty, empty, empty, empty>
{
  typedef typename boost::mpl::vector<T1, T2, T3, T4> Base;
};

template<typename T1, typename T2, typename T3, typename T4, typename T5>
struct Classes<T1, T2, T3, T4, T5, empty, empty, empty, empty>
{
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5> Base;
};

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct Classes<T1, T2, T3, T4, T5, T6, empty, empty, empty>
{
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6> Base;
};

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
struct Classes<T1, T2, T3, T4, T5, T6, T7, empty, empty>
{
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6, T7> Base;
};

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
struct Classes<T1, T2, T3, T4, T5, T6, T7, T8, empty>
{
  typedef typename boost::mpl::vector<T1, T2, T3, T4, T5, T6, T7, T8> Base;
};

/* build final visitor */
  template<typename C, typename T, typename IsModifier>
struct Filter
{
  struct _T : public T, public C
  {
    _T() : T() {};
    _T(typename T::arguments_type& args) : T(args) {};
    typedef _T Action;
  };

  typedef typename
  boost::mpl::fold<
    typename C::Base,
    _T,
    VisitMaker<boost::mpl::_2, boost::mpl::_1, IsModifier>
    >::type Make;
  };

  template<typename C, typename T, typename IsModifier=boost::mpl::true_>
struct Visitor
{
  typedef typename Filter<C, T, IsModifier>::Make LocalFilter;

  typedef typename GlobalVisitor<LocalFilter, IsModifier>::Make Make;
};


  typedef boost::mpl::true_ VisitorWriter;
  typedef boost::mpl::false_ VisitorReader;

}

#endif
