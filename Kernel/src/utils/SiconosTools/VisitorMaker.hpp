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


typedef Visitor < Classes < LagrangianDS, NewtonEulerDS >, 
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


namespace Alternative{

template<typename T, typename Action>
struct Call : public Action
{
  typedef Call<T, Action> type;
  
  using SiconosVisitor::visit;
  
  virtual void visit(const T& x)
  {
    (*this)(x);
  }
};

template<typename T, typename Action>
struct NoCall : public Action
{
  typedef NoCall type;

  using SiconosVisitor::visit;

  virtual void visit(const T& x)
  {
  }



};


template<typename T, typename Pred>
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
                      Call<T, typename Pred::Action>, 
                      NoCall<T, typename Pred::Action> >::type Action;

};


#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#define REGISTER(X) VisitMaker<X, 
#define REGISTER_STRUCT(X)
#define REGISTER_BASE(X, Y) REGISTER(X)
#define REGISTER_BASE_EXTERN(X, Y)

  template<typename T>
  struct GlobalVisitor
  {
    typedef typename 
    KERNEL_CLASSES()
#ifdef HAVE_SICONOS_MECHANICS
    MECHANICS_CLASSES()
#endif
#ifdef HAVE_BULLET
    BULLET_CLASSES()
#endif
    T
 
#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#define REGISTER(X) >
#define REGISTER_STRUCT(X) 
#define REGISTER_BASE(X, Y) REGISTER(X)
#define REGISTER_BASE_EXTERN(X, Y)
    KERNEL_CLASSES()
#ifdef HAVE_SICONOS_MECHANICS
    MECHANICS_CLASSES()
#endif
#ifdef HAVE_BULLET
    BULLET_CLASSES()
#endif
      ::Action Make;
  };


/* hide mpl::vector */
struct empty{};

template<typename T1 = empty, typename T2 = empty, typename T3 = empty, 
         typename T4 = empty>
struct Classes
{
  typedef typename boost::mpl::vector<T1, T2, T3, T4> Base;
};

template<typename T1>
struct Classes<T1, empty, empty, empty>
{
  typedef typename boost::mpl::vector<T1> Base;
};

template<typename T1, typename T2>
struct Classes<T1, T2, empty, empty>
{
  typedef typename boost::mpl::vector<T1, T2> Base;
};

template<typename T1, typename T2, typename T3>
struct Classes<T1, T2, T3, empty>
{
  typedef typename boost::mpl::vector<T1, T2, T3> Base;
};

/* build final visitor */
template<typename C, typename T>
struct Filter
{
  struct _T : public T, public C
  {
    typedef _T Action;
  };

  typedef typename 
  boost::mpl::fold<
    typename C::Base,
    _T,
    VisitMaker<boost::mpl::_2, boost::mpl::_1 >
    >::type Make;
};

template<typename C, typename T>
struct Visitor
{
  typedef typename Filter<C, T>::Make LocalFilter;
  
  typedef typename GlobalVisitor<LocalFilter>::Make Make;

};


}

/* dummy definitions for some types */

/* note one cannot really define is_complete<T>, cf 
 http://stackoverflow.com/questions/8449036/is-it-possible-to-deduce-whether-type-is-incomplete-without-compilation-failure?lq=1
*/


#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#define REGISTER(X)                             \
  struct X                                      \
  {                                             \
    ACCEPT_STD_VISITORS();                      \
  };                                            \

#define REGISTER_STRUCT(X) REGISTER(X)
#define REGISTER_BASE(X, Y) REGISTER(X)
#define REGISTER_BASE_EXTERN(X, Y) REGISTER(X)


#ifndef HAVE_SICONOS_MECHANICS

MECHANICS_CLASSES()

#endif

#ifndef HAVE_BULLET

BULLET_CLASSES()

#endif

#ifndef HAVE_LMGC

LMGC_CLASSES()

#endif

#endif
