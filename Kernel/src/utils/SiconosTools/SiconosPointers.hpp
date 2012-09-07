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

#ifndef SiconosPointers_hpp
#define SiconosPointers_hpp

/*! \file SiconosPointers.hpp
  Siconos interface to reference-counting pointers
*/

/** Siconos pointers

 \author SICONOS Development Team - copyright INRIA
 \version 3.0.0.
 \date (Creation) 2010

 Siconos pointers are reference counting pointers. Memory pointed by a
 Siconos pointer is automaticaly deallocated.

 Basic usage :

 - to declare a Siconos pointer on a Siconos Object:

 \verbatim
 SP::<Siconos Object> myobject;
 \endverbatim

 - to set a Siconos pointer :

 \verbatim
 myobject.reset(new <Siconos Object>(...));
 \endverbatim


Siconos pointers are boost smart pointers.

More documentation on smart pointers and reference counting:

 - http://en.wikipedia.org/wiki/Smart_pointer

 - http://en.wikipedia.org/wiki/Reference_counting

 - http://www.boost.org/doc/libs/release/libs/smart_ptr

 */

#include <SiconosNumerics.h>

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace SharedPointer {}
namespace SharedPointerConst {}
namespace SharedArray {}

/** \namespace SP Namespace for Siconos smart pointers : memory
    pointed by a Siconos smart pointers is automaticaly
    deallocated */
namespace SP = SharedPointer;

/** \namespace SA Namespace for Siconos shared arrays : Siconos shared
    arrays are automaticaly deallocated */
namespace SA = SharedArray;

/** \namespace SPC Namespace for const shared pointers : memory
    pointed by a const shared pointers is automaticaly deallocated */
namespace SPC = SharedPointerConst;



/** Using a shared_ptr to hold a pointer to a statically allocated
   object
   use create<type>SPtr(<type> &x)
   cf http://www.boost.org/doc/
*/
struct nullDeleter
{
  void operator()(void const *) const {}
};

/* template namespace : no */

#define NAME_SPACE_SPTR(X) \
  namespace SharedPointer \
  { \
    typedef SPtr##X X; \
  } \
  namespace SharedPointerConst \
  { \
    typedef SPtrConst##X X;\
  }

#define NAME_SPACE_SAPTR(X)                     \
  namespace SharedArray \
  { \
    typedef X##SAPtr X; \
  }


/* template typedef : no */

#define TYPEDEF_SPTR(X) \
  typedef boost::shared_ptr<X> SPtr##X; \
  typedef boost::shared_ptr<const X> SPtrConst##X; \
  inline SPtr##X create##SPtr##X(X &x) \
  { \
    boost::shared_ptr<X> px(&x, nullDeleter()); \
    return px; \
  } \
  inline SPtrConst##X create##SPtrConst##X(const X &x) \
  { \
    boost::shared_ptr<const X> px(&x, nullDeleter()); \
    return px; \
  } \
  NAME_SPACE_SPTR(X)

#define TYPEDEF_SAPTR(X) \
  typedef boost::shared_array<X> X##SAPtr ;\
  NAME_SPACE_SAPTR(X)

#define DEFINE_SPTR(X) \
  class X; \
  TYPEDEF_SPTR(X)


#define DEFINE_SAPTR(X) \
  class X; \
  TYPEDEF_SAPTR(X)


/* template with one argument */

#define NAME_SPACE_TPL1_SPTR(N,X,Y)              \
  namespace SharedPointer                        \
  {                                              \
    typedef SPtr##N N;                           \
  }                                              \
  namespace SharedPointerConst                   \
  {                                              \
    typedef SPtrConst##N N;                      \
  }

#define TYPEDEF_TPL1_SPTR(N,X,Y)                           \
  typedef boost::shared_ptr<X<Y> > SPtr##N;                \
  typedef boost::shared_ptr<const X<Y> > SPtrConst##N;     \
  inline SPtr##N create##SPtr##N(X<Y> &x)                  \
  {                                                        \
    boost::shared_ptr<X<Y> > px(&x, nullDeleter());        \
    return px;                                             \
  }                                                        \
  inline SPtrConst##N create##SPtrConst##N(const X<Y> &x)  \
  {                                                        \
    boost::shared_ptr<const X<Y> > px(&x, nullDeleter());  \
    return px;                                             \
  }                                                        \
  NAME_SPACE_TPL1_SPTR(N,X,Y)


/* *SPtr types definitions, outside class headers because of some race
   conditions */

DEFINE_SPTR(BlockCSRMatrix)

DEFINE_SPTR(DynamicalSystemXML)

DEFINE_SPTR(Interaction)
DEFINE_SPTR(InteractionXML)

DEFINE_SPTR(Model)


DEFINE_SPTR(NonSmoothDynamicalSystemXML)
DEFINE_SPTR(NonSmoothDynamicalSystem)

DEFINE_SPTR(NonSmoothLawXML)

DEFINE_SPTR(OneStepNSProblem)
DEFINE_SPTR(OneStepNSProblemXML)

DEFINE_SPTR(OneStepIntegrator)
DEFINE_SPTR(OneStepIntegratorXML)

DEFINE_SPTR(Relation)
DEFINE_SPTR(RelationXML)

DEFINE_SPTR(Simulation)
DEFINE_SPTR(SimulationXML)

DEFINE_SPTR(LCP)

DEFINE_SPTR(MixedComplementarityConditionNSL)
TYPEDEF_SPTR(MixedLinearComplementarityProblem)

DEFINE_SPTR(TimeDiscretisationXML)
DEFINE_SPTR(TimeDiscretisation)

TYPEDEF_SAPTR(integer)
TYPEDEF_SPTR(integer)
TYPEDEF_SAPTR(doublereal)
TYPEDEF_SPTR(doublereal)

#endif /* SiconosPointers_hpp */
