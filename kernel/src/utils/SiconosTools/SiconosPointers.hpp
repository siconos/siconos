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

#include <SiconosConfig.h>

// Case 1 : ref == shared pointers from c++ (>=11) standard
// SICONOS_STD_SHARED_PTR is automatically set by cmake (CXXCompilerSetup.cmake)
// while SICONOS_USE_BOOST_FOR_CXX11 is a user option.
#if defined(SICONOS_STD_SHARED_PTR) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#include <memory>
namespace std11 = std;
#else
// Case 2 : ref == boost shared pointers
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
namespace std11 = boost;
#endif

namespace SP {}
namespace SPC {}
namespace SA {}

/** \namespace SP Namespace for Siconos smart pointers : memory
    pointed by a Siconos smart pointers is automaticaly
    deallocated */
namespace SharedPointer = SP;

/** \namespace SA Namespace for Siconos shared arrays : Siconos shared
    arrays are automaticaly deallocated */
namespace SharedArray = SA;

/** \namespace SPC Namespace for const shared pointers : memory
    pointed by a const shared pointers is automaticaly deallocated */
namespace SharedPointerConst = SPC;



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
  namespace SP \
  { \
    typedef SPtr##X X; \
  } \
  namespace SPC \
  { \
    typedef SPtrConst##X X;\
  }

#define NAME_SPACE_SAPTR(X)                     \
  namespace SA \
  { \
    typedef X##SAPtr X; \
  }


/* template typedef : no */

#define TYPEDEF_SPTR(X) \
  typedef std11::shared_ptr<X> SPtr##X; \
  typedef std11::shared_ptr<const X> SPtrConst##X; \
  inline SPtr##X create##SPtr##X(X &x) \
  { \
    std11::shared_ptr<X> px(&x, nullDeleter()); \
    return px; \
  } \
  inline SPtrConst##X create##SPtrConst##X(const X &x) \
  { \
    std11::shared_ptr<const X> px(&x, nullDeleter()); \
    return px; \
  } \
  NAME_SPACE_SPTR(X)

// boost shared_arrays : at the time required
// only in HEM5 and LSodar
#include <boost/shared_array.hpp>
#define TYPEDEF_SAPTR(X)                        \
  typedef boost::shared_array<X> X##SAPtr ;     \
  NAME_SPACE_SAPTR(X)


#define DEFINE_SPTR(X)                          \
  class X; \
  TYPEDEF_SPTR(X)

#define DEFINE_SPTR_STRUCT(X) \
  struct X; \
  TYPEDEF_SPTR(X)

#define DEFINE_SAPTR(X) \
  class X; \
  TYPEDEF_SAPTR(X)


/* template with one argument */

#define NAME_SPACE_TPL1_SPTR(N,X,Y)              \
  namespace SP                                   \
  {                                              \
    typedef SPtr##N N;                           \
  }                                              \
  namespace SPC                                  \
  {                                              \
    typedef SPtrConst##N N;                      \
  }

#define TYPEDEF_TPL1_SPTR(N,X,Y)                           \
  typedef std11::shared_ptr<X<Y> > SPtr##N;                \
  typedef std11::shared_ptr<const X<Y> > SPtrConst##N;     \
  inline SPtr##N create##SPtr##N(X<Y> &x)                  \
  {                                                        \
    std11::shared_ptr<X<Y> > px(&x, nullDeleter());        \
    return px;                                             \
  }                                                        \
  inline SPtrConst##N create##SPtrConst##N(const X<Y> &x)  \
  {                                                        \
    std11::shared_ptr<const X<Y> > px(&x, nullDeleter());  \
    return px;                                             \
  }                                                        \
  NAME_SPACE_TPL1_SPTR(N,X,Y)

#endif /* SiconosPointers_hpp */
