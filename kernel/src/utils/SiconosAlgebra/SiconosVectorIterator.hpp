/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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

/*! \file SiconosVector.hpp
 */

#ifndef __SiconosVectorIterator__
#define __SiconosVectorIterator__

#include <iterator>

/** Iterator for SiconosVector covering both possible types. */
template< typename V, typename T >
struct SiconosVectorIteratorType
  : public std::iterator<std::forward_iterator_tag, T>
{
  V* v;
  size_t p;
  SiconosVectorIteratorType() : v(0), p(0) {};
  SiconosVectorIteratorType(V& _v, size_t _p) : v(&_v), p(_p) {};
  bool operator!=(const SiconosVectorIteratorType& it)
    { return v != it.v || p != it.p; }
  bool operator==(const SiconosVectorIteratorType& it)
    { return v == it.v && p == it.p; }
  SiconosVectorIteratorType& operator=(const SiconosVectorIteratorType& it)
    { v=it.v; p=it.p; return *this; }
  SiconosVectorIteratorType& operator++()
    { if (p<v->size()) p++; return *this; }
  SiconosVectorIteratorType operator++(int) {
    SiconosVectorIteratorType tmp(*this);
    if (p < v->size()) p++;
    return tmp;
  }
  T& operator*() { T &d((*v)(p)); return d; }
};

/* Note: Derived classes here instead of typedefs only because they
 *       are easier to deal with in SWIG */

/** Specialization for non-const SiconosVector */
struct SiconosVectorIterator
  : public SiconosVectorIteratorType<SiconosVector, double>
{
  SiconosVectorIterator() : SiconosVectorIteratorType<SiconosVector, double>() {}
  SiconosVectorIterator(SiconosVectorIteratorType<SiconosVector, double>& it)
    : SiconosVectorIteratorType<SiconosVector, double>(*it.v,it.p) {}
  SiconosVectorIterator(SiconosVector& _v, size_t _p)
    : SiconosVectorIteratorType<SiconosVector, double>(_v,_p) {}
};

/** Specialization for const SiconosVector */
struct SiconosVectorConstIterator
  : public SiconosVectorIteratorType<const SiconosVector, const double>
{
  SiconosVectorConstIterator()
    : SiconosVectorIteratorType<const SiconosVector, const double>() {}
  SiconosVectorConstIterator(
    SiconosVectorIteratorType<const SiconosVector, const double>& it)
    : SiconosVectorIteratorType<const SiconosVector, const double>(*it.v,it.p) {}
  SiconosVectorConstIterator(const SiconosVector& _v, size_t _p)
    : SiconosVectorIteratorType<const SiconosVector, const double>(_v,_p) {}
};

#endif
