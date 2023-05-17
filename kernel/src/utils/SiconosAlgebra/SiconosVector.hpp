/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifndef __SiconosVector__
#define __SiconosVector__

// #include <boost/numeric/ublas/fwd.hpp>
// #include <boost/numeric/ublas/storage.hpp>
// #include <boost/numeric/ublas/vector.hpp>
#include "EigenInclude.hpp"
#include <Eigen/Core>
// #include <vector>

// #include "SiconosAlgebraTypes.hpp"   // For UblasType
// #include "SiconosSerialization.hpp"  // For ACCEPT_SERIALIZATION
// #include "SiconosVectorIterator.hpp"

namespace siconos::algebra {

// struct SiconosVectorIterator;
// struct SiconosVectorConstIterator;
// class BlockVector;
// class SimpleMatrix;

// adaptor to allow construction of a boost vector from memory without
// copy
// Source
// https://stackoverflow.com/questions/1735841/initializing-a-ublas-vector-from-a-c-array
// Example:
//
// using vector_adaptor = boost::numeric::ublas::vector<double,shallow_array_adaptor<double> >;
// double a[size];
// vector_adaptor v(shallow_array_adaptor<double>(size, &tab[0]));
// std::vector w(size)
// vector_adaptor v2(shallow_array_adaptor<double>(size, w.data()));
//
// Experimental: works only if  BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR is defined (see
// kernel/CMakeLists.txt)
// and breaks many things in bindings, ublas ...
//
// template <typename T>
// class shallow_array_adaptor : public boost::numeric::ublas::shallow_array_adaptor<T> {
//  public:
//   typedef boost::numeric::ublas::shallow_array_adaptor<T> base_type;
//   typedef typename base_type::size_type size_type;
//   typedef typename base_type::pointer pointer;

//   shallow_array_adaptor(size_type n) : base_type(n) {}
//   shallow_array_adaptor(size_type n, pointer data) : base_type(n, data) {}
//   shallow_array_adaptor(const shallow_array_adaptor &c) : base_type(c) {}

//   void swap(shallow_array_adaptor &a)
//   {
//     if (base_type::begin() != a.begin())
//       std::swap_ranges(base_type::begin(), base_type::end(), a.begin());
//   }
// };


/**
   Vectors of double. (Interface to various types of Boost-Ublas vectors).

   Two possible types: siconos::algebra::DENSE (default) and Siconos:SPARSE.

*/

using SiconosVector = Eigen::VectorXd;

void concatenateVectors(SiconosVector& target, const SiconosVector& a, const SiconosVector& b);

}  // namespace siconos::algebra

#endif