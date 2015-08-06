//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cstdlib>
#include <complex>
#include <assert.h>

namespace Siconos {
  namespace algebra {

    template <typename T>
    T random_value() {
      assert( false );
      return 0;
    }

    template <>
    float random_value<float>() {return float(std::rand()) / float(RAND_MAX) - 0.5f;}

    template <>
    double random_value<double>() {return double(std::rand()) / double(RAND_MAX) - 0.5;}

    template <>
    std::complex<float> random_value< std::complex<float> >() {return std::complex<float>(random_value<float>(), random_value<float>() );}

    template <>
    std::complex<double> random_value< std::complex<double> >() {return std::complex<double>(random_value<double>(), random_value<double>() );}


  } // Algebra namespace
} // Siconos namespace
#endif
