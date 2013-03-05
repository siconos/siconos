//
// Copyright (c) 2010 Thomas Klimpel
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <iostream>
#include <vector>
#include <boost/numeric/bindings/detail/complex_utils.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

int test_inshuffle(std::size_t n) {
  std::vector<std::size_t> data(2*n);
  for (std::size_t i = 0; i < n; ++i) {
    data[i] = 2*i+1;
    data[i+n] = 2*i;
  }
  boost::numeric::bindings::detail::inshuffle(&data[0],n);
  for (std::size_t i = 1; i < 2*n; ++i)
    if (data[i-1]+1 != data[i]) {
      std::cout << "Test inshuffle for n=" << n << std::endl;
      for (std::size_t j = 0; j < 2*n; ++j)
        std::cout << " " << data[j];
      std::cout << std::endl;
      std::cout << "logic error" << std::endl;
      return 255;
    }
  return 0;
}

int test_interlace(std::size_t n) {
  std::vector<std::complex<double> > data(n);
  for (std::size_t i = 0; i < n; ++i)
    data[i] = std::complex<double>(2*i,2*i+1);
  boost::numeric::bindings::detail::interlace(data);
  for (std::size_t i = 0; i < n; ++i)
    if (data[i].real() != i || data[i].imag() != i+n) {
      std::cout << "Test interlace for n=" << n << std::endl;
      for (std::size_t j = 0; j < n; ++j)
        std::cout << " " << data[j];
      std::cout << std::endl;
      std::cout << "logic error" << std::endl;
      return 255;
    }
  return 0;
}

int main() {
  for (std::size_t n = 1; n <= 1000; ++n) {
    int success = test_inshuffle(n);
    if (success != 0)
      return success;
  }
  for (std::size_t n = 1; n <= 1000; ++n) {
    int success = test_interlace(n);
    if (success != 0)
      return success;
  }
  return 0;
}
