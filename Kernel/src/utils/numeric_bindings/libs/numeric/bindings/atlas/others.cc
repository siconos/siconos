
// BLAS level 1
// std::vector<>, std::valarray<>, boost::array<>, C array

// element type: float or double
#ifdef F_FLOAT
typedef float real_t; 
#else
typedef double real_t; 
#endif

#include <iostream>
#include <iterator>
#include <algorithm>
#include <complex>

#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/bindings/std/valarray.hpp>
#include <boost/numeric/bindings/boost/array.hpp>
#include <boost/numeric/bindings/blas/level1.hpp>

#include "utils.h"

namespace blas = boost::numeric::bindings::blas;

using std::cout;
using std::endl;
using std::size_t; 

int main() {

  int n = 10; 

  cout << endl; 
  cout << "std::vector" << endl; 
  std::vector<real_t> sv (n); 
  init_v (sv, kpp (1)); 
  print_v (sv, "sv"); 
  cout << "std::valarray" << endl; 
  std::valarray<real_t> va (n); 
  blas::set (0.1, va); 
  print_v (va, "va"); 
  cout << endl; 

  cout << "dot(): sv^t va: ";
  real_t d = 0;
  for (int i = 0; i < n; ++i)
    d += sv[i] * va[i]; 

  cout << "is " << d << " == " << blas::dot (sv, va) << " ?" << endl; 
  cout << endl; 

#ifdef F_FLOAT
  cout << "sdsdot(): 10 + sv^T va = " << blas::sdsdot (10, sv, va) << endl; 
  cout << endl;
#endif  

  blas::scal (real_t(2), sv);
  print_v (sv, "scal(): 2 sv"); 

  cout << endl; 

  std::random_shuffle (sv.begin(), sv.end());
  cout << "shuffled sv: "; 
  std::copy (sv.begin(), sv.end(), std::ostream_iterator<real_t> (cout, " ")); 
  cout << endl; 
  int i = blas::iamax (sv); 
  cout << "iamax():\n  index of max el = " << i 
       << "; max el = " << sv[i] << endl; 
  cout << endl; 

  cout << "asum():\n  ||sv||_1 =  " << blas::asum (sv) 
       << "; ||va||_1 = " << blas::asum (va) << endl; 
  cout << "nrm2():\n  ||sv||_2 = " << blas::nrm2 (sv) 
       << "; ||va||_2 = " << blas::nrm2 (va) << endl; 
  cout << endl; 

  cout << "boost::array" << endl;
  boost::array<double, 10> ba;
  blas::set (0.1, ba);
  print_v (ba, "ba");
  cout << "C array" << endl; 
  typedef double double_array[10]; 
  double_array ca; 
  blas::set (1., ca); 
  print_v (ca, "ca");
  cout << endl; 
  
  blas::axpy (0.1, ba, ca); 
  print_v (ca, "axpy(): 0.1 ba + ca"); 

//  blas::axpby (0.1, ba, 2., ca); 
//  print_v (ca, "axpby(): 0.1 ba + 2.0 ca"); 

  cout << endl;
}
