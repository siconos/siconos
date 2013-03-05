
// solving A * X = B
// A symmetric
// driver function sysv()

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/driver/sysv.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cin;
using std::cout;
using std::endl; 

typedef double real_t; 
typedef std::complex<real_t> cmplx_t; 

typedef ublas::matrix<real_t, ublas::column_major> m_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> cm_t;

typedef ublas::symmetric_adaptor<m_t, ublas::lower> symml_t; 
typedef ublas::symmetric_adaptor<m_t, ublas::upper> symmu_t; 

typedef ublas::symmetric_adaptor<cm_t, ublas::lower> csymml_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::upper> csymmu_t; 

template <typename M>
void init_symm2 (M& m) {
  for (int i = 0; i < m.size1(); ++i) 
    for (int j = i; j < m.size1(); ++j)
      m (i, j) = m (j, i) = 1 + j - i; 
}

int main (int argc, char **argv) {
  size_t n = 0;
  if (argc > 1) {
    n = atoi(argv [1]);
  }

  cout << endl; 

  // symmetric 
  cout << "real symmetric\n" << endl; 

  if (n <= 0) {
    cout << "n -> ";
    cin >> n;
  }
  if (n < 5) n = 5; 
  cout << "min n = 5" << endl << endl; 
  size_t nrhs = 2; 
  m_t al (n, n), au (n, n);  // matrices (storage)
  symml_t sal (al);   // symmetric adaptor
  symmu_t sau (au);   // symmetric adaptor
  m_t x (n, nrhs);
  m_t bl (n, nrhs), bu (n, nrhs);  // RHS matrices

  std::vector<fortran_int_t> ipiv (n);

  init_symm2 (al); 
  swap (row (al, 1), row (al, 4)); 
  swap (column (al, 1), column (al, 4)); 

  print_m (al, "al"); 
  cout << endl; 

  init_symm2 (au); 
  swap (row (au, 2), row (au, 3)); 
  swap (column (au, 2), column (au, 3)); 

  print_m (au, "au"); 
  cout << endl; 

  for (int i = 0; i < x.size1(); ++i) {
    x (i, 0) = 1.;
    x (i, 1) = 2.; 
  }
  bl = prod (sal, x); 
  bu = prod (sau, x); 

  print_m (bl, "bl"); 
  cout << endl; 
  print_m (bu, "bu"); 
  cout << endl; 

  m_t al1 (al), au1 (au);  // for part 2
  m_t bl1 (bl), bu1 (bu); 

//  lapack::sysv (sal, bl);
//  no ipiv less version is currently provided, so fall back to using ipiv
  lapack::sysv (sal, ipiv, bl);
  print_m (bl, "xl"); 
  cout << endl; 

//  lapack::sysv (sau, bu);
//  no ipiv less version is currently provided, so fall back to using ipiv
  lapack::sysv (sau, ipiv, bu);
  print_m (bu, "xu"); 
  cout << endl; 

  // part 2 

  std::vector<real_t> work (1); 

  symml_t sal1 (al1);
  int err = lapack::sysv (sal1, ipiv, bl1, lapack::workspace(work));
  print_m (al1, "al1 factored"); 
  cout << endl; 
  print_v (ipiv, "ipiv"); 
  cout << endl; 
  print_m (bl1, "xl1"); 
  cout << endl; 

  symmu_t sau1 (au1);
  err = lapack::sysv (sau1, ipiv, bu1, lapack::workspace(work));
  print_m (au1, "au1 factored"); 
  cout << endl; 
  print_v (ipiv, "ipiv"); 
  cout << endl; 
  print_m (bu1, "xu1"); 
  cout << endl; 
  cout << endl; 

  //////////////////////////////////////////////////////////
  cout << "\n==========================================\n" << endl; 
  cout << "complex symmetric\n" << endl; 

  cm_t cal (n, n), cau (n, n);   // matrices (storage)
  csymml_t scal (cal);   // symmetric adaptor 
  csymmu_t scau (cau);   // symmetric adaptor 
  cm_t cx (n, 1); 
  cm_t cbl (n, 1), cbu (n, 1);  // RHS

  init_symm2 (cal); 
  init_symm2 (cau); 
  cal *= cmplx_t (1, 1); 
  cau *= cmplx_t (1, -0.5); 

  print_m (cal, "cal"); 
  cout << endl; 
  print_m (cau, "cau"); 
  cout << endl; 

  for (int i = 0; i < cx.size1(); ++i) 
    cx (i, 0) = cmplx_t (1, -1); 
  print_m (cx, "cx"); 
  cout << endl; 
  cbl = prod (scal, cx);
  cbu = prod (scau, cx);
  print_m (cbl, "cbl"); 
  cout << endl; 
  print_m (cbu, "cbu"); 
  cout << endl; 

//  int ierr = lapack::sysv (scal, cbl);
//  no ipiv less version is currently provided, so fall back to using ipiv
  int ierr = lapack::sysv (scal, ipiv, cbl);
  if (ierr == 0)
    print_m (cbl, "cxl"); 
  else 
    cout << "matrix is not regular: ierr = " 
         << ierr << endl;
  cout << endl; 

  std::vector<cmplx_t> cwork (n); 

  ierr = lapack::sysv (scau, ipiv, cbu, lapack::workspace(cwork));
  if (ierr == 0) {
    print_v (ipiv, "ipiv"); 
    cout << endl; 
    print_m (cbu, "cxu"); 
  }
  else 
    cout << "matrix is not regular: ierr = " 
         << ierr << endl;
  cout << endl; 

}

