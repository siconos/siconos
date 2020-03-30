/*hseqr.cpp
 */

#include<iostream>
#include<complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
//#include <boost/numeric/bindings/detail/complex_utils.hpp>
#include <boost/numeric/bindings/vector_view.hpp>
#include <boost/numeric/bindings/lapack/computational/hseqr.hpp>
#include <boost/numeric/bindings/lapack/computational/trevc.hpp>

using std::cout;
using std::endl;
using std::vector;
using std::complex;

namespace ublas =  boost::numeric::ublas;
namespace lapack =  boost::numeric::bindings::lapack;
namespace bindings =  boost::numeric::bindings;
namespace tag =  boost::numeric::bindings::tag;

void hseqr(int);
template <typename T>
void Hessenberg(ublas::matrix<T, ublas::column_major>&);
template <typename T>
void Hessenberg(ublas::matrix<complex<T>, ublas::column_major>&);

int main()
{
  cout << "I'm testing uBlas." << endl;

  int n = 5;
  hseqr(n);

}

void hseqr(int n)
{
  cout << "\nCalculating eigenvalues using LAPACK's hseqr." << endl;
  ublas::matrix<double, ublas::column_major> H(n,n);
  Hessenberg(H);
  cout << "\nUpper Hessenberg matrix H:\n" << H << endl;

  ublas::vector<complex<double> > values(n);
  ublas::matrix<double, ublas::column_major> Z(n,n);

  cout << "\nHSEQR for only eigenvalues." << endl;
  ublas::matrix<double, ublas::column_major> Z_dummy(1,1);
  lapack::hseqr('E', 'N', 1, n, H, values, Z_dummy);
  /*
  lapack::hseqr('E', 'N', 1, n, H,
      bindings::detail::real_part_view(values),
      bindings::detail::imag_part_view(values),
      Z_dummy);
  bindings::detail::interlace(values);
  */
  cout << "\nH:\n" << H << endl;
  cout << "\nvalues: " << values << endl;

  cout << "\nHSEQR for eigenvalues and Schur vectors." << endl;
  Hessenberg(H);
  cout << "H:\n" << H << endl;
  lapack::hseqr('S', 'I', 1, n, H, values, Z);
  /*
  lapack::hseqr('S', 'I', 1, n, H,
      bindings::detail::real_part_view(values),
      bindings::detail::imag_part_view(values),
      Z);
  bindings::detail::interlace(values);
  */
  cout << "\nH: " << H << endl;
  cout << "\nvalues: " << values << endl;
  cout << "\nZ: " << Z << endl;

  cout << "\n==================================" << endl;
  cout << "Recalculating original matrix..." << endl;
  ublas::matrix<double, ublas::column_major> cH(n,n);
  cH = ublas::prod(H, ublas::herm(Z));
  cH = ublas::prod(Z, cH);
  cout << "'New' original matrix:\n" << cH << endl;
  cout << "==================================" << endl;


  cout << "\nHSEQR for only eigenvalues.  Complex version" << endl;
  ublas::matrix<complex<double>, ublas::column_major> G(n,n);
  Hessenberg(G);
  cout << "\nG:\n" << G << endl;
  ublas::matrix<complex<double>, ublas::column_major> cZ_dummy(1,1);
  lapack::hseqr('E', 'N', 1, n, G, values, cZ_dummy);
  cout << "\nG:\n" << G << endl;
  cout << "\nvalues: " << values << endl;

  cout << "\nHSEQR for eigenvalues and Schur vectors." << endl;
  Hessenberg(G);
  cout << "G:\n" << G << endl;
  ublas::matrix<complex<double>, ublas::column_major> cZ(n,n);
  lapack::hseqr('S', 'I', 1, n, G, values, cZ);
  cout << "\nG:\n " << G << endl;
  cout << "\nvalues: " << values << endl;
  cout << "\nZ:\n " << Z << endl;

  cout << "\n==================================" << endl;
  cout << "Recalculating original matrix..." << endl;
  ublas::matrix<complex<double>, ublas::column_major> origG(G);
  origG = ublas::prod(G, ublas::herm(cZ));
  origG = ublas::prod(cZ, origG);
  cout << "'New' original matrix:\n" << origG << endl;
  cout << "==================================" << endl;

  ublas::matrix<complex<double>, ublas::column_major> cVL(cZ);
  ublas::matrix<complex<double>, ublas::column_major> cVR(cZ);
  boost::numeric::bindings::detail::array<complex<double> > work_c(2*n);
  boost::numeric::bindings::detail::array<double> work_r(n);
  ublas::vector<fortran_bool_t> select_dummy(n);
  fortran_int_t m_info(n+1);
  lapack::trevc(tag::both(),'B',select_dummy,G,cVL,cVR,n,m_info,lapack::workspace(work_c,work_r));

  cout << "\n==================================" << endl;
  cout << "Testing left & right eigenvectors..." << endl;
  Hessenberg(G);
  cout << "Many 'zeros':" << endl;
  for(int i=0; i<n; ++i)
  {
    cout << (ublas::prod(G, column(cVR, i)) - values(i) * column(cVR, i)) << endl;
    cout << (ublas::prod(ublas::herm(G), column(cVL, i)) - conj(values(i)) * column(cVL, i)) << endl;
  }
  cout << "==================================" << endl;

  cout << "\n==================================" << endl;
  cout << "Verifying diagonal matrix..." << endl;
  ublas::matrix<complex<double>, ublas::column_major> cG(n,n,0);
  Hessenberg(G);
  G = ublas::prod(G, cVR);
  G = ublas::prod(ublas::herm(cVL), G);
  cout << "'diagonal' matrix:\n" << G << endl;
  cout << "==================================" << endl;

}


template <typename T>
void Hessenberg(ublas::matrix<T, ublas::column_major>& H)
{
  T k = 1;
  for(unsigned int i = 0; i < H.size1(); ++i)
  {
    for(unsigned int j = i; j <= H.size2(); ++j)
    {
      if(j > 0)
      {
        H(i,j-1) = k;
        k += 1;
      }
    }
  }
}

template <typename T>
void Hessenberg(ublas::matrix<complex<T>, ublas::column_major>& H)
{
  T k = 1.0;
  for(unsigned int i = 0; i < H.size1(); ++i)
  {
    for(unsigned int j = i; j <= H.size2(); ++j)
    {
      if(j > 0)
      {
        T real = k++;
        T imag = k++;
        H(i,j-1) = complex<T>(real, imag);
      }
    }
  }
}

