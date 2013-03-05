/*hseqr.cpp
 */

#include<iostream>
#include<complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/lapack/geev.hpp>
#include <boost/numeric/bindings/lapack/hseqr.hpp>
#include <boost/numeric/bindings/lapack/trevc.hpp>

using std::cout;
using std::endl;
using std::vector;
using std::complex;

namespace ublas =  boost::numeric::ublas;
namespace lapack =  boost::numeric::bindings::lapack;

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
  ublas::matrix<double, ublas::column_major> H(n, n);
  Hessenberg(H);
  cout << "\nUpper Hessenberg matrix H:\n" << H << endl;

  ublas::vector<complex<double> > values(n);
  ublas::matrix<double, ublas::column_major> Z(n, n);

  cout << "\nHSEQR for only eigenvalues." << endl;
  lapack::hseqr('E', H, values);
  cout << "\nH:\n" << H << endl;
  cout << "\nvalues: " << values << endl;

  cout << "\nHSEQR for eigenvalues and Schur vectors." << endl;
  Hessenberg(H);
  cout << "H:\n" << H << endl;
  lapack::hseqr('S', 'I', H, values, Z);
  cout << "\nH: " << H << endl;
  cout << "\nvalues: " << values << endl;
  cout << "\nZ: " << Z << endl;

  cout << "\n==================================" << endl;
  cout << "Recalculating original matrix..." << endl;
  ublas::matrix<double, ublas::column_major> cH(n, n);
  cH = ublas::prod(H, ublas::herm(Z));
  cH = ublas::prod(Z, cH);
  cout << "'New' original matrix:\n" << cH << endl;
  cout << "==================================" << endl;


  cout << "\nHSEQR for only eigenvalues.  Complex version" << endl;
  ublas::matrix<complex<double>, ublas::column_major> G(n, n);
  Hessenberg(G);
  cout << "\nG:\n" << G << endl;
  lapack::hseqr('E', G, values);
  cout << "\nG:\n" << G << endl;
  cout << "\nvalues: " << values << endl;

  cout << "\nHSEQR for eigenvalues and Schur vectors." << endl;
  Hessenberg(G);
  cout << "G:\n" << G << endl;
  ublas::matrix<complex<double>, ublas::column_major> cZ(n, n);
  lapack::hseqr('S', 'I', G, values, cZ);
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

}


template <typename T>
void Hessenberg(ublas::matrix<T, ublas::column_major>& H)
{
  T k = 1;
  for (unsigned int i = 0; i < H.size1(); ++i)
  {
    for (unsigned int j = i; j <= H.size2(); ++j)
    {
      if (j > 0)
      {
        H(i, j - 1) = k;
        k += 1;
      }
    }
  }
}

template <typename T>
void Hessenberg(ublas::matrix<complex<T>, ublas::column_major>& H)
{
  T k = 1.0;
  for (unsigned int i = 0; i < H.size1(); ++i)
  {
    for (unsigned int j = i; j <= H.size2(); ++j)
    {
      if (j > 0)
      {
        T real = k++;
        T imag = k++;
        H(i, j - 1) = complex<T>(real, imag);
      }
    }
  }
}

