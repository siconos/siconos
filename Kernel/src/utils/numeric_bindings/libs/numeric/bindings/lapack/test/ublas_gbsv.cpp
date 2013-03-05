
#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/bindings/ublas/banded.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/bandwidth.hpp>
#include <boost/numeric/bindings/lapack/computational/gbtrf.hpp>
#include <boost/numeric/bindings/lapack/computational/gbtrs.hpp>
#include <vector>
#include <stdexcept>

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;

// solves the equation Ax = B, and puts the solution in B
// A is mutated by this routine
template <typename MatrA, typename MatrB>
void InPlaceSolve(MatrA& a, MatrB& b)
{
  // if the matrix has kl lower and ku upper diagonals, then we should have
  // allocated kl lower and kl+ku upper diagonals
  fortran_int_t const kl = bindings::bandwidth_lower (a);
  fortran_int_t const ku = bindings::bandwidth_upper (a) - kl;
  std::vector<fortran_int_t> piv(a.size1());
  int ret = lapack::gbtrf(/*kl, ku, */a, piv);
  if (ret < 0) {
    //CStdString err;
    //err.Format("banded::Solve: argument %d in DGBTRF had an illegal value", -ret);
    //throw RuntimeError(err);
    throw std::runtime_error("banded::Solve: argument %d in DGBTRF had an illegal value");
  }
  if (ret > 0) {
    //CStdString err;
    //err.Format("banded::Solve: the (%d,%d) diagonal element is 0 after DGBTRF", ret, ret);
    //throw RuntimeError(err);
    throw std::runtime_error("banded::Solve: the (%d,%d) diagonal element is 0 after DGBTRF");
  }

  ret = lapack::gbtrs(/*kl, ku, */a, piv, b);
  if (ret < 0) {
    //CStdString err;
    //err.Format("banded::Solve: argument %d in DGBTRS had an illegal value", -ret);
    //throw RuntimeError(err);
    throw std::runtime_error("banded::Solve: argument %d in DGBTRS had an illegal value");
  }
}

template<typename T>
void do_typename()
{
  using namespace boost::numeric::ublas;
  // if the matrix has kl lower and ku upper diagonals, then we should
  // allocate kl lower and kl+ku upper diagonals
  size_t sz = 1000, kl = 1, ku = 1;
  ublas::banded_matrix<T> a(sz, sz, kl, kl+ku);
  ublas::vector<T> b(sz);
  // fill values in a and b
  for (size_t i = 0; i < sz; ++i) {
    a(i,i) = i;
    b(i) = i;
  }
  for (size_t i = 1; i < sz; ++i) {
    a(i,i-1) = i;
    a(i-1,i) = 1;
  }
  InPlaceSolve(a, b);
}

int main()
{
  do_typename<float>();
  do_typename<double>();
  do_typename<std::complex<float> >();
  do_typename<std::complex<double> >();
  return 0;
}
