#include "SiconosVectorOperators.hpp"

#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace bindings = boost::numeric::bindings::blas;

template<>
void Fill::operator()<DenseVect> (DenseVect& v)
{
  bindings::set(this->param(), v);
}

template<>
void Zero::operator()<DenseVectStorage> (DenseVectStorage& dv)
{
  bindings::scal(0.0, dv.internal_data);
}

template<>
void Scal::operator()<DenseVectStorage> (DenseVectStorage& dv)
{
  bindings::scal(this->param(), dv.internal_data);
};

double inner_prod(DenseVect& v1, DenseVect& v2)
{
  return bindings::dot(v1, v2);
}

template<>
void OpCopy::operator()(const DenseVect& dv2, DenseVect& dv1)
{
  dv1.resize(dv2.size());
  bindings::copy(dv2, dv1);
}

template<>
void OpAxpby::operator()<DenseVect, DenseVect>(const DenseVect& x,
                                               const double& a,
                                               const double& b,
                                               DenseVect& y)
{
  bindings::scal(b, y);
  bindings::axpy(a, x, y);
}

template<>
void OpAxpby::operator()<SparseVect, SparseVect>(const SparseVect& x,
                                                 const double& a,
                                                 const double& b,
                                                 SparseVect& y)
{
  y *= b;
  if (&y != &x)
  {
    noalias(y) += a*x;
  }
  else
  {
    y += a*x;
  }
}

template<>
void GetRValue::operator()<SparseVectStorage> (SparseVectStorage& storage)
{
  assert (this->param() <= storage.internal_data.size());

  answer = &(storage.internal_data(this->param()).ref());
}
