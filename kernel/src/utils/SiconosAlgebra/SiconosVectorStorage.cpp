#include "SiconosVectorStorage.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/blas.hpp>

#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>




namespace bindings = boost::numeric::bindings::blas;

DenseVectStorage::DenseVectStorage() : internal_data() {};

DenseVectStorage::DenseVectStorage(size_t r) : internal_data(r) {};

DenseVectStorage::DenseVectStorage(size_t r, double val) : internal_data(r)
{
  bindings::set(val, internal_data);
};

SparseVectStorage::SparseVectStorage() : internal_data() {};

SparseVectStorage::SparseVectStorage(size_t r) : internal_data(r) {};
