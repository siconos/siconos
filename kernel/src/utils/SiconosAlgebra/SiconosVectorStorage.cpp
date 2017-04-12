#include "SiconosVectorStorage.hpp"

namespace bindings = boost::numeric::bindings::blas;

DenseVectStorage::DenseVectStorage() : internal_data() {};

DenseVectStorage::DenseVectStorage(size_t r) : internal_data(r) {};

DenseVectStorage::DenseVectStorage(size_t r, double val) : internal_data(r)
{
  bindings::set(val, internal_data);
};


SparseVectStorage::SparseVectStorage() : internal_data() {};

SparseVectStorage::SparseVectStorage(size_t r) : internal_data(r) {};
