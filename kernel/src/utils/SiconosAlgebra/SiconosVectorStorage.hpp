#ifndef SiconosVectorStorage_hpp
#define SiconosVectorStorage_hpp



#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosVisitor.hpp"

#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include <boost/array.hpp>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

struct SiconosVectorStorage
{
  ACCEPT_STD_VISITORS();

  virtual ~SiconosVectorStorage() {};

};

struct DenseVectStorage : public SiconosVectorStorage
{
  typedef DenseVect storage_type;
  DenseVect internal_data;
  DenseVectStorage();
  DenseVectStorage(size_t r);

  DenseVectStorage(size_t r, double val);
  DenseVectStorage(DenseVectStorage& st);

  virtual ~DenseVectStorage() {};

  ACCEPT_STD_VISITORS();
};


struct SparseVectStorage : public SiconosVectorStorage
{
  typedef SparseVect storage_type;
  SparseVect internal_data;
  SparseVectStorage();
  SparseVectStorage(size_t r);

  virtual ~SparseVectStorage() {};

  ACCEPT_STD_VISITORS();
};

#endif
