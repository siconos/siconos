#ifndef BoundedVector_hpp
#define BoundedVector_hpp

#include "SiconosVector.hpp"
#include "SiconosVectorStorage.hpp"

#include "SiconosVectorOperators.hpp"

template<size_t N>
class BoundedVector : public SiconosVector
{
public:
  BoundedVector() : SiconosVector(_bounded_storage) {};

  BoundedVector(const SiconosVector& svect) : SiconosVector(_bounded_storage)
  {
    apply_visitor<Copy>(storage(svect), storage(*this));
  };
  BoundedVector(const BoundedVector<N>& svect) : SiconosVector(_bounded_storage)
  {
    apply_visitor<Copy>(storage(svect), storage(*this));
  };

  ~BoundedVector() { this->_storage = NULL; };


protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosVector);

  BoundedVectStorage<N> _bounded_storage;

};

typedef BoundedVector<3> SiconosVector3;
typedef BoundedVector<4> SiconosVector4;
typedef BoundedVector<6> SiconosVector6;
typedef BoundedVector<7> SiconosVector7;

#endif
