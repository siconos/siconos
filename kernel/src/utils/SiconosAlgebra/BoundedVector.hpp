#ifndef BoundedVector_hpp
#define BoundedVector_hpp

#include "SiconosVector.hpp"
#include "SiconosVectorStorage.hpp"

#include "SiconosVectorOperators.hpp"

#define MAKE_BOUNDED_VECTOR(N)                                \
  class SiconosVector##N : public SiconosVector               \
  {                                                           \
  public:                                                     \
    SiconosVector##N() : SiconosVector(_bounded_storage) {};            \
                                                                        \
    SiconosVector##N(const SiconosVector& svect) : SiconosVector(_bounded_storage) \
    {                                                                   \
      apply_visitor<Copy>(storage(svect), storage(*this));              \
    };                                                                  \
    SiconosVector##N(const SiconosVector##N& svect) : SiconosVector(_bounded_storage) \
    {                                                                   \
      apply_visitor<Copy>(storage(svect), storage(*this));              \
    };                                                                  \
                                                                        \
    ~SiconosVector##N() { this->_storage = NULL; };                     \
                                                                        \
  protected:                                                            \
    ACCEPT_SERIALIZATION(SiconosVector);                                \
                                                                        \
    BoundedVectStorage##N _bounded_storage;                             \
                                                                        \
  }

MAKE_BOUNDED_VECTOR(3);
MAKE_BOUNDED_VECTOR(4);
MAKE_BOUNDED_VECTOR(6);
MAKE_BOUNDED_VECTOR(7);

#endif
