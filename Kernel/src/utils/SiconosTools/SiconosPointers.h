#ifndef SiconosPointers_h
#define SiconosPointers_h

#ifdef WithSmartPtr

#include <boost/shared_ptr.hpp>

#define DEFINE_SPTR(X) \
  class X; \
  typedef boost::shared_ptr<X> X##SPtr

#else /* WithSmartPtr */

#define DEFINE_SPTR(X) \
  class X; \
  typedef X * X##SPtr

#endif /* WithSmartPtr */

DEFINE_SPTR(SiconosVector);
DEFINE_SPTR(SiconosMatrix);
DEFINE_SPTR(OSNSMatrix);
DEFINE_SPTR(NonSmoothLaw);
DEFINE_SPTR(MixedComplementarityConditionNSL);

#endif /* SiconosPointers_h */
