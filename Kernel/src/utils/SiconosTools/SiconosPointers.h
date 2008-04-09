#ifndef SiconosPointers_h
#define SiconosPointers_h

#ifdef WithSmartPtr

#include <boost/shared_ptr.hpp>

#define DEFINE_SPTR(X) \
  typedef boost::shared_ptr<X> X##SPtr

#else /* WithSmartPtr */

#define DEFINE_SPTR(X) \
  typedef X * X##SPtr

#endif /* WithSmartPtr */

#include "SimpleVector.h"
DEFINE_SPTR(SiconosVector);

#include "SimpleMatrix.h"
DEFINE_SPTR(SiconosMatrix);

#include "OSNSMatrix.h"
DEFINE_SPTR(OSNSMatrix);

#include "NumericsMatrix.h"
DEFINE_SPTR(NumericsMatrix);


#endif /* SiconosPointers_h */
