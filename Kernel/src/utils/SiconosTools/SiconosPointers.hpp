#ifndef SiconosPointers_hpp
#define SiconosPointers_hpp

/* *SPtr types are smart pointers */

/* http://www.boost.org/doc/libs/release/libs/smart_ptr */

#include <SiconosNumerics.h>

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace SharedPointer {};
namespace SharedPointerConst {};
namespace SharedArray {};

/** Using a shared_ptr to hold a pointer to a statically allocated
   object
   use create<type>SPtr(<type> &x)
   cf http://www.boost.org/doc/
*/
struct nullDeleter
{
  void operator()(void const *) const {}
};

/* template namespace : no */

#define NAME_SPACE_SPTR(X) \
  namespace SharedPointer \
  { \
    typedef SPtr##X X; \
  }; \
  namespace SharedPointerConst \
  { \
    typedef SPtrConst##X X;\
  }

#define NAME_SPACE_SAPTR(X)                     \
  namespace SharedArray \
  { \
    typedef X##SAPtr X; \
  }


/* template typedef : no */

#define TYPEDEF_SPTR(X) \
  typedef boost::shared_ptr<X> SPtr##X; \
  typedef boost::shared_ptr<const X> SPtrConst##X; \
  inline SPtr##X create##SPtr##X(X &x) \
  { \
    boost::shared_ptr<X> px(&x, nullDeleter()); \
    return px; \
  }; \
  inline SPtrConst##X create##SPtrConst##X(const X &x) \
  { \
    boost::shared_ptr<const X> px(&x, nullDeleter()); \
    return px; \
  } ;\
  NAME_SPACE_SPTR(X)


#define TYPEDEF_SAPTR(X) \
  typedef boost::shared_array<X> X##SAPtr ;\
  NAME_SPACE_SAPTR(X)

#define DEFINE_SPTR(X) \
  class X; \
  TYPEDEF_SPTR(X)


#define DEFINE_SAPTR(X) \
  class X; \
  TYPEDEF_SAPTR(X)

/* *SPtr types definitions, outside class headers because of some race
   conditions */

DEFINE_SPTR(BlockCSRMatrix);

DEFINE_SPTR(DynamicalSystemXML);

DEFINE_SPTR(Interaction);
DEFINE_SPTR(InteractionXML);

DEFINE_SPTR(Model);


DEFINE_SPTR(NonSmoothDynamicalSystemXML);
DEFINE_SPTR(NonSmoothDynamicalSystem);

DEFINE_SPTR(NonSmoothLawXML);

DEFINE_SPTR(OneStepNSProblem);
DEFINE_SPTR(OneStepNSProblemXML);

DEFINE_SPTR(OneStepIntegrator);
DEFINE_SPTR(OneStepIntegratorXML);

DEFINE_SPTR(Relation);
DEFINE_SPTR(RelationXML);

DEFINE_SPTR(Simulation);
DEFINE_SPTR(SimulationXML);


DEFINE_SPTR(MixedComplementarityConditionNSL);
TYPEDEF_SPTR(MixedLinearComplementarity_Problem);

DEFINE_SPTR(TimeDiscretisationXML);

TYPEDEF_SAPTR(integer);
TYPEDEF_SPTR(integer);
TYPEDEF_SAPTR(integer);
TYPEDEF_SPTR(doublereal);
TYPEDEF_SAPTR(doublereal);

namespace SP = SharedPointer;
namespace SA = SharedArray;
namespace SPC = SharedPointerConst;

#endif /* SiconosPointers_hpp */
