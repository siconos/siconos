#ifndef SiconosPointers_h
#define SiconosPointers_h

/* *SPtr types are smart pointers */

/* http://www.boost.org/doc/libs/release/libs/smart_ptr */

#include <SiconosNumerics.h>

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace SharedPointer {};
namespace SharedPointerConst {};
namespace SharedArray {};

/* Using a shared_ptr to hold a pointer to a statically allocated
   object */
/* use create<type>SPtr(<type> &x)
/* cf http://www.boost.org/doc/ */
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

/* *SPtr types definitions */

DEFINE_SPTR(SiconosVector);
DEFINE_SPTR(SimpleVector);
DEFINE_SPTR(BlockVector);
DEFINE_SPTR(SiconosMatrix);
DEFINE_SPTR(SimpleMatrix);
DEFINE_SPTR(BlockMatrix);
DEFINE_SPTR(OSNSMatrix);
DEFINE_SPTR(NonSmoothLaw);
DEFINE_SPTR(MixedComplementarityConditionNSL);
DEFINE_SPTR(NonSmoothDynamicalSystem);
DEFINE_SPTR(SiconosMemory);
DEFINE_SPTR(DynamicalSystem);
DEFINE_SPTR(InteractionXML);
DEFINE_SPTR(RelationXML);
DEFINE_SPTR(NonSmoothDynamicalSystemXML)
DEFINE_SPTR(Topology);
DEFINE_SPTR(Simulation);
DEFINE_SPTR(SiconosModelXML);
DEFINE_SPTR(SimulationXML);
DEFINE_SPTR(Model);
DEFINE_SPTR(TimeDiscretisation);
DEFINE_SPTR(OneStepNSProblem);
DEFINE_SPTR(TimeStepping);
DEFINE_SPTR(OneStepIntegrator);
DEFINE_SPTR(EventsManager);
DEFINE_SPTR(UnitaryRelation);
DEFINE_SPTR(DynamicalSystemXML);
DEFINE_SPTR(FirstOrderLinearDS);
DEFINE_SPTR(FirstOrderLinearTIDS);
DEFINE_SPTR(FirstOrderLinearDSXML);
DEFINE_SPTR(NonSmoothLawXML);
DEFINE_SPTR(LagrangianLinearTIDS);
DEFINE_SPTR(LagrangianLinearTIDSXML);
DEFINE_SPTR(LagrangianDS);
DEFINE_SPTR(LagrangianDSXML);
DEFINE_SPTR(FirstOrderNonLinearDS);
DEFINE_SPTR(FirstOrderNonLinearDSXML);
DEFINE_SPTR(NewtonImpactFrictionNSL);
DEFINE_SPTR(NewtonImpactFrictionNSLXML);
DEFINE_SPTR(SparseBlockMatrix);
DEFINE_SPTR(OneStepIntegratorXML);
DEFINE_SPTR(Lsodar);
DEFINE_SPTR(Moreau);
DEFINE_SPTR(MoreauXML);
DEFINE_SPTR(SiconosMemoryXML);
DEFINE_SPTR(OneStepNSProblemXML);
DEFINE_SPTR(NonSmoothSolver);
DEFINE_SPTR(TimeDiscretisationXML);
DEFINE_SPTR(EventDriven);
DEFINE_SPTR(FrictionContact);
DEFINE_SPTR(FrictionContactXML);
DEFINE_SPTR(Sensor)
DEFINE_SPTR(Actuator)
DEFINE_SPTR(Event);
DEFINE_SPTR(LsodarXML);
DEFINE_SPTR(QPXML);
DEFINE_SPTR(LCPXML);
DEFINE_SPTR(Interaction);
DEFINE_SPTR(Relation);

TYPEDEF_SPTR(MixedLinearComplementarity_Problem);
TYPEDEF_SPTR(Numerics_Options);
TYPEDEF_SPTR(NumericsMatrix);
TYPEDEF_SPTR(SparseBlockStructuredMatrix);

TYPEDEF_SAPTR(integer);
TYPEDEF_SPTR(integer);
TYPEDEF_SAPTR(integer);
TYPEDEF_SPTR(doublereal);
TYPEDEF_SAPTR(doublereal);

namespace SP = SharedPointer;
namespace SA = SharedArray;
namespace SPC = SharedPointerConst;

#endif /* SiconosPointers_h */
