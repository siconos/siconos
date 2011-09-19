#define  BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR
#include <boost/numeric/ublas/vector.hpp>



#include "KernelTest.hpp"


#define DEBUG_MESSAGES 1
#include "../Register.hpp"



#include <SiconosKernel.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/typeof/typeof.hpp>

#define NVP(X) BOOST_SERIALIZATION_NVP(X)

CPPUNIT_TEST_SUITE_REGISTRATION(KernelTest);

BOOST_SERIALIZATION_ASSUME_ABSTRACT(SiconosVector)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(SiconosMatrix)

SICONOS_IO_REGISTER(SiconosMatrix,
                    (dimRow)
                    (dimCol)
                    (num))

SICONOS_IO_REGISTER_WITH_BASES(BlockMatrix, (SiconosMatrix),
                               (_mat)
                               (_tabRow)
                               (_tabCol))
SICONOS_IO_REGISTER_WITH_BASES(BlockVector, (SiconosVector),
                               (_sizeV)
                               (vect)
                               (_tabIndex))
SICONOS_IO_REGISTER_WITH_BASES(SensorPosition, (Sensor),
                               (_nSteps)
                               (_dataPlot)
                               (_k))
SICONOS_IO_REGISTER(NonSmoothDynamicalSystem,
                    (BVP)
                    (_topology)
                    (mIsLinear))
SICONOS_IO_REGISTER(Relation,
                    (_pluginh)
                    (_pluginJachx)
                    (_pluginJachlambda)
                    (_pluging)
                    (_pluginJacLg)
                    (_pluginf)
                    (_plugine)
                    (relationType)
                    (subType)
                    (_interaction)
                    (data)
                    (_workR)
                    (_workX)
                    (_workXdot)
                    (_workZ)
                    (_workY)
                    (_workL)
                    (_Residuy)
                    (_h_alpha)
                    (_jachlambda))
SICONOS_IO_REGISTER(ControlManager,
                    (_allSensors)
                    (_allActuators)
                    (_model))
SICONOS_IO_REGISTER_WITH_BASES(NewtonImpactNSL, (NonSmoothLaw),
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(NewtonEulerDS, (DynamicalSystem),
                               (_v)
                               (_v0)
                               (_vMemory)
                               (_qMemory)
                               (_fLMemory)
                               (_dotqMemory)
                               (_ndof)
                               (_qDim)
                               (_q)
                               (_deltaq)
                               (_q0)
                               (_dotq)
                               (_MObjToAbs)
                               (_I)
                               (_centerOfMass)
                               (_mass)
                               (_M)
                               (_luW)
                               (_T)
                               (_p)
                               (_mExt)
                               (_jacobianqmExt)
                               (_jacobianqDotmExt)
                               (_fExt)
                               (_NNL)
                               (_jacobianNNLq)
                               (_jacobianNNLqDot)
                               (_fL)
                               (_jacobianvFL)
                               (_jacobianqDotFL)
                               (computeFIntPtr)
                               (computeFExtPtr)
                               (computeMExtPtr)
                               (computeNNLPtr)
                               (computeJacobianFIntqPtr)
                               (computeJacobianFIntqDotPtr)
                               (computeJacobianNNLqPtr)
                               (computeJacobianNNLqDotPtr))
SICONOS_IO_REGISTER_WITH_BASES(PrimalFrictionContact, (LinearOSNS),
                               (_contactProblemDim)
                               (_sizeLocalOutput)
                               (_localVelocity)
                               (_localReaction)
                               (_tildeLocalVelocity)
                               (H)
                               (_mu)
                               (primalFrictionContact_driver))
SICONOS_IO_REGISTER(OSNSMatrix,
                    (dimRow)
                    (dimColumn)
                    (storageType)
                    (DSBlocksPositions)
                    //  (numericsMat)
                    (M1)
                    (Mt)
                    (M2)
                   )
SICONOS_IO_REGISTER_WITH_BASES(RelayNSL, (NonSmoothLaw),
                               (_lb)
                               (_ub))
SICONOS_IO_REGISTER_WITH_BASES(MixedComplementarityConditionNSL, (NonSmoothLaw),
                               (EqualitySize))
SICONOS_IO_REGISTER_WITH_BASES(Lsodar, (OneStepIntegrator),
                               (intData)
                               (rtol)
                               (atol)
                               (rwork)
                               (iwork)
                               (jroot)
                               (xWork))
SICONOS_IO_REGISTER(UnitaryRelation,
                    (_mainInteraction)
                    (_relativePosition)
                    (_number)
                    (_absolutePosition)
                    (_absolutePositionProj)
                    (_workX)
                    (_workXq)
                    (_workFree)
                    (_workYp)
                    (_workZ))
SICONOS_IO_REGISTER_WITH_BASES(SensorEvent, (Event),
                               (_sensor))
SICONOS_IO_REGISTER(SiconosSharedLibrary,
                    (isPlugged))
SICONOS_IO_REGISTER(OneStepIntegrator,
                    (integratorType)
                    (OSIDynamicalSystems)
                    (OSIInteractions)
                    (sizeMem)
                    (simulationLink))
SICONOS_IO_REGISTER_WITH_BASES(MLCP, (LinearOSNS),
                               (_n)
                               (_m)
                               (_curBlock)
                               //  (_numerics_problem)
                              )
SICONOS_IO_REGISTER_WITH_BASES(FrictionContact, (LinearOSNS),
                               (_contactProblemDim)
                               (_mu)
                               (_frictionContact_driver)
                               //  (_numerics_problem)
                              )
SICONOS_IO_REGISTER_WITH_BASES(NewtonEulerRImpact, (NewtonEulerR),
                               (_isOnContact)
                               (_Pc1)
                               (_Pc2)
                               (_Nc)
                               (_Mabs_C)
                               (_NPG1)
                               (_NPG2)
                               (_AUX1)
                               (_AUX2))
SICONOS_IO_REGISTER_WITH_BASES(QP, (OneStepNSProblem),
                               (Q)
                               (_p))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianR, (Relation),
                               (_jachq)
                               (_jachqDot))
SICONOS_IO_REGISTER(BlockCSRMatrix,
                    (nr)
                    (nc)
                    //  (numericsMatSparse)
                    //  (MBlockCSR)     ===> double* CSR : boost serialization fail!
                    (diagSizes)
                    (rowPos)
                    (colPos)
                   )
SICONOS_IO_REGISTER_WITH_BASES(LagrangianLinearTIR, (LagrangianR),
                               (_F)
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(NewtonImpactFrictionNSL, (NonSmoothLaw),
                               (_en)
                               (_et)
                               (_mu))
SICONOS_IO_REGISTER_WITH_BASES(LCP, (LinearOSNS),
                               //  (_numerics_problem)
                              )
SICONOS_IO_REGISTER_WITH_BASES(NewtonEulerR, (Relation),
                               (_ysize)
                               (_xsize)
                               (_qsize)
                               (_workQ)
                               (_jachq)
                               (_jachqDot)
                               (_jachlambda)
                               (_e)
                               (_yProj)
                               (_contactForce)
                               (_jachqT))
SICONOS_IO_REGISTER_WITH_BASES(EventDriven, (Simulation),
                               (istate)
                               (TOL_ED))
SICONOS_IO_REGISTER_WITH_BASES(TimeStepping, (Simulation),
                               (_computeResiduY)
                               (_computeResiduR)
                               (_newtonTolerance)
                               (_newtonMaxIteration)
                               (_newtonNbSteps)
                               (_newtonResiduDSMax)
                               (_newtonResiduYMax)
                               (_newtonResiduRMax)
                               (_newtonOptions))
SICONOS_IO_REGISTER(Interaction,
                    (_initialized)
                    (_id)
                    (_number)
                    (_relativeDegree)
                    (_interactionSize)
                    (_numberOfRelations)
                    (_sizeOfDS)
                    (_sizeZ)
                    (_y)
                    (_yOld)
                    (_y_k)
                    (_lambda)
                    (_lambdaOld)
                    (_involvedDS)
                    (_nslaw)
                    (_relation))
SICONOS_IO_REGISTER(DynamicalSystem,
                    (_number)
                    (_n)
                    (_x0)
                    (_residuFree)
                    (_r)
                    (_normRef)
                    (_x)
                    (_jacxRhs)
                    (_jacgx)
                    (_jacxDotG)
                    (_z)
                    (_g)
                    (_pluging)
                    (_pluginJacgx)
                    (_pluginJacxDotG)
                    (_xMemory)
                    (_stepsInMemory)
                    (_workV)
                    (_workMatrix)
                    (_workFree)
                    (count))
SICONOS_IO_REGISTER(Sensor,
                    (_type)
                    (_id)
                    (_data)
                    (_model)
                    (_timeDiscretisation)
                    (_eSensor))
SICONOS_IO_REGISTER(Topology,
                    (_minRelativeDegree)
                    (_maxRelativeDegree)
                    (_allInteractions)
                    (_DSG)
                    (_URG)
                    (_isTopologyUpToDate)
                    (_hasChanged)
                    (_numberOfConstraints)
                    (_symmetric))
SICONOS_IO_REGISTER(Actuator,
                    (_type)
                    (_id)
                    (_allSensors)
                    (_allDS)
                    (_model)
                    (_timeDiscretisation)
                    (_eActuator))
SICONOS_IO_REGISTER(SiconosException,
                    (reportMsg))
SICONOS_IO_REGISTER(SiconosMemory,
                    (maxSize)
                    (nbVectorsInMemory)
                    (_vectorMemory))
SICONOS_IO_REGISTER(Simulation,
                    (_name)
                    (_timeDiscretisation)
                    (_eventsManager)
                    (_tinit)
                    (_tend)
                    (_tout)
                    (_allOSI)
                    (_osiMap)
                    (_allNSProblems)
                    (_model)
                    (_levelMin)
                    (_levelMax)
                    (_tolerance)
                    (_printStat)
                    //  (statOut)
                    (_useRelativeConvergenceCriterion)
                    (_relativeConvergenceCriterionHeld)
                    (_relativeConvergenceTol))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianLinearTIDS, (LagrangianDS),
                               (_K)
                               (_C))
SICONOS_IO_REGISTER_WITH_BASES(GenericMechanical, (LinearOSNS),
                               (_pnumerics_GMP))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianScleronomousR, (LagrangianR),
                               (_pluginjqh)
                               (_pluginjqhdot)
                               (_NLh2dot))
SICONOS_IO_REGISTER(Model,
                    (_t)
                    (_t0)
                    (_T)
                    (_strat)
                    (_nsds)
                    (_title)
                    (_author)
                    (_description)
                    (_date))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderNonLinearDS, (DynamicalSystem),
                               (_M)
                               (_f)
                               (_fold)
                               (_jacobianfx)
                               (_pluginf)
                               (_pluginJacxf)
                               (_pluginM)
                               (_rMemory)
                               (_residur)
                               (_g_alpha)
                               (_xp)
                               (_xq)
                               (_invM))
SICONOS_IO_REGISTER_WITH_BASES(Relay, (LinearOSNS),
                               (_lb)
                               (_ub)
                               //  (_numerics_problem)
                              )
SICONOS_IO_REGISTER(OneStepNSProblem,
                    (_numerics_solver_id)
                    (_numerics_solver_options)
                    (_id)
                    (_sizeOutput)
                    (_DSBlocks)
                    (_simulation)
                    (_OSNSInteractions)
                    (_levelMin)
                    (_levelMax)
                    (_maxSize)
                    (_CPUtime)
                    (_nbIter)
                    (_numerics_options)
                    (_hasBeUpdated))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderLinearDS, (FirstOrderNonLinearDS),
                               (_A)
                               (_b)
                               (_pluginA)
                               (_pluginb))
SICONOS_IO_REGISTER_WITH_BASES(MLCP2, (MLCP),
                               (mFirstCall))
SICONOS_IO_REGISTER_WITH_BASES(LinearOSNS, (OneStepNSProblem),
                               (_w)
                               (_z)
                               (_M)
                               (_q)
                               (_MStorageType)
                               (_keepLambdaAndYState))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderType2R, (FirstOrderR),
                               (jacgx))
SICONOS_IO_REGISTER_WITH_BASES(TimeSteppingProjectOnConstraints, (TimeStepping),
                               (_constraintTol)
                               (_constraintTolUnilateral)
                               (_doProj)
                               (_doOnlyProj))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianRheonomousR, (LagrangianR),
                               (_hDot)
                               (_pluginhDot)
                               (_pluginJachq))
SICONOS_IO_REGISTER_WITH_BASES(MultipleImpactNSL, (NonSmoothLaw),
                               (_ResCof)
                               (_Stiff)
                               (_ElasCof))
SICONOS_IO_REGISTER_WITH_BASES(LagrangianCompliantR, (LagrangianR),
                               (_pluginJachq)
                               (_pluginJachlambda))
SICONOS_IO_REGISTER(TimeDiscretisation,
                    (h)
                    (k)
                    (tk)
                    (tdCase)
                    (pos))
SICONOS_IO_REGISTER(NonSmoothLaw,
                    (_size)
                    (_sizeProjectOnConstraints))
SICONOS_IO_REGISTER(BoundaryCondition,
                    (_velocityIndices)
                    (_prescribedVelocity)
                    (_prescribedVelocityOld)
                    (_pluginPrescribedVelocity))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderLinearR, (FirstOrderR),
                               (_F)
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderLinearTIR, (FirstOrderR),
                               (_F)
                               (_e))
SICONOS_IO_REGISTER_WITH_BASES(Equality, (LinearOSNS),
                               //  (_numerics_problem)
                              )
SICONOS_IO_REGISTER_WITH_BASES(FirstOrderR, (Relation),
                               (Jachx)
                               (Jacglambda))
SICONOS_IO_REGISTER_WITH_BASES(Moreau, (OneStepIntegrator),
                               (WMap)
                               (_WBoundaryConditionsMap)
                               (_theta)
                               (_gamma)
                               (_useGamma)
                               (_useGammaForRelation))
SICONOS_IO_REGISTER(Event,
                    //  (timeOfEvent)  => need __mpz_struct (gmp) serialization
                    (type)
                    (dTime)
                    (tick))
SICONOS_IO_REGISTER_WITH_BASES(ActuatorEvent, (Event),
                               (_actuator))
SICONOS_IO_REGISTER_WITH_BASES(OSNSMultipleImpact, (LinearOSNS),
                               (Impulse_variable)
                               (Time_variable)
                               (Ncontact)
                               (NstepEst)
                               (NstepMax)
                               (TOL_IMPACT)
                               (TypeCompLaw)
                               (VelContact)
                               (OldVelContact)
                               (EnerContact)
                               (WcContact)
                               (DistriVector)
                               (StateContact)
                               (Kcontact)
                               (ResContact)
                               (ElasCoefContact)
                               (DelImpulseContact)
                               (TolImpulseContact)
                               (ImpulseContact_update)
                               (ForceContact)
                               (SelectPrimaConInVel)
                               (IdPrimaContact)
                               (IsPrimaConEnergy)
                               (VelAtPrimaCon)
                               (EnerAtPrimaCon)
                               (DeltaP)
                               (OutputFile)
                               (NameFile)
                               (YesSaveData)
                               (NstepSave)
                               (IsNumberOfStepsEst)
                               (_DataMatrix)
                               (YesSaveByMatrix)
                               (SizeDataSave)
                               (_IsImpactEnd))
SICONOS_IO_REGISTER(PluggedObject,
                    //  (fPtr)
                    (pluginName))
SICONOS_IO_REGISTER(EventsManager,
                    //  (_allEvents)        => segfault (save)
                    //  (_currentEvent)
                    //  (_nextEvent)
                    //  (_ETD)
                    //  (_ENonSmooth)
                    //  (_simulation)
                    //  (_hasNS)
                    //  (_hasCM)
                    // (GapLimit2Events)
                   )
SICONOS_IO_REGISTER_WITH_BASES(LagrangianDS, (DynamicalSystem),
                               (_ndof)
                               (_q)
                               (_q0)
                               (_velocity0)
                               (_qMemory)
                               (_velocityMemory)
                               (_p)
                               (_mass)
                               (_fInt)
                               (_jacobianFIntq)
                               (_jacobianFIntqDot)
                               (_fExt)
                               (_NNL)
                               (_jacobianNNLq)
                               (_jacobianNNLqDot)
                               (_fL)
                               (_jacobianqFL)
                               (_jacobianqDotFL)
                               (_boundaryConditions)
                               (_reactionToBoundaryConditions)
                               (_pluginMass)
                               (_pluginFInt)
                               (_pluginFExt)
                               (_pluginNNL)
                               (_pluginJacqFInt)
                               (_pluginJacqDotFInt)
                               (_pluginJacqNNL)
                               (_pluginJacqDotNNL))
SICONOS_IO_REGISTER(UnitaryRelationsSet,)
//  (fpt)
//  (setOfT))
SICONOS_IO_REGISTER(UnitaryRelationsGraph,
                    //  (vertex_descriptor)
                    (g))
SICONOS_IO_REGISTER(InteractionsSet,)
//  (fpt)
//  (setOfT))
SICONOS_IO_REGISTER(DynamicalSystemsGraph,
                    //  (vertex_descriptor)
                    (g))

/* hand written */

SICONOS_IO_REGISTER(NumericsOptions, (verboseMode));

SICONOS_IO_REGISTER(RelationData, (block)(blockProj)(source)(target));

SICONOS_IO_REGISTER(SystemData, (upper_block)(lower_block)(upper_blockProj)(lower_blockProj));

BOOST_TYPEOF_REGISTER_TYPE(_SolverOptions);

BOOST_TYPEOF_REGISTER_TYPE(LinearComplementarityProblem);

#include <boost/preprocessor/tuple/elem.hpp>

#define SERIALIZE_I(r,T,M)                                              \
  BOOST_PP_TUPLE_ELEM(2,0,T) & ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(M), BOOST_PP_TUPLE_ELEM(2,1,T) . M);

#define SERIALIZE(S,MEMBERS, ARCHIVE)                       \
  BOOST_PP_SEQ_FOR_EACH(SERIALIZE_I, (ARCHIVE, S), MEMBERS)

#define SERIALIZE_C_ARRAY(DIM, STRUCT, ARRAY, ARCHIVE)                   \
  if (Archive::is_loading::value)                                       \
  {                                                                     \
    STRUCT . ARRAY = (BOOST_TYPEOF(STRUCT . ARRAY)) malloc(DIM * sizeof(BOOST_TYPEOF(* (STRUCT . ARRAY)))); \
  };                                                                    \
  {                                                                     \
    boost::serialization::array<BOOST_TYPEOF(*(STRUCT . ARRAY))>        \
      wrapper = boost::serialization::make_array(STRUCT . ARRAY,DIM);   \
    ARCHIVE & boost::serialization::make_nvp(BOOST_PP_STRINGIZE(ARRAY),wrapper); \
  }                                                                     \
 
template <class Archive>
void siconos_io(Archive& ar, _SolverOptions&v, unsigned int version)
{
  SERIALIZE(v, (solverId)(isSet)(iSize)(dSize)(filterOn)(numberOfInternalSolvers), ar);

  SERIALIZE_C_ARRAY(v.iSize, v, iparam, ar);
  SERIALIZE_C_ARRAY(v.dSize, v, dparam, ar);
  SERIALIZE_C_ARRAY(v.numberOfInternalSolvers, v, internalSolvers, ar);
}

template <class Archive>
void siconos_io(Archive& ar, LinearComplementarityProblem& v, unsigned int version)
{
  SERIALIZE(v, (size)(M), ar);
  SERIALIZE_C_ARRAY(v.size, v, q, ar);
}

template <class Archive>
void siconos_io(Archive& ar, SparseBlockStructuredMatrix& v, unsigned int version)
{
  SERIALIZE(v, (nbblocks)(blocknumber0)(blocknumber1)(filled1)(filled2), ar);
  SERIALIZE_C_ARRAY(v.filled1, v, index1_data, ar);
  SERIALIZE_C_ARRAY(v.filled2, v, index2_data, ar);
}

template <class Archive>
void siconos_io(Archive&ar, NumericsMatrix& v, unsigned int version)
{
  SERIALIZE(v, (storageType)(size0)(size1)(matrix1), ar);
  SERIALIZE_C_ARRAY(v.size0 * v.size1, v, matrix0, ar);
}

template <class Archive>
void siconos_io(Archive& ar, DynamicalSystemsSet& v, unsigned int version)
{
  ar &  boost::serialization::make_nvp("ThisShouldNotBeASetAnyMore",
                                       boost::serialization::base_object< std::vector<SP::DynamicalSystem> >(v));
}

template <class Archive>
void siconos_io(Archive & ar, SimpleVector & v, unsigned int version)
{
  ar & boost::serialization::make_nvp("_dense", v._dense);
  if (v._dense)
  {
    ar & boost::serialization::make_nvp("vect", v.vect.Dense);
  }
  if (!v._dense)
  {
    ar & boost::serialization::make_nvp("vect", v.vect.Sparse);
  }
  ar &  boost::serialization::make_nvp("SiconosVector",
                                       boost::serialization::base_object<SiconosVector>(v));
}

template <class Archive>
void siconos_io(Archive & ar, SiconosVector & v, unsigned int version)
{

}


template <class Archive>
void siconos_io(Archive & ar, SimpleMatrix & m, unsigned int version)
{
  ar & boost::serialization::make_nvp("num", m.num);
  ar & boost::serialization::make_nvp("dimRow", m.dimRow);
  ar & boost::serialization::make_nvp("dimCol", m.dimCol);
  ar & boost::serialization::make_nvp("ipiv", m.ipiv);
  ar & boost::serialization::make_nvp("isPLUFactorized", m.isPLUFactorized);
  ar & boost::serialization::make_nvp("isPLUInversed", m.isPLUInversed);
  switch (m.num)
  {
  case 1:
  {
    ar & boost::serialization::make_nvp("mat", m.mat.Dense);
    break;
  }
  case 2:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Triang);
    break;
  }
  case 3:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Sym);
    break;
  }
  case 4:
  {
    ar & boost::serialization::make_nvp("mat", m.mat.Sparse);
    break;
  }
  case 5:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Banded);
    break;
  }
  case 6:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Zero);
    break;
  }
  case 7:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Identity);
    break;
  }
  }
  ar &  boost::serialization::make_nvp("SiconosMatrix",
                                       boost::serialization::base_object<SiconosMatrix>(m));
}


namespace boost
{
namespace serialization
{

template <class Archive>
void serialize(Archive& ar, NumericsMatrix& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, SparseBlockStructuredMatrix& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, LinearComplementarityProblem& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, _SolverOptions& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, DynamicalSystemsSet& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, SimpleVector& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, SimpleMatrix& m, unsigned int version)
{
  siconos_io(ar, m, version);
}

template <class Archive>
void serialize(Archive& ar, SiconosVector& v, unsigned int version)
{
  siconos_io(ar, v, version);
}
}
}

void KernelTest::setUp() {};
void KernelTest::tearDown() {};

void KernelTest::t0()
{
  SP::SiconosVector q(new SimpleVector(3));
  SP::SiconosVector q0(new SimpleVector(3));

  (*q)(0) = 1.0;
  (*q)(1) = 2.0;
  (*q)(2) = 2.0;


  std::ofstream ofs("Kernelt0.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<SimpleVector*>(NULL));
    oa << NVP(q);
  }

  std::ifstream ifs("Kernelt0.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<SimpleVector*>(NULL));
    ia >> NVP(q0);
  }

  CPPUNIT_ASSERT(*q0 == *q);
}


void KernelTest::t1()
{
  SP::SiconosMatrix m1(new SimpleMatrix(3, 3));
  SP::SimpleMatrix m2(new SimpleMatrix(3, 3));

  m1->eye();
  (*m1)(1, 0) = 3.0;
  (*m1)(2, 1) = -7;

  std::ofstream ofs("Kernelt1.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<SimpleMatrix*>(NULL));
    oa << NVP(m1);
  }

  std::ifstream ifs("Kernelt1.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<SimpleMatrix*>(NULL));
    ia >> NVP(m2);
  }

  m1->display();
  m2->display();

  CPPUNIT_ASSERT(*m1 == *m2);

}

void KernelTest::t2()
{
  SP::SiconosMatrix m(new SimpleMatrix(3, 3));
  SP::SiconosVector v(new SimpleVector(3));
  SP::SiconosVector q(new SimpleVector(3));

  m->eye();


  SP::DynamicalSystem ds1(new LagrangianDS(q, v, m));
  SP::DynamicalSystem ds2(new LagrangianDS(q, v, m));

  std::ofstream ofs("Kernelt2.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<SimpleMatrix*>(NULL));
    oa.register_type(static_cast<SimpleVector*>(NULL));
    oa.register_type(static_cast<LagrangianDS*>(NULL));
    oa << NVP(ds1);
  }

  std::ifstream ifs("Kernelt2.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<SimpleMatrix*>(NULL));
    ia.register_type(static_cast<SimpleVector*>(NULL));
    ia.register_type(static_cast<LagrangianDS*>(NULL));
    ia >> NVP(ds2);
  }

  CPPUNIT_ASSERT(*(boost::static_pointer_cast<LagrangianDS>(ds1)->mass())
                 == *(boost::static_pointer_cast<LagrangianDS>(ds2)->mass()));
  CPPUNIT_ASSERT(*(boost::static_pointer_cast<LagrangianDS>(ds1)->q())
                 == *(boost::static_pointer_cast<LagrangianDS>(ds2)->q()));
  CPPUNIT_ASSERT(*(boost::static_pointer_cast<LagrangianDS>(ds1)->velocity())
                 == *(boost::static_pointer_cast<LagrangianDS>(ds2)->velocity()));

}


void KernelTest::t3()
{
  SP::SolverOptions so(new SolverOptions);
  SP::SolverOptions sor(new SolverOptions);
  so->solverId = SICONOS_FRICTION_3D_GLOBALAC;
  so->isSet = 36;
  so->iSize = 10;
  so->iparam = (int*) malloc(sizeof(int) * so->iSize);
  so->dSize = 10;
  so->dparam = (double *)malloc(sizeof(double) * so->dSize);
  so->filterOn = 1;
  so->numberOfInternalSolvers = 1;
  so->internalSolvers = (_SolverOptions *) malloc(sizeof(_SolverOptions) * so->numberOfInternalSolvers);


  std::ofstream ofs("SolverOptions.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa << NVP(so);
  }

  std::ifstream ifs("SolverOptions.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia >> NVP(sor);
  }

  CPPUNIT_ASSERT((so->iSize == sor->iSize));

}

void KernelTest::t4()
{
  SP::SiconosMatrix m(new SimpleMatrix(3, 3));
  SP::SiconosVector v(new SimpleVector(3));
  SP::SiconosVector q(new SimpleVector(3));

  m->eye();


  SP::DynamicalSystem ds1(new LagrangianDS(q, v, m));
  SP::DynamicalSystem ds2(new LagrangianDS(q, v, m));

  SP::DynamicalSystemsSet dsset(new DynamicalSystemsSet());

  dsset->insert(ds1);
  dsset->insert(ds2);

  std::ofstream ofs("t4.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<SimpleMatrix*>(NULL));
    oa.register_type(static_cast<SimpleVector*>(NULL));
    oa.register_type(static_cast<LagrangianDS*>(NULL));
    oa << NVP(dsset);
  }

  SP::DynamicalSystemsSet dssetfromfile(new DynamicalSystemsSet());

  std::ifstream ifs("t4.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<SimpleMatrix*>(NULL));
    ia.register_type(static_cast<SimpleVector*>(NULL));
    ia.register_type(static_cast<LagrangianDS*>(NULL));
    ia >> NVP(dssetfromfile);
  }

}


using namespace std;
void KernelTest::t5()
{

  // ================= Creation of the model =======================

  // User-defined main parameters
  unsigned int nDof = 3;           // degrees of freedom for the ball
  double t0 = 0;                   // initial computation time
  double T = 10;                  // final computation time
  double h = 0.005;                // time step
  double position_init = 1.0;      // initial position for lowest bead.
  double velocity_init = 0.0;      // initial velocity for lowest bead.
  double theta = 0.5;              // theta for Moreau integrator
  double R = 0.1; // Ball radius
  double m = 1; // Ball mass
  double g = 9.81; // Gravity
  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  cout << "====> Model loading ..." << endl << endl;

  SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
  (*Mass)(0, 0) = m;
  (*Mass)(1, 1) = m;
  (*Mass)(2, 2) = 3. / 5 * m * R * R;

  // -- Initial positions and velocities --
  SP::SimpleVector q0(new SimpleVector(nDof));
  SP::SimpleVector v0(new SimpleVector(nDof));
  (*q0)(0) = position_init;
  (*v0)(0) = velocity_init;

  // -- The dynamical system --
  SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));

  // -- Set external forces (weight) --
  SP::SimpleVector weight(new SimpleVector(nDof));
  (*weight)(0) = -m * g;
  ball->setFExtPtr(weight);

  // --------------------
  // --- Interactions ---
  // --------------------

  // -- nslaw --
  double e = 0.9;

  // Interaction ball-floor
  //
  SP::SiconosMatrix H(new SimpleMatrix(1, nDof));
  (*H)(0, 0) = 1.0;

  SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
  SP::Relation relation(new LagrangianLinearTIR(H));

  SP::Interaction inter(new Interaction(1, nslaw, relation));

  // -------------
  // --- Model ---
  // -------------
  SP::Model bouncingBall(new Model(t0, T));

  // add the dynamical system in the non smooth dynamical system
  bouncingBall->nonSmoothDynamicalSystem()->insertDynamicalSystem(ball);

  // link the interaction and the dynamical system
  bouncingBall->nonSmoothDynamicalSystem()->link(inter, ball);


  // ------------------
  // --- Simulation ---
  // ------------------

  // -- (1) OneStepIntegrators --
  SP::Moreau OSI(new Moreau(ball, theta));

  // -- (2) Time discretisation --
  SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

  // -- (3) one step non smooth problem
  SP::OneStepNSProblem osnspb(new LCP());

  // -- (4) Simulation setup with (1) (2) (3)
  SP::TimeStepping s(new TimeStepping(t, OSI, osnspb));

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Initialisation ..." << endl << endl;
  bouncingBall->initialize(s);


  std::ofstream ofs("BouncingBall1.xml");
  {
    boost::archive::xml_oarchive oa(ofs);

    SP::SolverOptions so;


    oa.register_type(static_cast<BlockVector*>(NULL));
    //    oa.register_type(static_cast<SensorPosition*>(NULL));
    oa.register_type(static_cast<NewtonImpactNSL*>(NULL));
    //    oa.register_type(static_cast<NewtonEulerDS*>(NULL));
    //    oa.register_type(static_cast<PrimalFrictionContact*>(NULL));
    //   oa.register_type(static_cast<RelayNSL*>(NULL));
    //   oa.register_type(static_cast<MixedComplementarityConditionNSL*>(NULL));
    //    oa.register_type(static_cast<SensorEvent*>(NULL));
    //    oa.register_type(static_cast<MLCP*>(NULL));

    //    oa.register_type(static_cast<NewtonEulerRImpact*>(NULL));
    //    oa.register_type(static_cast<QP*>(NULL));
    //    oa.register_type(static_cast<LagrangianR*>(NULL));
    oa.register_type(static_cast<LagrangianLinearTIR*>(NULL));
    oa.register_type(static_cast<SimpleVector*>(NULL));
    //    oa.register_type(static_cast<NewtonImpactFrictionNSL*>(NULL));
    //    oa.register_type(static_cast<NewtonEulerR*>(NULL));
    //    oa.register_type(static_cast<EventDriven*>(NULL));
    oa.register_type(static_cast<TimeStepping*>(NULL));
    oa.register_type(static_cast<LagrangianLinearTIDS*>(NULL));
    //    oa.register_type(static_cast<GenericMechanical*>(NULL));
    oa.register_type(static_cast<LagrangianScleronomousR*>(NULL));
    //    oa.register_type(static_cast<FirstOrderNonLinearDS*>(NULL));
    //    oa.register_type(static_cast<Lsodar*>(NULL));
    //    oa.register_type(static_cast<Relay*>(NULL));
    //    oa.register_type(static_cast<FirstOrderLinearDS*>(NULL));
    //    oa.register_type(static_cast<MLCP2*>(NULL));
    //    oa.register_type(static_cast<OneStepNSProblem*>(NULL));
    oa.register_type(static_cast<LCP*>(NULL));
    //    oa.register_type(static_cast<LinearOSNS*>(NULL));
    //    oa.register_type(static_cast<FirstOrderType2R*>(NULL));
    //   oa.register_type(static_cast<TimeSteppingProjectOnConstraints*>(NULL));
    oa.register_type(static_cast<LagrangianRheonomousR*>(NULL));
    //    oa.register_type(static_cast<MultipleImpactNSL*>(NULL));
    //    oa.register_type(static_cast<LagrangianCompliantR*>(NULL));
    //    oa.register_type(static_cast<FirstOrderLinearR*>(NULL));
    oa.register_type(static_cast<SimpleMatrix*>(NULL));
    oa.register_type(static_cast<BlockMatrix*>(NULL));
    //    oa.register_type(static_cast<FirstOrderLinearTIR*>(NULL));
    //    oa.register_type(static_cast<Equality*>(NULL));
    //    oa.register_type(static_cast<FirstOrderR*>(NULL));
    oa.register_type(static_cast<Moreau*>(NULL));
    //    oa.register_type(static_cast<ActuatorEvent*>(NULL));
    //    oa.register_type(static_cast<Event*>(NULL));
    //    oa.register_type(static_cast<OSNSMultipleImpact*>(NULL));
    oa.register_type(static_cast<LagrangianDS*>(NULL));

    // dump
    DEBUG_PRINT("saving\n");
    oa << NVP(bouncingBall);
  }



  SP::Model bouncingBallFromFile(new Model());

  std::ifstream ifs("BouncingBall1.xml");
  {
    boost::archive::xml_iarchive ia(ifs);

    ia.register_type(static_cast<BlockVector*>(NULL));
    //    ia.register_type(static_cast<SensorPosition*>(NULL));
    ia.register_type(static_cast<NewtonImpactNSL*>(NULL));
    //    ia.register_type(static_cast<NewtonEulerDS*>(NULL));
    //    ia.register_type(static_cast<PrimalFrictionContact*>(NULL));
    //   ia.register_type(static_cast<RelayNSL*>(NULL));
    //   ia.register_type(static_cast<MixedComplementarityConditionNSL*>(NULL));
    //    ia.register_type(static_cast<SensorEvent*>(NULL));
    //    ia.register_type(static_cast<MLCP*>(NULL));

    //    ia.register_type(static_cast<NewtonEulerRImpact*>(NULL));
    //    ia.register_type(static_cast<QP*>(NULL));
    //    ia.register_type(static_cast<LagrangianR*>(NULL));
    ia.register_type(static_cast<LagrangianLinearTIR*>(NULL));
    ia.register_type(static_cast<SimpleVector*>(NULL));
    //    ia.register_type(static_cast<NewtonImpactFrictionNSL*>(NULL));
    //    ia.register_type(static_cast<NewtonEulerR*>(NULL));
    //    ia.register_type(static_cast<EventDriven*>(NULL));
    ia.register_type(static_cast<TimeStepping*>(NULL));
    ia.register_type(static_cast<LagrangianLinearTIDS*>(NULL));
    //    ia.register_type(static_cast<GenericMechanical*>(NULL));
    ia.register_type(static_cast<LagrangianScleronomousR*>(NULL));
    //    ia.register_type(static_cast<FirstOrderNonLinearDS*>(NULL));
    //    ia.register_type(static_cast<Lsodar*>(NULL));
    //    ia.register_type(static_cast<Relay*>(NULL));
    //    ia.register_type(static_cast<FirstOrderLinearDS*>(NULL));
    //    ia.register_type(static_cast<MLCP2*>(NULL));
    //    ia.register_type(static_cast<OneStepNSProblem*>(NULL));
    ia.register_type(static_cast<LCP*>(NULL));
    //    ia.register_type(static_cast<LinearOSNS*>(NULL));
    //    ia.register_type(static_cast<FirstOrderType2R*>(NULL));
    //   ia.register_type(static_cast<TimeSteppingProjectOnConstraints*>(NULL));
    ia.register_type(static_cast<LagrangianRheonomousR*>(NULL));
    //    ia.register_type(static_cast<MultipleImpactNSL*>(NULL));
    //    ia.register_type(static_cast<LagrangianCompliantR*>(NULL));
    //    ia.register_type(static_cast<FirstOrderLinearR*>(NULL));
    ia.register_type(static_cast<SimpleMatrix*>(NULL));
    ia.register_type(static_cast<BlockMatrix*>(NULL));
    //    ia.register_type(static_cast<FirstOrderLinearTIR*>(NULL));
    //    ia.register_type(static_cast<Equality*>(NULL));
    //    ia.register_type(static_cast<FirstOrderR*>(NULL));
    ia.register_type(static_cast<Moreau*>(NULL));
    //    ia.register_type(static_cast<ActuatorEvent*>(NULL));
    //    ia.register_type(static_cast<Event*>(NULL));
    //    ia.register_type(static_cast<OSNSMultipleImpact*>(NULL));
    ia.register_type(static_cast<LagrangianDS*>(NULL));

    DEBUG_PRINT("loading\n");
    ia >> NVP(bouncingBallFromFile);
  }


  CPPUNIT_ASSERT((bouncingBallFromFile->t0() == bouncingBall->t0()));
  // in depth comparison?

  // now we should try to run the bouncing ball from file

  // BUT: non serialized members => must be initialized or serialized

}
