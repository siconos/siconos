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

#include "../../tools/Full.hpp"

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
void siconos_io(Archive& ar, InteractionsSet& v, unsigned int)
{
  if (Archive::is_loading::value)
  {
    v.fpt = &Interaction::getSort;
  }
  ar & boost::serialization::make_nvp("setOfT", v.setOfT);
}



template <class Archive>
void siconos_io(Archive& ar, DynamicalSystemsGraph& v, unsigned int version)
{

  ar & boost::serialization::make_nvp("g", v.g);

  if (Archive::is_loading::value)
  {
    DynamicalSystemsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = v.vertices(); vi != viend; ++vi)
    {
      v.vertex_descriptor[v.bundle(*vi)] = *vi;
    }
  }

}

template <class Archive>
void siconos_io(Archive& ar, UnitaryRelationsGraph& v, unsigned int version)
{

  ar & boost::serialization::make_nvp("g", v.g);

  if (Archive::is_loading::value)
  {
    DynamicalSystemsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = v.vertices(); vi != viend; ++vi)
    {
      v.vertex_descriptor[v.bundle(*vi)] = *vi;
    }
  }

}





template <class Archive>
void siconos_io(Archive& ar, std::basic_ofstream<char>&v , unsigned int version)
{
  // do nothing
}


template <class Archive>
void siconos_io(Archive& ar, __mpz_struct& v, unsigned int version)
{
  SERIALIZE(v, (_mp_alloc)(_mp_size), ar);
  SERIALIZE_C_ARRAY(v._mp_alloc, v, _mp_d, ar);
}


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
void serialize(Archive& ar, __mpz_struct& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, InteractionsSet& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, PluginHandle& v, unsigned int version)
{

}

template <class Archive>
void serialize(Archive& ar, UnitaryRelationsGraph& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, DynamicalSystemsGraph& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, std::basic_ofstream<char>& v, unsigned int version)
{
  siconos_io(ar, v, version);
}



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
    oa.register_type(static_cast<SensorEvent*>(NULL));
    oa.register_type(static_cast<ActuatorEvent*>(NULL));
    oa.register_type(static_cast<NonSmoothEvent*>(NULL));
    oa.register_type(static_cast<TimeDiscretisationEvent*>(NULL));
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
    ia.register_type(static_cast<SensorEvent*>(NULL));
    ia.register_type(static_cast<ActuatorEvent*>(NULL));
    ia.register_type(static_cast<NonSmoothEvent*>(NULL));
    ia.register_type(static_cast<TimeDiscretisationEvent*>(NULL));
    ia.register_type(static_cast<LagrangianDS*>(NULL));

    DEBUG_PRINT("loading\n");
    ia >> NVP(bouncingBallFromFile);
  }


  CPPUNIT_ASSERT((bouncingBallFromFile->t0() == bouncingBall->t0()));
  // in depth comparison?

  // now we should try to run the bouncing ball from file

  // BUT: non serialized members => must be initialized or serialized

}


void KernelTest::t6()
{
  SP::Model bouncingBall(new Model());

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
    ia.register_type(static_cast<SensorEvent*>(NULL));
    ia.register_type(static_cast<ActuatorEvent*>(NULL));
    ia.register_type(static_cast<NonSmoothEvent*>(NULL));
    ia.register_type(static_cast<TimeDiscretisationEvent*>(NULL));
    ia.register_type(static_cast<LagrangianDS*>(NULL));

    DEBUG_PRINT("loading\n");
    ia >> NVP(bouncingBall);

    try
    {
      double T = bouncingBall->finalT();
      double t0 = bouncingBall->t0();
      double h = bouncingBall->simulation()->timeStep();
      int N = (int)((T - t0) / h); // Number of time steps



      SP::LagrangianDS ball = boost::static_pointer_cast<LagrangianDS>
                              (bouncingBall->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0));

      SP::Interaction inter = *(bouncingBall->nonSmoothDynamicalSystem()->interactions()->begin());
      SP::TimeStepping s = boost::static_pointer_cast<TimeStepping>(bouncingBall->simulation());


      // --- Get the values to be plotted ---
      // -> saved in a matrix dataPlot
      unsigned int outputSize = 5;
      SimpleMatrix dataPlot(N + 1, outputSize);



      SP::SiconosVector q = ball->q();
      SP::SiconosVector v = ball->velocity();
      SP::SiconosVector p = ball->p(1);
      SP::SiconosVector lambda = inter->lambda(1);

      dataPlot(0, 0) = bouncingBall->t0();
      dataPlot(0, 1) = (*q)(0);
      dataPlot(0, 2) = (*v)(0);
      dataPlot(0, 3) = (*p)(0);
      dataPlot(0, 4) = (*lambda)(0);
      // --- Time loop ---
      cout << "====> Start computation ... " << endl << endl;
      // ==== Simulation loop - Writing without explicit event handling =====
      int k = 1;
      boost::progress_display show_progress(N);

      boost::timer time;
      time.restart();

      while (s->nextTime() < T)
      {
        s->computeOneStep();

        // --- Get values to be plotted ---
        dataPlot(k, 0) =  s->nextTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*lambda)(0);
        s->nextStep();
        ++show_progress;
        k++;
      }
      cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
      cout << "Computation Time " << time.elapsed()  << endl;

      // --- Output files ---
      cout << "====> Output file writing ..." << endl;
      ioMatrix io("result.dat", "ascii");
      dataPlot.resize(k, outputSize);
      io.write(dataPlot, "noDim");
      // Comparison with a reference file
      SimpleMatrix dataPlotRef(dataPlot);
      dataPlotRef.zero();
      ioMatrix ref("result.ref", "ascii");
      ref.read(dataPlotRef);

      if ((dataPlot - dataPlotRef).normInf() > 1e-12)
      {
        std::cout << "Warning. The results is rather different from the reference file." << std::endl;
        CPPUNIT_ASSERT(false);
      }

    }

    catch (SiconosException e)
    {
      cout << e.report() << endl;
      CPPUNIT_ASSERT(false);
    }
    catch (...)
    {
      cout << "Exception caught in BouncingBallTS.cpp" << endl;
      CPPUNIT_ASSERT(false);

    }


  }
}
