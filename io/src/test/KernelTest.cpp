#include "SiconosConfig.h"
#define WITH_SERIALIZATION

#ifdef HAVE_SICONOS_MECHANICS
#include "MechanicsIO.hpp"
#endif

#include "KernelTest.hpp"
#include "SiconosKernel.hpp"

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_sparse.hpp>
#include <boost/numeric/bindings/ublas/matrix_sparse.hpp>

#define DEBUG_MESSAGES 1
#include "../SiconosFull.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

CPPUNIT_TEST_SUITE_REGISTRATION(KernelTest);

void KernelTest::setUp()
{
  BBxml = "BouncingBall1.xml";
}

void KernelTest::tearDown() {};

void KernelTest::t0()
{
  SP::SiconosVector q(new SiconosVector(3));
  SP::SiconosVector q0(new SiconosVector(3));

  (*q)(0) = 1.0;
  (*q)(1) = 2.0;
  (*q)(2) = 2.0;


  std::ofstream ofs("Kernelt0.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<SiconosVector*>(NULL));
    oa << NVP(q);
  }

  std::ifstream ifs("Kernelt0.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<SiconosVector*>(NULL));
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
  SP::SiconosVector v(new SiconosVector(3));
  SP::SiconosVector q(new SiconosVector(3));

  m->eye();


  SP::DynamicalSystem ds1(new LagrangianDS(q, v, m));
  SP::DynamicalSystem ds2(new LagrangianDS(q, v, m));

  std::ofstream ofs("Kernelt2.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<SimpleMatrix*>(NULL));
    oa.register_type(static_cast<SiconosVector*>(NULL));
    oa.register_type(static_cast<LagrangianDS*>(NULL));
    oa << NVP(ds1);
  }

  std::ifstream ifs("Kernelt2.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<SimpleMatrix*>(NULL));
    ia.register_type(static_cast<SiconosVector*>(NULL));
    ia.register_type(static_cast<LagrangianDS*>(NULL));
    ia >> NVP(ds2);
  }

  CPPUNIT_ASSERT(*(std11::static_pointer_cast<LagrangianDS>(ds1)->mass())
                 == *(std11::static_pointer_cast<LagrangianDS>(ds2)->mass()));
  CPPUNIT_ASSERT(*(std11::static_pointer_cast<LagrangianDS>(ds1)->q())
                 == *(std11::static_pointer_cast<LagrangianDS>(ds2)->q()));
  CPPUNIT_ASSERT(*(std11::static_pointer_cast<LagrangianDS>(ds1)->velocity())
                 == *(std11::static_pointer_cast<LagrangianDS>(ds2)->velocity()));

}


void KernelTest::t3()
{
  SP::SolverOptions so(new SolverOptions);
  SP::SolverOptions sor(new SolverOptions);
  solver_options_fill(so.get(), SICONOS_FRICTION_3D_NSN_AC, 10, 10, 0, 0.);
  so->numberOfInternalSolvers = 1;
  so->internalSolvers = (SolverOptions *) malloc(sizeof(SolverOptions) * so->numberOfInternalSolvers);
  solver_options_fill(so->internalSolvers, SICONOS_FRICTION_3D_NSN_AC, 10, 10, 0, 0.);

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

  solver_options_delete(so.get());
  solver_options_delete(sor.get());

}

// void KernelTest::t4()
// {
//   SP::SiconosMatrix m(new SimpleMatrix(3, 3));
//   SP::SiconosVector v(new SiconosVector(3));
//   SP::SiconosVector q(new SiconosVector(3));

//   m->eye();


//   SP::DynamicalSystem ds1(new LagrangianDS(q, v, m));
//   SP::DynamicalSystem ds2(new LagrangianDS(q, v, m));

//   SP::DynamicalSystemsSet dsset(new DynamicalSystemsSet());

//   dsset->insert(ds1);
//   dsset->insert(ds2);

//   std::ofstream ofs("t4.xml");
//   {
//     boost::archive::xml_oarchive oa(ofs);
//     oa.register_type(static_cast<SimpleMatrix*>(NULL));
//     oa.register_type(static_cast<SiconosVector*>(NULL));
//     oa.register_type(static_cast<LagrangianDS*>(NULL));
//     oa << NVP(dsset);
//   }

//   SP::DynamicalSystemsSet dssetfromfile(new DynamicalSystemsSet());

//   std::ifstream ifs("t4.xml");
//   {
//     boost::archive::xml_iarchive ia(ifs);
//     ia.register_type(static_cast<SimpleMatrix*>(NULL));
//     ia.register_type(static_cast<SiconosVector*>(NULL));
//     ia.register_type(static_cast<LagrangianDS*>(NULL));
//     ia >> NVP(dssetfromfile);
//   }

// }


#include "SiconosRestart.hpp"

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
  double theta = 0.5;              // theta for MoreauJeanOSI integrator
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
  SP::SiconosVector q0(new SiconosVector(nDof));
  SP::SiconosVector v0(new SiconosVector(nDof));
  (*q0)(0) = position_init;
  (*v0)(0) = velocity_init;

  // -- The dynamical system --
  SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));

  // -- Set external forces (weight) --
  SP::SiconosVector weight(new SiconosVector(nDof));
  (*weight)(0) = -m * g;
  ball->setFExtPtr(weight);

  // --------------------
  // --- Interactions ---
  // --------------------

  // -- nslaw --
  double e = 0.9;

  // Interaction ball-floor
  //
  SP::SimpleMatrix H(new SimpleMatrix(1, nDof));
  (*H)(0, 0) = 1.0;

  SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
  SP::Relation relation(new LagrangianLinearTIR(H));

  SP::Interaction inter(new Interaction(nslaw, relation));

  // -------------
  // --- Model ---
  // -------------
  SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));

  // add the dynamical system in the non smooth dynamical system
  bouncingBall->insertDynamicalSystem(ball);

  // link the interaction and the dynamical system
  bouncingBall->link(inter, ball);


  // ------------------
  // --- Simulation ---
  // ------------------

  // -- (1) OneStepIntegrators --
  SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));

  // -- (2) Time discretisation --
  SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

  // -- (3) one step non smooth problem
  SP::OneStepNSProblem osnspb(new LCP());

  // -- (4) Simulation setup with (1) (2) (3)
  SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));
  s->associate(OSI, ball);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  Siconos::save(s, BBxml);

  SP::Simulation simFromFile = Siconos::load(BBxml);
  SP::NonSmoothDynamicalSystem bouncingBallFromFile =
    simFromFile->nonSmoothDynamicalSystem();

  CPPUNIT_ASSERT((bouncingBallFromFile->t0() == bouncingBall->t0()));
  // in depth comparison?

  // now we should try to run the bouncing ball from file

  // BUT: non serialized members => must be initialized or serialized

}


void KernelTest::t6()
{
  SP::Simulation s = Siconos::load(BBxml);

  try
  {
    SP::NonSmoothDynamicalSystem bouncingBall = s->nonSmoothDynamicalSystem();

    double T = bouncingBall->finalT();
    double t0 = bouncingBall->t0();
    double h = s->timeStep();
    int N = (int)((T - t0) / h); // Number of time steps

    SP::DynamicalSystemsGraph dsg =
      bouncingBall->topology()->dSG(0);

    SP::LagrangianDS ball = std11::static_pointer_cast<LagrangianDS>
      (dsg->bundle(*(dsg->begin())));

    SP::TimeStepping s = std11::static_pointer_cast<TimeStepping>(s);
    SP::Interaction inter;
    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet0 = bouncingBall->topology()->indexSet(0);
    for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
      inter = indexSet0->bundle(*ui);


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

    while (s->hasNextEvent())
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
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("result.ref", "ascii", dataPlotRef);

    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout <<
        "Warning. The results is rather different from the reference file :"
                <<
        (dataPlot - dataPlotRef).normInf()
                <<
        std::endl;
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

#ifdef HAVE_SICONOS_MECHANICS
#include <Disk.hpp>

void KernelTest::t7()
{

  SP::DynamicalSystem ds1, ds2;

  // Must be size=1, cannot deserialize a LagrangianDS with _ndof==0
  SP::SiconosVector q(new SiconosVector(1));
  SP::SiconosVector v(new SiconosVector(1));

  ds1.reset(new Disk(1, 1, q, v));

  ds2.reset(new Disk(2, 2, q, v));

  std::ofstream ofs("Kernelt7.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    siconos_io_register_Mechanics(oa);
    oa << NVP(ds1);
  }

  std::ifstream ifs("Kernelt7.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    siconos_io_register_Mechanics(ia);
    ia >> NVP(ds2);
  }

  CPPUNIT_ASSERT(std11::static_pointer_cast<Disk>(ds1)->getRadius() ==
                 std11::static_pointer_cast<Disk>(ds2)->getRadius());
}

void KernelTest::t8()
{
  SP::DynamicalSystem ds1, ds2;

  SP::SiconosVector q(new SiconosVector(3));
  SP::SiconosVector v(new SiconosVector(3));

  (*q)(0) = 0.;
  (*q)(1) = 1.;
  (*q)(2) = 1.;

  (*v)(0) = 0;
  (*v)(1) = 0;
  (*v)(2) = 10.;

  SP::NonSmoothDynamicalSystem nsds(new NonSmoothDynamicalSystem(0,10));

  ds1.reset(new Disk(1, 1, q, v));
  ds2.reset(new Disk(2, 2, q, v));

  nsds->insertDynamicalSystem(ds1);
  nsds->insertDynamicalSystem(ds2);

  MechanicsIO IO;

  SP::SimpleMatrix positions = IO.positions(*nsds);
  SP::SimpleMatrix velocities = IO.velocities(*nsds);

  //ids
  CPPUNIT_ASSERT((*positions)(0,0) == 1);
  CPPUNIT_ASSERT((*velocities)(0,0) == 1);

  CPPUNIT_ASSERT((*positions)(1,0) == 2);
  CPPUNIT_ASSERT((*velocities)(1,0) == 2);

  CPPUNIT_ASSERT((*positions)(0,1) == 0.);
  CPPUNIT_ASSERT((*velocities)(0,1) == 0.);
  CPPUNIT_ASSERT((*positions)(0,2) == 1.);
  CPPUNIT_ASSERT((*positions)(1,2) == 1.);
  CPPUNIT_ASSERT((*velocities)(0,3) == 10.);
  CPPUNIT_ASSERT((*velocities)(1,3) == 10.);

}
#endif
