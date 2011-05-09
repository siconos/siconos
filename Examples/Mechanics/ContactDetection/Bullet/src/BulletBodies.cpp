#include "BulletBodies.hpp"
#include <BulletSpaceFilter.hpp>
#include <BulletTimeStepping.hpp>
#include <BulletDS.hpp>
#include <BulletWeightedShape.hpp>

#ifdef WITH_GLOBALAC
#include <mpi.h>
#endif

#define NDOF 3

void BulletBodies::init()
{
  SP::TimeDiscretisation timedisc;
  SP::BulletTimeStepping simulation;
  SP::FrictionContact osnspb;

  // User-defined main parameters


  double t0 = 0;                   // initial computation time

  double T = 0.02;

  double h = 0.01;                // time step


  // -----------------------------------------
  // --- Dynamical systems && interactions ---
  // -----------------------------------------

  try
  {

    // ------------
    // --- Init ---
    // ------------

    std::cout << "====> Model loading ..." << std::endl << std::endl;


    double theta = 0.;

    double a = 1;
    double b = 0;
    double c = 0;
    double k = (sin(theta / 2)) / sqrt(a * a + b * b + c * c);

    SP::SimpleVector position(new SimpleVector(7));
    SP::SimpleVector velocity(new SimpleVector(6));
    position->zero();
    (*position)(1) = 0.;
    (*position)(2) = 10.;
    (*position)(3) = cos(theta / 2);
    (*position)(4) = a * k;
    (*position)(5) = b * k;
    (*position)(6) = c * k;

    velocity->zero();
    (*velocity)(3) = 0.;
    (*velocity)(4) = 0.;
    (*velocity)(5) = 0.;

    SP::btCollisionShape box(new btBoxShape(btVector3(1, 1, 1)));
    SP::BulletWeightedShape box1(new BulletWeightedShape(box, 1.0));

    SP::BulletDS body(new BulletDS(box1, position, velocity));
    SP::SimpleVector FExt;
    FExt.reset(new SimpleVector(3)); //
    FExt->zero();
    FExt->setValue(2, -9.81);
    body->setFExtPtr(FExt);

    double theta2 = acos(1 / sqrt(3)) + 0.01;

    double a2 = 1;
    double b2 = 0;
    double c2 = 0;
    double k2 = (sin(theta2 / 2)) / sqrt(a2 * a2 + b2 * b2 + c2 * c2);

    SP::SimpleVector position2(new SimpleVector(7));
    SP::SimpleVector velocity2(new SimpleVector(6));
    position2->zero();
    (*position2)(1) = .3;
    (*position2)(2) = 20.;
    (*position2)(3) = cos(theta2 / 2);
    (*position2)(4) = a2 * k2;
    (*position2)(5) = b2 * k2;
    (*position2)(6) = c2 * k2;

    velocity2->zero();
    (*velocity2)(3) = 0.;
    (*velocity2)(4) = 0.;
    (*velocity2)(5) = 0.;

    SP::BulletDS body2(new BulletDS(box1, position2, velocity2));
    body2->setFExtPtr(FExt);

    double theta3 = acos(1 / sqrt(3)) + 0.10;

    double a3 = 1;
    double b3 = 0;
    double c3 = 0;
    double k3 = (sin(theta3 / 2)) / sqrt(a3 * a3 + b3 * b3 + c3 * c3);

    SP::SimpleVector position3(new SimpleVector(7));
    SP::SimpleVector velocity3(new SimpleVector(6));
    position3->zero();
    (*position3)(1) = 0.;
    (*position3)(2) = 30.;
    (*position3)(3) = cos(theta3 / 2);
    (*position3)(4) = a3 * k3;
    (*position3)(5) = b3 * k3;
    (*position3)(6) = c3 * k3;

    velocity3->zero();
    (*velocity3)(3) = 0.;
    (*velocity3)(4) = 0.;
    (*velocity3)(5) = 0.;

    SP::BulletDS body3(new BulletDS(box1, position3, velocity3));
    body3->setFExtPtr(FExt);

    SP::btCollisionObject ground(new btCollisionObject());
    SP::btCollisionShape groundShape(new btBoxShape(btVector3(50, 50, 3)));
    btMatrix3x3 basis;
    basis.setIdentity();
    ground->getWorldTransform().setBasis(basis);
    ground->setCollisionShape(&*groundShape);

    ground->getWorldTransform().getOrigin().setZ(-3);

    SP::btCollisionObject fbox(new btCollisionObject());
    SP::btCollisionShape fboxShape(new btBoxShape(btVector3(1, 1, 1)));
    fbox->getWorldTransform().setBasis(basis);
    fbox->setCollisionShape(&*fboxShape);

    fbox->getWorldTransform().getOrigin().setZ(5);

    // -------------
    // --- Model ---
    // -------------
    _model.reset(new Model(t0, T));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator osi;
    osi.reset(new Moreau(0.5));


    _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);
    _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body2);
    //    _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body3);

    osi->insertDynamicalSystem(body);
    osi->insertDynamicalSystem(body2);
    //    osi->insertDynamicalSystem(body3);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    timedisc.reset(new TimeDiscretisation(t0, h));

    simulation.reset(new BulletTimeStepping(timedisc));

#ifdef WITH_GLOBALAC
    osnspb.reset(new FrictionContact(3, SICONOS_FRICTION_3D_GLOBALAC));
#else
    osnspb.reset(new FrictionContact(3));
#endif

    osnspb->numericsSolverOptions()->iparam[0] = 10000; // Max number of
    // iterations

#ifdef WITH_GLOBALAC
    osnspb->numericsSolverOptions()->iparam[1] = 1; // compute error with GLOBALAC
#else
    osnspb->numericsSolverOptions()->iparam[1] = 10; // with NSGS
#endif


    // iterations
    osnspb->numericsSolverOptions()->dparam[0] = 1e-7; // Tolerance

    osnspb->setMaxSize(16384);
    osnspb->setMStorageType(0);
    osnspb->setNumericsVerboseMode(0);
    osnspb->setKeepLambdaAndYState(true);
    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);

#ifdef WITH_GLOBALAC
    int myid, ierr;
    int argc = 0;
    char **argv = 0;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    osnspb->numericsSolverOptions()->iparam[4] = myid;
    osnspb->numericsSolverOptions()->iparam[5] = 1;

    frictionContact3D_sparseGlobalAlartCurnierInit(osnspb->numericsSolverOptions().get());
#endif
    //simulation_->setCheckSolverFunction(localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0, 0, 0.2, 3));

    _playground.reset(new BulletSpaceFilter(_model->nonSmoothDynamicalSystem(),
                                            nslaw));

    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*ground);

    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*fbox);

    ask<ForStaticObjects>(*_playground)->push_back(ground);

    ask<ForStaticObjects>(*_playground)->push_back(fbox);

    ask<ForStaticShapes>(*_playground)->push_back(groundShape);

    ask<ForStaticShapes>(*_playground)->push_back(fboxShape);

    SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();


    _model->nonSmoothDynamicalSystem()->setSymmetric(true);

    _model->initialize(simulation);

    simulation->updateWorldFromDS();
    _playground->buildInteractions(_model->currentTime());

    std::cout << "====> End of initialisation ..." << std::endl << std::endl;

  }

  catch (SiconosException e)
  {
    std::cout << e.report() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cout << "Exception caught in BulletBodies::init()" << std::endl;
    exit(1);
  }

}

