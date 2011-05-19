#include "BulletBodies.hpp"
#include <BulletSpaceFilter.hpp>
#include <BulletTimeStepping.hpp>
#include <BulletDS.hpp>
#include <BulletWeightedShape.hpp>

#ifdef WITH_GLOBALAC
#include <mpi.h>
#endif

#include <stdlib.h>

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

    // -------------
    // --- Model ---
    // -------------
    _model.reset(new Model(t0, T));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator osi;
    osi.reset(new Moreau(0.5));

    SP::btCollisionShape box(new btBoxShape(btVector3(1, 1, 1)));
    SP::BulletWeightedShape box1(new BulletWeightedShape(box, 1.0));

    SP::btCollisionShape bbox(new btBoxShape(btVector3(.5, 1, 1)));
    SP::BulletWeightedShape box2(new BulletWeightedShape(bbox, 1.0));

    SP::btCollisionShape capsule(new btCapsuleShape(.5, .5));
    SP::BulletWeightedShape capsule1(new BulletWeightedShape(capsule, 1.0));

    SP::btCollisionShape bcapsule(new btCapsuleShape(1, 2));
    SP::BulletWeightedShape capsule2(new BulletWeightedShape(bcapsule, 2.0));

    SP::btCollisionShape cylinder(new btCylinderShape(btVector3(.5, .5, 1)));
    SP::BulletWeightedShape cylinder1(new BulletWeightedShape(cylinder, 1.0));

    std::vector<SP::BulletWeightedShape> shapes;
    shapes.push_back(box1);
    shapes.push_back(box1);
    shapes.push_back(box1);
    shapes.push_back(box1);
    //    shapes.push_back(box1);
    //    shapes.push_back(box2);
    //    shapes.push_back(capsule1);
    //    shapes.push_back(cylinder1);

    unsigned int N = 3;

    srand(1.);

    SP::SimpleVector FExt;
    FExt.reset(new SimpleVector(3)); //
    FExt->zero();
    FExt->setValue(2, -9.81 * box1->mass());

    for (unsigned int i = 0; i < N; ++i)
    {
      for (unsigned int j = 0; j < N; ++j)
      {
        for (unsigned int k = 0; k < N; ++k)
        {

          SP::SimpleVector position(new SimpleVector(7));
          SP::SimpleVector velocity(new SimpleVector(6));
          velocity->zero();
          position->zero();

          double theta = i + j + k; //acos(1/sqrt(3));

          double a = 1;
          double b = 1;
          double c = 0;
          double n = (sin(theta / 2)) / sqrt(a * a + b * b + c * c);
          (*position)(0) = 4.01 * i - 4.01 * (N - 1) / 2 + (double) rand() / (10.*RAND_MAX);
          (*position)(1) = 4.01 * j - 4.01 * (N - 1) / 2;
          (*position)(2) = 4.01 * k + 10;
          (*position)(3) = cos(theta / 2);
          (*position)(4) = a * n;
          (*position)(5) = b * n;
          (*position)(6) = c * n;

          velocity->zero();
          (*velocity)(3) = 0.;
          (*velocity)(4) = 0.;
          (*velocity)(5) = 0.;

          SP::BulletDS body(new BulletDS(shapes[(i + j + k) % 4],
                                         boost::shared_ptr<SimpleVector>(position),
                                         boost::shared_ptr<SimpleVector>(velocity)));
          body->setFExtPtr(FExt);

          _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);
          osi->insertDynamicalSystem(body);
        }
      }
    }


    SP::btCollisionObject ground(new btCollisionObject());
    SP::btCollisionShape groundShape(new btBoxShape(btVector3(30, 30, 3)));
    btMatrix3x3 basis;
    basis.setIdentity();
    ground->getWorldTransform().setBasis(basis);
    ground->setCollisionShape(&*groundShape);

    ground->getWorldTransform().getOrigin().setZ(-3.01);

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

    osnspb->numericsSolverOptions()->iparam[0] = 1000; // Max number of
    // iterations

#ifdef WITH_GLOBALAC
    osnspb->numericsSolverOptions()->iparam[1] = 1; // compute error with GLOBALAC
#else
    osnspb->numericsSolverOptions()->iparam[1] = 10; // with NSGS
#endif


    // iterations
    osnspb->numericsSolverOptions()->dparam[0] = 1e-5; // Tolerance

    osnspb->setMaxSize(16384);
    osnspb->setMStorageType(1);
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

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0, 0, .7, 3));


    SP::btVector3 aabbmax(new btVector3(100, 100, 100));
    SP::btVector3 aabbmin(new btVector3(-100, -100, -100));

    _playground.reset(new BulletSpaceFilter(_model->nonSmoothDynamicalSystem(),
                                            nslaw, aabbmin, aabbmax));

    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*ground);

    ask<ForStaticObjects>(*_playground)->push_back(ground);

    ask<ForStaticShapes>(*_playground)->push_back(groundShape);

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

