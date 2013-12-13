#include "BulletBodies.hpp"
#include <BulletSpaceFilter.hpp>
#include <BulletTimeStepping.hpp>
#include <BulletDS.hpp>
#include <BulletWeightedShape.hpp>

#include <stdlib.h>

#include <limits>


#define NDOF 3

void BulletBodies::init()
{
  SP::TimeDiscretisation timedisc;
  SP::BulletTimeStepping simulation;
  SP::FrictionContact osnspb;

  // User-defined main parameters


  double t0 = 0;                   // initial computation time

  double T = std::numeric_limits<double>::infinity();;

  double h = 0.005;


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

    SP::btCollisionShape capsule(new btCapsuleShape(.5, 2));
    SP::BulletWeightedShape capsule1(new BulletWeightedShape(capsule, 1.0));

    SP::btCollisionShape bcapsule(new btCapsuleShape(1, 2));
    SP::BulletWeightedShape capsule2(new BulletWeightedShape(bcapsule, 2.0));

    SP::btCollisionShape cylinder(new btCylinderShapeX(btVector3(1, .1, .1)));
    SP::BulletWeightedShape cylinder1(new BulletWeightedShape(cylinder, 1.0));


    SP::btCollisionShape sphere(new btSphereShape(0.5));
    //    sphere->setLocalScaling(btVector3(.5,10.,.2));
    SP::BulletWeightedShape sphere1(new BulletWeightedShape(sphere, 1.0));

    SP::btCollisionShape cone(new btConeShape(1, .1));
    SP::BulletWeightedShape cone1(new BulletWeightedShape(cone, 1.0));


    SP::btCollisionShape mshape(new btConvexHullShape());

    {
      double height = 1;
      double r = 1;
      std11::static_pointer_cast<btConvexHullShape>(mshape)->addPoint(btVector3(0.0,  0.75 * height, 0.0));
      std11::static_pointer_cast<btConvexHullShape>(mshape)->addPoint(btVector3(-r, -0.25 * height,    r));
      std11::static_pointer_cast<btConvexHullShape>(mshape)->addPoint(btVector3(r, -0.25 * height,    r));
      std11::static_pointer_cast<btConvexHullShape>(mshape)->addPoint(btVector3(r, -0.25 * height,   -r));
      std11::static_pointer_cast<btConvexHullShape>(mshape)->addPoint(btVector3(-r, -0.25 * height,   -r));
    }
    SP::BulletWeightedShape mshape1(new BulletWeightedShape(mshape, 1.0));


    const int numSpheres = 2;
    btVector3 positions[numSpheres] = {btVector3(0, 1, 0), btVector3(0, -1, 0)};
    btScalar radi[numSpheres] = {1, 1};
    SP::btCollisionShape mspheres(new btMultiSphereShape(positions, radi, 2));
    SP::BulletWeightedShape mspheres1(new BulletWeightedShape(mspheres, 1.0));


    std::vector<SP::BulletWeightedShape> shapes;
    shapes.push_back(box1);
    shapes.push_back(box1);
    shapes.push_back(box1);
    shapes.push_back(box1);
    //    shapes.push_back(box1);
    //    shapes.push_back(box2);
    //    shapes.push_back(capsule1);
    //    shapes.push_back(cylinder1);

    unsigned int N = 1;

    srand(1.);

    SP::SiconosVector FExt;
    FExt.reset(new SiconosVector(3)); //
    FExt->zero();
    FExt->setValue(2, -9.81 * box1->mass());

    for (unsigned int i = 0; i < N; ++i)
    {
      for (unsigned int j = 0; j < N; ++j)
      {
        for (unsigned int k = 0; k < 3; ++k)
        {

          SP::SiconosVector position(new SiconosVector(7));
          SP::SiconosVector velocity(new SiconosVector(6));
          velocity->zero();
          position->zero();

          double theta = 0.;//i+j+k;//acos(1/sqrt(3));

          double a = 1;
          double b = 1;
          double c = 0;
          double n = (sin(theta / 2)) / sqrt(a * a + b * b + c * c);
          (*position)(0) = 2.01 * i - 2.01 * (N - 1) / 2; // + (double) rand()/(10.*RAND_MAX);
          (*position)(1) = 2.01 * j - 2.01 * (N - 1) / 2;
          (*position)(2) = 2.01 * k + 2.;
          (*position)(3) = cos(theta / 2);
          (*position)(4) = a * n;
          (*position)(5) = b * n;
          (*position)(6) = c * n;

          velocity->zero();
          (*velocity)(3) = 0.;
          (*velocity)(4) = 0.;
          (*velocity)(5) = 0.;

          SP::BulletDS body(new BulletDS(shapes[(i + j + k) % 4],
                                         boost::shared_ptr<SiconosVector>(position),
                                         boost::shared_ptr<SiconosVector>(velocity)));
          body->setFExtPtr(FExt);

          _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);
          osi->insertDynamicalSystem(body);
        }
      }
    }

    ////////////////////////////////////////////////////////////////////////////////

    SP::btCollisionObject ground(new btCollisionObject());
    ground->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    SP::btCollisionObject wall1(new btCollisionObject());
    SP::btCollisionObject wall2(new btCollisionObject());
    SP::btCollisionObject wall3(new btCollisionObject());
    SP::btCollisionObject wall4(new btCollisionObject());
    ground->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    wall1->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    wall2->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    wall3->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    wall4->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);

    SP::btCollisionShape groundShape(new btBoxShape(btVector3(30, 30, .5)));
    SP::btCollisionShape wallShape(new btBoxShape(btVector3(3, 30, .5)));
    //    groundShape->setMargin(1.);
    //    wallShape->setMargin(1.);
    btMatrix3x3 basis;
    basis.setIdentity();

    ground->getWorldTransform().setBasis(basis);
    wall1->getWorldTransform().setBasis(basis);
    wall2->getWorldTransform().setBasis(basis);
    wall3->getWorldTransform().setBasis(basis);
    wall4->getWorldTransform().setBasis(basis);

    ground->setCollisionShape(&*groundShape);
    wall1->setCollisionShape(&*wallShape);
    wall2->setCollisionShape(&*wallShape);
    wall3->setCollisionShape(&*wallShape);
    wall4->setCollisionShape(&*wallShape);

    double s = sqrt(2.) / 2.;

    ground->getWorldTransform().getOrigin().setZ(-.51);
    wall1->getWorldTransform().getOrigin().setX(-30.01);
    wall1->getWorldTransform().getOrigin().setZ(-1.01);
    wall1->getWorldTransform().getBasis().setRotation(btQuaternion(0, s, 0, s));
    wall2->getWorldTransform().getOrigin().setY(-30.01);
    wall2->getWorldTransform().getOrigin().setZ(-1.01);
    wall2->getWorldTransform().getBasis().setRotation(btQuaternion(0, 0, s, s)*btQuaternion(0, s, 0, s));
    wall3->getWorldTransform().getOrigin().setX(30.01);
    wall3->getWorldTransform().getOrigin().setZ(-1.01);
    wall3->getWorldTransform().getBasis().setRotation(btQuaternion(0, s, 0, s));
    wall4->getWorldTransform().getOrigin().setY(30.01);
    wall4->getWorldTransform().getOrigin().setZ(-1.01);
    wall4->getWorldTransform().getBasis().setRotation(btQuaternion(0, 0, s, s)*btQuaternion(0, s, 0, s));
    // ------------------
    // --- Simulation ---
    // ------------------

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(.0, .0, .3, 3));


    SP::btVector3 aabbmax(new btVector3(100, 100, 100));
    SP::btVector3 aabbmin(new btVector3(-100, -100, -100));

    _playground.reset(new BulletSpaceFilter(_model,
                                            nslaw));

    // -- Time discretisation --
    timedisc.reset(new TimeDiscretisation(t0, h));

    simulation.reset(new BulletTimeStepping(timedisc, std11::static_pointer_cast<BulletSpaceFilter>(_playground)));

#ifdef WITH_SOLVER
    osnspb.reset(new FrictionContact(3, WITH_SOLVER));
#else
    osnspb.reset(new FrictionContact(3));
#endif

    osnspb->numericsSolverOptions()->iparam[0] = 1000; // Max number of
    // iterations

    osnspb->numericsSolverOptions()->dparam[0] = 1e-5; // Tolerance

    osnspb->setMaxSize(16384);
    osnspb->setMStorageType(1);                       // Sparse
    osnspb->setNumericsVerboseMode(0);
    osnspb->setKeepLambdaAndYState(true);
    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;



    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*ground);
    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*wall1);
    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*wall2);
    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*wall3);
    ask<ForCollisionWorld>(*_playground)->addCollisionObject(&*wall4);


    ask<ForStaticObjects>(*_playground)->push_back(ground);
    ask<ForStaticObjects>(*_playground)->push_back(wall1);
    ask<ForStaticObjects>(*_playground)->push_back(wall2);
    ask<ForStaticObjects>(*_playground)->push_back(wall3);
    ask<ForStaticObjects>(*_playground)->push_back(wall4);

    ask<ForStaticShapes>(*_playground)->push_back(groundShape);
    ask<ForStaticShapes>(*_playground)->push_back(wallShape);

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

