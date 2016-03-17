
#include "testContact.hpp"

#include "Contactor.hpp"
#include "SiconosShape.hpp"
#include "BulletBroadphase.hpp"
#include "BodyDS.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ContactTest);

void ContactTest::setUp() {}
void ContactTest::tearDown() {}

void ContactTest::t1()
{
  try
  {
    // Set up a Siconos Mechanics environment:
    // A BodyDS with a contactor consisting of a single sphere.
    SP::BodyDS body(new BodyDS());
    SP::Contactor contactor(new Contactor());
    SP::SiconosSphere sphere(new SiconosSphere(0,0,10,1.0));
    contactor->addShape(sphere);
    body->setContactor(contactor);

    // A BodyDS with a contactor consisting of a plane
    SP::BodyDS body2(new BodyDS());
    SP::Contactor contactor2(new Contactor());
    SP::SiconosPlane plane(new SiconosPlane(0,0,1,0.0));
    contactor2->addShape(plane);
    body2->setContactor(contactor2);

    // Object to manage the Bullet implementation of broadphase
    SP::BulletBroadphase broadphase(new BulletBroadphase());

    // Build Bullet representational mirror of contactors
    broadphase->buildGraph(contactor);
    broadphase->buildGraph(contactor2);

    // Perform Bullet broadphase, generates IndexSet1
    broadphase->performBroadphase();

    // Update a property
    sphere->setRadius(0.5);

    // Check for dirty objects and update the graph
    broadphase->updateGraph();

    // Perform Bullet broadphase, generates IndexSet1
    broadphase->performBroadphase();
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}
