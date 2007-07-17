/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "PlatformTest.h"

#include "LinearDS.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// Class test registration
CPPUNIT_TEST_SUITE_REGISTRATION(PlatformTest);

void PlatformTest::setUp()
{
  sizePb = 3;
  x0 = new SimpleVector(sizePb);
  for (unsigned int i = 0; i < sizePb; i++)
    (*x0)(i) = i;
}

void PlatformTest::tearDown()
{
  delete x0;
}

void PlatformTest::testManualCreation1()
{
  // Create a new model
  double t0 = 0.0;
  double T = 10.0
             Model * m = new Model(t0, T);
  m.display();
  cout << " Model creation ok " << endl;
  cout << "------------------------------------------------------ " << endl;

  // Create a NonSmoothDynamicalSystem
  NonSmoothDynamicalSystem* nsds = new NonSmoothDynamicalSystem(true);
  cout << " NonSmoothDynamicalSystem creation ok " << endl;
  cout << "------------------------------------------------------ " << endl;

  // Create dynamical systems

  // general type
  DynamicalSystem *ds1 = new DynamicalSystem(1, sizePb, *x0, "DefaultPlugin:vectorField");

  // linear
  DynamicalSystem *ds2 = new LinearDS(2, 5, sv);

  // Lagrangian
  DynamicalSystem *ds3 = new LagrangianDS(3, sizePbLag, *q0, *v0, *mass);

  // Lagrangian linear and time invariant
  DynamicalSystem *ds4 = new LagrangianLinearTIDS(4, sizePbLag, *q0, *v0, *mass,
      "DefaultPlugin:computeFExt", *K, *C);

  nsds->addDynamicalSystem(ds1);

  nsds->addDynamicalSystem(ds2);

  nsds->addDynamicalSystem(ds3);

  nsds->addDynamicalSystem(ds4);

  m.setNonSmoothDynamicalSystemPtr(nsds);

  vector<int> vect(2);
  vector<DynamicalSystem*> dsVect;
  dsVect.push_back(ds);
  dsVect.push_back(ds2);

  Interaction* inter;
  inter = nsds->addInteraction(1, 2, &vect, &dsVect);
  inter->createComplementarityConditionNSL();

  inter = nsds->addInteraction(2, 2, &vect, &dsVect);
  inter->createNewtonImpactLawNSL(9.81);

  inter = nsds->addInteraction(3, 2, &vect, &dsVect);
  inter->createRelayNSL(1.2, 5.2);

  Strategy* str;
  str = m.createTimeStepping();

  TimeDiscretisation* td;

  OneStepIntegrator* osi;
  osi = str->addAdams(td, ds);
  osi = str->addLsodar(td, ds2);
  osi = str->addMoreau(td, ds3, 1.0);

  //OneStepNSProblem* onspb;
  //  onspb = str->createLCP();

  m.saveToDOMTree();
  if (!m.checkXMLDOMTree()) exit(0);
  cout << "PlatformTest >>> testManualCreation - manual construction of the platform ..................... OK\n ";
}

void PlatformTest::testManualCreation2()
{
  Model m(0, 10);
  m.display();

  NonSmoothDynamicalSystem* nsds;
  nsds = m.createNonSmoothDynamicalSystem(false);
  nsds->display();
  m.setNonSmoothDynamicalSystemPtr(nsds);

  cout << "=== creation des DSs ===" << endl;
  DynamicalSystem *ds1, *ds2;
  /*SiconosVector*/
  SimpleVector sv;

  SimpleVector q0(3);
  q0.zero();
  q0(0) = 1.0;

  SimpleVector v0(3);
  v0.zero();
  v0(0) = 1.0;

  SiconosMatrix mass(3, 3);
  mass(0, 0) = 1.0;
  mass(0, 1) = 0.0;
  mass(0, 2) = 0.0;
  mass(1, 0) = 0.0;
  mass(1, 1) = 1.0;
  mass(1, 2) = 0.0;
  mass(2, 0) = 0.0;
  mass(2, 1) = 0.0;
  mass(2, 2) = 1.0;


  SiconosMatrix K(3, 3);
  K.zero();
  SiconosMatrix mC(3, 3);
  mC.zero();

  //  ds1 = nsds->addLagrangianLinearTIDS(1, 3, &q0, &v0, &mass,"BallPlugin:ballFExt", &K, &mC);
  ds1 = new LagrangianLinearTIDS(1, 3, q0, v0, mass, "BallPlugin:ballFExt", K, mC);
  nsds->addDynamicalSystem(ds1);
  //ds1->setQ()


  SimpleVector q0_(1);
  q0_.zero();

  SimpleVector v0_(1);
  v0_.zero();

  SiconosMatrix mass_(1, 1);
  mass_(0, 0) = 1000.0;

  SiconosMatrix K_(1, 1);
  K_.zero();
  SiconosMatrix C_(1, 1);
  C_.zero();

  //  ds2 = nsds->addLagrangianLinearTIDS(2, 1, &q0_, &v0_, &mass_,"BallPlugin:groundFExt", &K_, &C_);
  ds2 = new LagrangianLinearTIDS(2, 1, q0_, v0_, mass_, "BallPlugin:groundFExt", K_, C_);
  nsds->addDynamicalSystem(ds2);

  cout << "=== creation des Interactions ===" << endl;
  vector<int> vect(1);
  vect[0] = 0;
  vector<DynamicalSystem*> dsVect;
  dsVect.push_back(ds1);
  dsVect.push_back(ds2);

  SiconosMatrix H(1, 4);
  H(0, 0) = 1.0;
  H(0, 1) = 0.0;
  H(0, 2) = 0.0 ;
  H(0, 3) = -0.0;
  SimpleVector b(1);
  b.zero();
  Interaction* inter;
  inter = nsds->addInteraction(1, 1, &vect, &dsVect);
  //inter->createLagrangianLinearR(&H, &b);
  //inter->createNewtonImpactLawNSL(0.9);


  Strategy* str;
  str = m.createTimeStepping();

  vector<int> tk(2);
  tk[0] = 0;
  tk[1] = 1;
  TimeDiscretisation* td;
  //td = str->createTimeDiscretisation(0.005, 1, &sv, 0.0, 0.0, true);

  OneStepIntegrator* osi;
  osi = str->addMoreau(td, ds1, 1.0);

  osi = str->addMoreau(td, ds2, 1.0);

  SiconosMatrix M(1, 1);
  M(0, 0) = 0.0;
  b(0) = 1.0;

  OneStepNSProblem* onspb;
  //onspb = str->createLCP();
  LCP* lcp = static_cast<LCP*>(onspb);
  lcp->setN(1);
  lcp->setM(M);
  lcp->setQ(b);

  m.saveToDOMTree();
  if (!m.checkXMLDOMTree()) exit(0);
  cout << "PlatformTest >>> testManualCreation2 - manual construction of the platform ..................... OK\n ";
}

void PlatformTest::testManualCreation3()
{
  Model mo(0, 10);
  mo.display();

  NonSmoothDynamicalSystem* nsds;
  nsds = mo.createNonSmoothDynamicalSystem(false);
  mo.setNonSmoothDynamicalSystemPtr(nsds);

  DynamicalSystem *ds, *ds2, *ds3, *ds4;
  SimpleVector sv;


  ds = new DynamicalSystem(1, 2, sv, "DefaultPlugin:vectorField");
  ds2 = new LinearDS(2, 5, sv);
  /*    ds3 = new LagrangianDS(3, 2, sv, sv, "DefaultPlugin:computeMass",
    "DefaultPlugin:computeFInt","DefaultPlugin:computeFExt","DefaultPlugin:computeJacobianQFInt",
    "DefaultPlugin:computeJacobianVelocityFInt","DefaultPlugin:computeJacobianQQNLInertia",
    "DefaultPlugin:computeJacobianVelocityQNLInertia", "DefaultPlugin:computeQNLInertia");
  */
  SiconosMatrix sm;
  ds4 = new LagrangianLinearTIDS(4, 8, sv, sv, sm, "DefaultPlugin:computeFExt", sm, sm);

  nsds->addDynamicalSystem(ds);
  nsds->addDynamicalSystem(ds2);
  nsds->addDynamicalSystem(ds3);
  nsds->addDynamicalSystem(ds4);

  vector<int> vect(2);
  vector<DynamicalSystem*> dsVect;
  dsVect.push_back(ds);
  dsVect.push_back(ds2);

  Strategy* str;
  str = mo.createTimeStepping();

  TimeDiscretisation* td;
  //td = str->createTimeDiscretisation(1.0, 1, &sv, 1.0, 1.0, true);

  OneStepIntegrator* osi;
  osi = str->addAdams(td, ds);
  osi = str->addLsodar(td, ds2);
  osi = str->addMoreau(td, ds3, 1.0);

  mo.saveToDOMTree();
  if (!mo.checkXMLDOMTree()) exit(0);
  cout << "PlatformTest >>> testManualCreation3 - manual construction of the platform ..................... OK\n ";
}

void PlatformTest::testPlatformException()
{
  Model m(10, 0);
  cout << "PlatformTest >>> testPlatformException - creation that must launch an exception ..................... OK\n ";
}

