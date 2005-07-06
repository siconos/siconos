#include "PlatformTest.h"

#include "LinearSystemDS.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// Class test registration
CPPUNIT_TEST_SUITE_REGISTRATION(PlatformTest);

void PlatformTest::setUp()
{}

void PlatformTest::tearDown()
{}

//______________________________________________________________________________

void PlatformTest::testOptionalAttributes1()
{
  /*
   * test with all the attributes
   */
  Model m("xml_test.xml");
  m.saveToXMLFile("xml_test.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes1 - XML file with all the attributes ......................... OK\n ";
}

void PlatformTest::testOptionalAttributes2()
{
  /*
   * test without unnecessary memory objects
   */
  Model m2("xml_uncomplete1.xml");
  m2.saveToXMLFile("xml_uncomplete1.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes2 - XML file with no memory objects ........................ OK\n ";
}

void PlatformTest::testOptionalAttributes3()
{
  /*
   * test with no attribute 'id' for DynamicalSystem and Interaction
   */
  Model m3("xml_uncomplete2.xml");
  m3.saveToXMLFile("xml_uncomplete2.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes3 - XML file with no object's id ........................ OK\n ";
}

void PlatformTest::testOptionalAttributes4()
{
  /*
   * test without Strategy
   */
  Model m4("xml_uncomplete3.xml");
  m4.saveToXMLFile("xml_uncomplete3.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes4 - XML file without Strategy ....................... OK\n ";
}

void PlatformTest::testOptionalAttributes5()
{
  /*
   * test without 'computeInput' and 'computeOutput' for the Relations
   */
  Model m5("xml_uncomplete4.xml");
  m5.saveToXMLFile("xml_uncomplete4.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes5 - XML without computeInput and computeOutput ................... OK\n ";
}

void PlatformTest::testOptionalAttributes6()
{
  /*
   * test with empty SiconosMemory
   */
  Model m6("xml_uncomplete5.xml");
  m6.saveToXMLFile("xml_uncomplete5.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes6 - XML file with empty SiconosMemory ....................... OK\n ";
}

void PlatformTest::testOptionalAttributes7()
{
  /*
   * test without t (current time), and T (final time for the simulation)
   */
  Model m7("xml_uncomplete6.xml");
  //m7.saveToXMLFile("xml_uncomplete6.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes7 - XML file without t and T ........................ OK\n ";
}

void PlatformTest::testOptionalAttributes8()
{
  /*
   * test without interactions
   */
  Model m8("xml_uncomplete7.xml");
  //m8.saveToXMLFile("xml_uncomplete7.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes8 - XML file without Interaction ........................ OK\n ";
}

void PlatformTest::testOptionalAttributes9()
{
  /*
   * test with no OneStepNSProblem
   */
  Model m9("xml_uncomplete8.xml");
  m9.saveToXMLFile("xml_uncomplete8.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes9 - XML without OneStepNSProblem ........................... OK\n ";
}

void PlatformTest::testOptionalAttributes10()
{
  /*
   * test with a NonLinearSystemDS
   */
  Model m10("xml_uncomplete9.xml");
  m10.saveToXMLFile("xml_uncomplete9.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes10 - XML file with a NonLinearSystemDS ..................... OK\n ";
}

void PlatformTest::testOptionalAttributes11()
{
  /*
   * test with no XML, but createModel with parameters
   */
  Model m11(0, 10);

  cout << "PlatformTest >>> testOptionalAttributes11 - createModel with parameters ....................... OK\n ";
}

void PlatformTest::testOptionalAttributes12()
{
  /*
   * test with no x, xDot, q and velocity. Only initial data (data at t0)
   */
  Model m12("xml_uncomplete10.xml");
  m12.saveToXMLFile("xml_uncomplete10.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes12 - Only initial data (data at t0) ....................... OK\n ";
}

void PlatformTest::testOptionalAttributes13()
{
  /*
   * test with no x definition for Lagrangian Systems
   */
  Model m13("xml_uncomplete11.xml");
  m13.saveToXMLFile("xml_uncomplete11.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes13 - no x definition for Lagrangian Systems ....................... OK\n ";
}

void PlatformTest::testOptionalAttributes14()
{
  /*
   * test with no optional attributes for the DynamicalSystems
   */
  Model m14("xml_uncomplete12.xml");
  m14.saveToXMLFile("xml_uncomplete12.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes14 - no optional attributes for the DynamicalSystems ....................... OK\n ";
}

void PlatformTest::testOptionalAttributes15()
{
  /*
   * test with no optional attributes for the Interactions
   */
  Model m15("xml_uncomplete13.xml");
  m15.saveToXMLFile("xml_uncomplete13.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes15 - no optional attributes for the Interactions ....................... OK\n ";
}

void PlatformTest::testOptionalAttributes16()
{
  /*
   * test with no optional attributes for the Strategy
   */
  Model m16("xml_uncomplete14.xml");
  m16.saveToXMLFile("xml_uncomplete14.save.xml");

  cout << "PlatformTest >>> testOptionalAttributes16 - no optional attributes for the Strategy ....................... OK\n ";
}

void PlatformTest::testManualCreation()
{
  Model m(0, 10);
  m.display();

  NonSmoothDynamicalSystem* nsds;
  nsds = new NonSmoothDynamicalSystem(true);
  m.setNonSmoothDynamicalSystemPtr(nsds);

  DynamicalSystem *ds, *ds2, *ds3, *ds4;
  SimpleVector sv;

  ds = new DynamicalSystem(1, 2, &sv, "BasicPlugin:vectorField");
  //ds = nsds->addNonLinearSystemDS( 1, 2, &sv, "BasicPlugin:vectorField" );
  ds->createPeriodicBC();
  nsds->addDynamicalSystem(ds);

  //ds2 = nsds->addLinearSystemDS( 2, 5, &sv );
  ds2 = new LinearSystemDS(2, 5, &sv);
  ds2->createNLinearBC();
  nsds->addDynamicalSystem(ds2);

  /*  ds3 = nsds->addLagrangianDS(3, 2, &sv, &sv, "BasicPlugin:computeMass",
        "BasicPlugin:computeFInt","BasicPlugin:computeFExt","BasicPlugin:computeJacobianQFInt",
        "BasicPlugin:computeJacobianVelocityFInt","BasicPlugin:computeJacobianQQNLInertia",
        "BasicPlugin:computeJacobianVelocityQNLInertia", "BasicPlugin:computeQNLInertia");
  */
  ds3 = new LagrangianDS(3, 2, &sv, &sv, "BasicPlugin:computeMass",
                         "BasicPlugin:computeFInt", "BasicPlugin:computeFExt", "BasicPlugin:computeJacobianQFInt",
                         "BasicPlugin:computeJacobianVelocityFInt", "BasicPlugin:computeJacobianQQNLInertia",
                         "BasicPlugin:computeJacobianVelocityQNLInertia", "BasicPlugin:computeQNLInertia");
  ds3->createPeriodicBC();
  nsds->addDynamicalSystem(ds3);

  SiconosMatrix sm;
  //  ds4 = nsds->addLagrangianLinearTIDS(4, 8, &sv, &sv, &sm,"BasicPlugin:computeFExt",
  //      &sm, &sm);
  ds4 = new LagrangianLinearTIDS(4, 8, &sv, &sv, &sm, "BasicPlugin:computeFExt", &sm, &sm);
  ds4->createLinearBC(&sv, &sm, &sm);
  nsds->addDynamicalSystem(ds4);

  vector<int> vect(2);
  vector<DynamicalSystem*> dsVect;
  dsVect.push_back(ds);
  dsVect.push_back(ds2);

  Interaction* inter;
  inter = nsds->addInteraction(1, 2, &vect, &dsVect);
  inter->createLagrangianLinearR(&sm, &sv);
  inter->createComplementarityConditionNSL();

  inter = nsds->addInteraction(2, 2, &vect, &dsVect);
  inter->createLagrangianNonLinearR("BasicPlugin:computeInput", "BasicPlugin:computeOutput");
  inter->createNewtonImpactLawNSL(9.81);

  inter = nsds->addInteraction(3, 2, &vect, &dsVect);
  inter->createLinearTIR(&sm, &sm, &sm, &sv);
  inter->createRelayNSL(1.2, 5.2);

  Strategy* str;
  str = m.createTimeStepping();

  TimeDiscretisation* td;
  td = str->createTimeDiscretisation(1.0, 1, &sv, 1.0, 1.0, true);

  OneStepIntegrator* osi;
  osi = str->addAdams(td, ds);
  osi = str->addLsodar(td, ds2);
  osi = str->addMoreau(td, ds3, 1.0);

  OneStepNSProblem* onspb;
  onspb = str->createLCP();

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
  m.setNonSmoothDynamicalSystem(nsds);

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
  ds1 = new LagrangianLinearTIDS(1, 3, &q0, &v0, &mass, "BallPlugin:ballFExt", &K, &mC);
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
  ds2 = new LagrangianLinearTIDS(2, 1, &q0_, &v0_, &mass_, "BallPlugin:groundFExt", &K_, &C_);
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
  inter->createLagrangianLinearR(&H, &b);
  inter->createNewtonImpactLawNSL(0.9);


  Strategy* str;
  str = m.createTimeStepping();

  vector<int> tk(2);
  tk[0] = 0;
  tk[1] = 1;
  TimeDiscretisation* td;
  td = str->createTimeDiscretisation(0.005, 1, &sv, 0.0, 0.0, true);

  OneStepIntegrator* osi;
  osi = str->addMoreau(td, ds1, 1.0);

  osi = str->addMoreau(td, ds2, 1.0);

  SiconosMatrix M(1, 1);
  M(0, 0) = 0.0;
  b(0) = 1.0;

  OneStepNSProblem* onspb;
  onspb = str->createLCP();
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
  mo.setNonSmoothDynamicalSystem(nsds);

  DynamicalSystem *ds, *ds2, *ds3, *ds4;
  SimpleVector sv;


  ds = new DynamicalSystem(1, 2, &sv, "BasicPlugin:vectorField");
  ds2 = new LinearSystemDS(2, 5, &sv);
  ds3 = new LagrangianDS(3, 2, &sv, &sv, "BasicPlugin:computeMass",
                         "BasicPlugin:computeFInt", "BasicPlugin:computeFExt", "BasicPlugin:computeJacobianQFInt",
                         "BasicPlugin:computeJacobianVelocityFInt", "BasicPlugin:computeJacobianQQNLInertia",
                         "BasicPlugin:computeJacobianVelocityQNLInertia", "BasicPlugin:computeQNLInertia");
  SiconosMatrix sm;
  ds4 = new LagrangianLinearTIDS(4, 8, &sv, &sv, &sm, "BasicPlugin:computeFExt", &sm, &sm);

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
  td = str->createTimeDiscretisation(1.0, 1, &sv, 1.0, 1.0, true);

  OneStepIntegrator* osi;
  osi = str->addAdams(td, ds);
  osi = str->addLsodar(td, ds2);
  osi = str->addMoreau(td, ds3, 1.0);

  mo.saveToDOMTree();
  if (!mo.checkXMLDOMTree()) exit(0);
  cout << "PlatformTest >>> testManualCreation3 - manual construction of the platform ..................... OK\n ";
}

void PlatformTest::testManualCreation4()
{
  cout << "PlatformTest >>> testManualCreation4 - manual construction of the platform ..................... OK\n ";
}

void PlatformTest::testCheckXMLPlatform()
{
  cout << "PlatformTest >>> testCheckXMLPlatform - with manual construction of the platform ................. OK\n ";
}

void PlatformTest::testMixteCreation()
{
  Model m("xml_test.xml");

  NonSmoothDynamicalSystem* nsds;
  DynamicalSystem *ds, *ds2;
  SimpleVector sv;
  SiconosMatrix sm;
  nsds = m.getNonSmoothDynamicalSystem();
  //ds = nsds->addNonLinearSystemDS( 11, 2, &sv, "BasicPlugin:vectorField" );
  ds = new DynamicalSystem(11, 2, &sv, "BasicPlugin:vectorField");
  ds->createPeriodicBC();
  nsds->addDynamicalSystem(ds);

  Interaction* inter;
  vector<int> vect(2);
  vector<DynamicalSystem*> dsVect;
  dsVect.push_back(ds);
  ds2 = nsds->getDynamicalSystem(1);
  dsVect.push_back(ds2);

  inter = nsds->addInteraction(11, 11, &vect, &dsVect);
  inter->createLagrangianLinearR(&sm, &sv);
  inter->createComplementarityConditionNSL();

  //  m.saveToDOMTree();
  //  if( !m.checkXMLDOMTree() ) exit(0);
  m.saveToXMLFile("mixte.save.xml");
  cout << "PlatformTest >>> testMixteCreation - construction of the platform with XML and manual adds ..................... OK\n ";
}

void PlatformTest::testXMLSchemaAttributeGood()
{
  /*
   * test with a good xml schema attribute
   */
  Model m("xml_schemaAttribute.xml");

  cout << "PlatformTest >>> testXMLSchemaAttributeGood - XML file with a good xml schema attribute file name ..................... OK\n ";
}

void PlatformTest::testXMLSchemaAttributeGood2()
{
  /*
   * test with a good xml schema attribute
   */
  Model m("xml_schemaAttributeA.xml");

  cout << "PlatformTest >>> testXMLSchemaAttributeGood2 - XML file with a good xml schema attribute file name (with spaces in the string) ..................... OK\n ";
}

void PlatformTest::testXMLSchemaAttributeBad1()
{
  /*
   * test with a bad xml schema attribute
   */
  Model m("xml_schemaAttribute2.xml");

  cout << "PlatformTest >>> testXMLSchemaAttributeGBad1 - XML file with a bad xml schema attribute (empty field) ..................... OK\n ";
}

void PlatformTest::testXMLSchemaAttributeBad2()
{
  /*
   * test with a bad xml schema attribute
   */
  Model m("xml_schemaAttribute3.xml");

  cout << "PlatformTest >>> testXMLSchemaAttributeBad2 - XML file with a bad xml schema attribute 2 (unexisting file) ..................... OK\n ";
}

void PlatformTest::testPlatformException()
{
  Model m(10, 0);
  cout << "PlatformTest >>> testPlatformException - creation that must launch an exception ..................... OK\n ";
}

void PlatformTest::testOSNSP1()
{
  /*
   * test with OneStepNSProblem : LCP + LcpSolving + gsnl
   */
  Model m("xml_osnsp1.xml");

  cout << "PlatformTest >>> testOSNSP1 - XML file with OneStepNSProblem : LCP + LcpSolving + gsnl ..................... OK\n ";
}

void PlatformTest::testOSNSP2()
{
  /*
   * test with OneStepNSProblem : LCP + RelayPrimalSolving + gcp
   */
  Model m("xml_osnsp2.xml");

  cout << "PlatformTest >>> testOSNSP2 - XML file with OneStepNSProblem : LCP + RelayPrimalSolving + gcp ..................... OK\n ";
}

void PlatformTest::testOSNSP3()
{
  /*
   * test with OneStepNSProblem : LCP + RelayDualSolving + lemke
   */
  Model m("xml_osnsp3.xml");

  cout << "PlatformTest >>> testOSNSP3 - XML file with OneStepNSProblem : LCP + RelayDualSolving + lemke ..................... OK\n ";
}

void PlatformTest::testOSNSP4()
{
  /*
   * test with OneStepNSProblem : LCP + ContactFrictionPrimalSolving + gcp
   */
  Model m("xml_osnsp4.xml");

  cout << "PlatformTest >>> testOSNSP4 - XML file with OneStepNSProblem : LCP + ContactFrictionPrimalSolving + gcp ..................... OK\n ";
}

void PlatformTest::testOSNSP5()
{
  /*
   * test with OneStepNSProblem : LCP + ContactFrictionDualSolving + latin
   */
  Model m("xml_osnsp5.xml");

  cout << "PlatformTest >>> testOSNSP5 - XML file with OneStepNSProblem : LCP + ContactFrictionDualSolving + latin ..................... OK\n ";
}

void PlatformTest::testLmgc90()
{
  /*
   * test with OneStepNSProblem : LCP + ContactFrictionDualSolving + latin
   */
  Model m("lmgc90.xml");

  cout << "PlatformTest >>> testLmgc90 - XML file with lmgc90 joker flags : ..................... OK\n ";
}

void PlatformTest::testMainSiconos()
{
  /*
   *  unused tests
   */

  /*  int i = system ("cd ../../..; ./siconos/SICONOS");
    CPPUNIT_ASSERT_MESSAGE("test SICONOS ", i != -1);

    i = system ("cd ../../..; ./siconos/SICONOS creation2");
    CPPUNIT_ASSERT_MESSAGE("test SICONOS creation2 ", i != -1);

    i = system ("cd ../../..; ./siconos/SICONOS mixte");
    CPPUNIT_ASSERT_MESSAGE("test SICONOS mixte ", i != -1);

    i = system ("cd ../../..; ./siconos/SICONOS lmgc90");
    CPPUNIT_ASSERT_MESSAGE("test SICONOS lmgc90 ", i != -1);
  */
  printf("PlatformTest >>> testMainSiconos ................................ OK\n");
}

