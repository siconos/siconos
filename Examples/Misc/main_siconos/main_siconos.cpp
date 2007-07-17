/* Siconos-sample version 2.1.1, Copyright INRIA 2005-2006.
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
#include "main_siconos.h"


extern "C" void essai_plugin(Model *m)
{
  OUT("--------- Tests d'integration : plateforme <=> plugins ----------------\n");
  m->getNonSmoothDynamicalSystem()->getDynamicalSystemOnNumber(1)->vectorField(123.3541);
  //  static_cast<LagrangianDS*>(m->getNonSmoothDynamicalSystem()->getDynamicalSystemOnNumber(1))->computeQNLInertia();
  OUT("-------------------------\n");
}

extern "C" void essai_model()
{
  OUT("-------------------------\n");
  OUT("OUT::creation Model m2\n");

  Model m2("xml_test.xml");
  OUT("-----------end createModel--------------\n");
  //m2.createModel("/home0/barbier/siconos/SICONOS/src/model/test/xml_schemaAttributeA.xml");
  m2.checkModelCoherency();

  //static methode_lcp meth_lcp  = {"gcp",101, 0.0001,0.6};

  OUT("-------------------------\n");
  m2.saveToXMLFile("xml_save.xml");
  OUT("-------------------------\n");
  cout << endl << "===================== THE END ========================" << endl;
}


extern "C" void essai_model_XML(char  xmlFile[])
{
  IN("OUT::creation Model m\n");
  //char xmlFile[50] = "../sample/xml_test.xml";

  //Model m;

  if ((strcmp(xmlFile, "creation") != 0) && (strcmp(xmlFile, "creation2") != 0) && (strcmp(xmlFile, "creation3") != 0)
      && (strcmp(xmlFile, "mixte") != 0) && (strcmp(xmlFile, "check") != 0) && (strcmp(xmlFile, "lmgc90") != 0))
  {
    Model m(xmlFile);
    char* saveFile = xmlFile;
    strcat(saveFile, ".save.xml");
    cout << "saveFile == " << saveFile << endl;
    m.saveToXMLFile(saveFile);
  }
  else if (strcmp(xmlFile, "lmgc90") == 0)
  {
    Model m("lmgc90.xml");
    char* saveFile = xmlFile;
    strcat(saveFile, ".save.xml");
    cout << "saveFile == " << saveFile << endl;
    m.saveToXMLFile(saveFile);
  }
  else if (strcmp(xmlFile, "mixte") == 0)
  {
    cout << "*** essai_model - creation de la plateforme avec fichier XML ET fichier de commande ***" << endl;
    Model m("xml_test.xml");

    NonSmoothDynamicalSystem* nsds;
    DynamicalSystem *ds, *ds2;
    SimpleVector sv(2);
    SiconosMatrix sm;
    cout << "=== Add of new objects to the Model ===" << endl;
    nsds = m.getNonSmoothDynamicalSystem();
    //ds = nsds->addNonLinearSystemDS( 11, 2, &sv, "BasicPlugin:vectorField" );
    ds = new DynamicalSystem(11, 2, &sv, "BasicPlugin:vectorField");
    ds->createPeriodicBC();

    Interaction* inter;
    vector<int> vect(2);
    vector<DynamicalSystem*> dsVect;
    dsVect.push_back(ds);
    ds2 = nsds->getDynamicalSystem(1);
    dsVect.push_back(ds2);

    SimpleVector sv_;
    inter = nsds->addInteraction(11, 11, &vect, &dsVect);
    inter->createLagrangianLinearR(&sm, &sv_);
    inter->createComplementarityConditionNSL();

    OneStepNSProblem *onspb = m.getStrategy()->getOneStepNSProblem();
    onspb->setLatinAlgorithm("LcpSolving", 0.00001, "max", 101, 0.6);

    m.saveToXMLFile("mixte.save.xml");

  }
  else if (strcmp(xmlFile, "creation") == 0)
  {
    clock_t t1 = clock();

    cout << "*** essai_model - creation complete de la plateforme depuis le fichier de commande ***" << endl;
    Model m(0, 10);
    m.display();

    NonSmoothDynamicalSystem* nsds;

    //nsds = m.createNonSmoothDynamicalSystem( false );
    nsds = new NonSmoothDynamicalSystem(false);
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
    //v0(0) = 0.0;

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
    K(0, 0) = 0.0;
    K(0, 1) = 0.0;
    K(0, 2) = 0.0;
    K(1, 0) = 0.0;
    K(1, 1) = 0.0;
    K(1, 2) = 0.0;
    K(2, 0) = 0.0;
    K(2, 1) = 0.0;
    K(2, 2) = 0.0;

    SiconosMatrix mC(3, 3);
    mC(0, 0) = 0.0;
    mC(0, 1) = 0.0;
    mC(0, 2) = 0.0;
    mC(1, 0) = 0.0;
    mC(1, 1) = 0.0;
    mC(1, 2) = 0.0;
    mC(2, 0) = 0.0;
    mC(2, 1) = 0.0;
    mC(2, 2) = 0.0;


    //ds1 = nsds->addLagrangianLinearTIDS(1, 3, &q0, &v0, &mass,"BallPlugin:ballFExt",&K, &mC);
    ds1 = new LagrangianLinearTIDS(1, 3, &q0, &v0, &mass, "BallPlugin:ballFExt", &K, &mC);
    nsds->addDynamicalSystem(ds1);
    //ds1->setQ()

    SimpleVector q0_(1);
    q0_.zero();

    SimpleVector v0_(1);
    v0_.zero();

    SiconosMatrix mass_(1, 1);
    mass_(0, 0) = 1.0;
    SiconosMatrix K_(1, 1);
    K_(0, 0) = 0.0;
    SiconosMatrix C_(1, 1);
    C_(0, 0) = 0.0;

    //    ds2 = nsds->addLagrangianLinearTIDS(2, 1, &q0_, &v0_, &mass_,"sample/BouncingBall/BallPlugin:groundFExt",
    //        &K_, &C_);
    ds2 = new LagrangianLinearTIDS(2, 1,
                                   &q0_, &v0_, &mass_,
                                   "BallPlugin:groundFExt",
                                   &K_, &C_);
    /*    static_cast<LagrangianLinearTIDS*>(ds2)->createDynamicalSystem(NULL, 2, 1,
                      &q0_, &v0_, &mass_,
                      "BallPlugin:groundFExt",
                      &K_, &C_); */
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
    cout << "# " << endl;

    OneStepNSProblem* onspb;
    onspb = str->createLCP();
    cout << "# " << endl;

    LCP* lcp = static_cast<LCP*>(onspb);
    lcp->setGsnlAlgorithm("LcpSolving", 0.0001, "max", 101);
    cout << "# " << endl;

    lcp->setN(1);
    lcp->setM(M);
    lcp->setQ(b);

    //============== comuptations ===============================
    try
    {
      cout << "\n *** formalisation completed ***" << endl;
      cout << "the strategy will be initialized" << endl;
      str->initialize();
      cout << "\n **** the strategy is ready ****" << endl;

      int k = td->getK();
      int N = td->getN();

      /*      // Trace Values
            SiconosMatrix mat(N+1, 6);
            //time
            mat(k, 0) = k*td->getH();
            // position
            LagrangianDS* ball = static_cast<LagrangianDS*> (m.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
            mat(k, 1) = (ball->getQ())(0);
                  // position
            mat(k, 2) = (ball->getVelocity())(0);
            LagrangianDS* ground = static_cast<LagrangianDS*> (m.getNonSmoothDynamicalSystem()->getDynamicalSystem(1));
            mat(k, 3) = (ground->getQ())(0);
                  // position
            mat(k, 4) = (ground->getVelocity())(0);

            SimpleVector vv;
            Interaction* I =  m.getNonSmoothDynamicalSystem()->getInteraction(0);
            I->display();

            mat(k,5) =0.0;
      */

      while (k < N)
      {
        str->nextStep();
        cout << "NextStep done" << endl;
        k = td->getK();

        cout << "iteration : " << k << endl;
        str->computeFreeState();
        str->formaliseOneStepNSProblem();
        str->computeOneStepNSProblem();
        str->updateState();

        /*        // Trace Values
                //time
                mat(k, 0) = k*td->getH();
                // position
                //LagrangianDS* ball = static_cast<LagrangianDS*> (bouncingBall.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
                mat(k, 1) = (ball->getQ())(0);
                      // position
                mat(k, 2) = (ball->getVelocity())(0);

                mat(k, 3) = (ground->getQ())(0);
                      // position
                mat(k, 4) = (ground->getVelocity())(0);

                mat(k, 5) = (m.getNonSmoothDynamicalSystem()->getInteraction(0)->getLambda())(0);
        */
        //ball->display();
        //ground->display();
        //        ball->getX()->display();
        //        cout<<"--- end of the step : "<<k<<" ---"<<endl;
        //        getchar();
      }
      cout << "iterations  done: " << k << endl;

      //      mat.write("result_test.dat", "ascii");

      //bouncingBall.saveToXMLFile("./BouncingBall_TIDS_test.xml.output");
      m.saveToXMLFile("save_totally_created_Model_test.xml");

      clock_t t2 = clock();
      printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);


    }
    catch (SiconosException e)
    {
      cout << e.report() << endl;
    }
    catch (...)
    {
      cout << "Exception caught in \'main_siconos\'" << endl;
    }
    //===========================================================

    //m.saveToXMLFile("save_totally_created_BB_test.xml");
  }
  else if (strcmp(xmlFile, "check") == 0)
  {
    //    cout<<"##############################"<<endl;
    //    m.xmlSchemaValidated("xml_test.xml", "config/xmlschema/SiconosModelSchema-V1.2.xsd");
    //    cout<<"##############################"<<endl;
  }
  else if (strcmp(xmlFile, "creation3") == 0)
  {
    Model m(0, 10);
    m.display();

    NonSmoothDynamicalSystem* nsds;
    nsds = m.createNonSmoothDynamicalSystem(false);
    m.setNonSmoothDynamicalSystem(nsds);

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
    str = m.createTimeStepping();

    TimeDiscretisation* td;
    td = str->createTimeDiscretisation(1.0, 1, &sv, 1.0, 1.0, true);

    OneStepIntegrator* osi;
    osi = str->addAdams(td, ds);
    osi = str->addLsodar(td, ds2);
    osi = str->addMoreau(td, ds3, 1.0);

    m.saveToDOMTree();
    if (!m.checkXMLDOMTree()) exit(0);
  }
  else
  {
    Model m(0, 10);
    m.display();
    m.setFileStorage(N_BINARY);

    NonSmoothDynamicalSystem* nsds;

    nsds = m.createNonSmoothDynamicalSystem(true);
    DynamicalSystem *ds, *ds2, *ds3, *ds4;
    SimpleVector sv;
    vector<double> test(20);
    for (int i = 0; i < 20; i++) test[i] = 1.23456;
    SimpleVector svv(test);


    //ds = nsds->addNonLinearSystemDS( 1, 2, &svv, "BasicPlugin:vectorField" );
    ds = new DynamicalSystem(1, 2, &svv, "BasicPlugin:vectorField");
    ds->createPeriodicBC();

    //ds2 = nsds->addLinearSystemDS( 2, 5, &sv );
    ds2 = new LinearSystemDS(2, 5, &sv);
    ds2->createNLinearBC();

    /*ds3 = nsds->addLagrangianDS(3, 2, &sv, &sv, "BasicPlugin:computeMass",
        "BasicPlugin:computeFInt","BasicPlugin:computeFExt","BasicPlugin:computeJacobianQFInt",
        "BasicPlugin:computeJacobianVelocityFInt","BasicPlugin:computeJacobianQQNLInertia",
        "BasicPlugin:computeJacobianVelocityQNLInertia", "BasicPlugin:computeQNLInertia");
    */
    ds3 = new LagrangianDS(3, 2, &sv, &sv, "BasicPlugin:computeMass",
                           "BasicPlugin:computeFInt", "BasicPlugin:computeFExt", "BasicPlugin:computeJacobianQFInt",
                           "BasicPlugin:computeJacobianVelocityFInt", "BasicPlugin:computeJacobianQQNLInertia",
                           "BasicPlugin:computeJacobianVelocityQNLInertia", "BasicPlugin:computeQNLInertia");
    ds3->createPeriodicBC();

    SiconosMatrix sm_(2, 2);
    sm_(0, 0) = 1.0;
    sm_(1, 0) = 1.0;
    sm_(0, 1) = 1.0;
    sm_(1, 1) = 1.0;
    DSInputOutput *dsio_;
    dsio_ = new LagrangianDSIO();
    static_cast<LagrangianDSIO*>(dsio_)->createDSInputOutput(NULL, 1);//, &sm_);
    ds3->addDSInputOutput(dsio_);

    SiconosMatrix sm;
    //    ds4 = nsds->addLagrangianLinearTIDS(4, 8, &sv, &sv, &sm,"BasicPlugin:computeFExt",
    //        &sm, &sm);
    ds4 = new LagrangianLinearTIDS(4, 8,
                                   &sv, &sv, &sm, "BasicPlugin:computeFExt", &sm, &sm);
    //    static_cast<LagrangianLinearTIDS*>(ds4)->createDynamicalSystem(NULL, 4, 8,
    //                &sv, &sv, &sm,"BasicPlugin:computeFExt",&sm, &sm);
    ds4->createLinearBC(&sv, &sm, &sm);

    DSInputOutput *dsio;
    dsio = new DSInputOutput();
    dsio->createDSInputOutput(NULL, 2);//, &sm);
    ds4->addDSInputOutput(dsio);

    vector<DSInputOutput*> dsioVec(2);
    dsioVec[0] = dsio_;
    dsioVec[1] = dsio;

    nsds->addDynamicalSystem(ds);
    nsds->addDynamicalSystem(ds2);
    nsds->addDynamicalSystem(ds3);
    nsds->addDynamicalSystem(ds4);

    EqualityConstraint *ec;
    ec = new EqualityConstraint();
    ec->createEqualityConstraint(NULL, 1, &sm, &dsioVec);
    nsds->addEqualityConstraint(ec);

    LinearTIEC *tiec;
    tiec = new LinearTIEC();
    tiec->createEqualityConstraint(NULL, 2, &sm, &dsioVec);
    nsds->addEqualityConstraint(tiec);

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

    //onspb->setSolvingMethod( "LcpSolving" );
    onspb->setLatinAlgorithm("LcpSolving", 0.0, "max", 10, 0.0);

    m.saveToXMLFile("save_totally_created_Model.xml");
  }
  cout << endl << "===================== THE END ========================" << endl;
}

extern "C" void essai_model2()
{
  OUT("-------------------------\n");
  IN("IN::creation Model m2\n");
  Model m2;
  //  m2.createModel(NULL, 1,0);
  OUT("OUT::creation Model m2\n");
  cout << endl << "===================== THE END ========================" << endl;
}

extern "C" void test_schema()
{
  IN("-------------------------\n");
  OUT("OUT::creation Model m\n");
  Model m2("test_schema.xml");
  OUT("-------------------------\n");
  m2.saveToXMLFile("xml_save.xml");
  OUT("-------------------------\n");
  cout << endl << "===================== THE END ========================" << endl;
}

void test_xmlfile(char * xmlFile)
{
  Model m(xmlFile);
  strcat(xmlFile, ".save.xml");
  m.saveToXMLFile(xmlFile);
  cout << " #_#_# file tested saved: " << xmlFile << endl;
}


void bench()
{
  /** declarations */
  const int size = 1000;

  /** time management */
  double t1, t2, elapsed;
  struct timeval tp;
  int rtn;
  rtn = gettimeofday(&tp, NULL);
  t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

  /** treatement */

  SiconosMatrix M(size, size);
  M.zero();
  SimpleVector v;
  for (int i = 0; i < 1000; i++)
  {
    // cout<<"i = "<<i<<endl;
    v = M.getRow(size / 2);
    v.zero();
  }
  /** time management */

  rtn = gettimeofday(&tp, NULL);
  t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
  elapsed = t2 - t1;

  cout << "time = " << elapsed << endl;
}

int main(int argc, char* argv[])
{
  try
  {


    cout << "====================================================== argc == " << argc << endl << endl;

    if (argc == 2)
    {
      essai_model_XML(argv[1]);
    }
    else if (argc == 3)
    {
      test_xmlfile(argv[2]);
    }
    else essai_model();

    cout << "======================================================" << endl << endl;

    bench();

    cout << "======================================================" << endl << endl;

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  //  catch(...)
  //  {
  //    cout << "Exception caught in \'main_siconos\'" << endl;
  //  }
}

