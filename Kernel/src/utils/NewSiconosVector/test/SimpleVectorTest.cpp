//$Id: SimpleVectorTest.cpp,v 1.12 2005/01/18 11:00:26 jbarbier Exp $

#include "SimpleVectorTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SimpleVectorTest);


void SimpleVectorTest::setUp()
{
  int i;
  vector<double> vq(5);
  vector<double> vdotq(5);

  for (i = 0; i < 5; i++)
  {
    vq.at(i) = 1;
    vdotq.at(i) = 2;
  }
}

void SimpleVectorTest::tearDown()
{ }

//______________________________________________________________________________

void SimpleVectorTest::testBuildSimpleVector()
{
  SimpleVector v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", v.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", v.size() == 0, true);

  cout << "SimpleVectorTest >>> testBuildSimpleVector ............................... OK\n ";
}

void SimpleVectorTest::testBuildSimpleVector1()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  SimpleVector v(vq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v.size() == vq.size(), true);

  vq.clear();
  SimpleVector(v1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v1.size() == 0, true);

  cout << "SimpleVectorTest >>> testBuildSimpleVector1 ............................... OK\n ";
}


void SimpleVectorTest::testBuildSimpleVector2()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq);
  SimpleVector sv(vq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv == v, true);

  SimpleVector sv1(sv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv1.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv1 == v, true);

  // with pointers
  SimpleVector *pv;
  pv = new SimpleVector(v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", pv->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", pv->size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", *pv == v, true);

  SimpleVector sv2(*pv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv2.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv2.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv2 == v, true);
  delete pv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv2.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv2.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", sv2 == v, true);

  SiconosVector *pv1 = new SimpleVector(v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", pv1->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", pv1->size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", *pv1 == v, true);
  delete pv1;

  cout << "SimpleVectorTest >>> testBuildSimpleVector2 ............................... OK\n ";
}



void SimpleVectorTest::testBuildSimpleVector3()
{
  const int SIZE = 10;
  SimpleVector v(SIZE);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v.size() == SIZE, true);

  SimpleVector sv(0);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", sv.size() == 0 , true);

  cout << "SimpleVectorTest >>> testBuildSimpleVector3 ............................... OK\n ";
}


void SimpleVectorTest::testOperatorComp()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  //test vecteurs vide
  SimpleVector sv, sv1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv.isComposite() ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv.size() == 0 ", sv.size() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv1 ", sv == sv1, true);

  // test auto comparaison
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv", sv == sv, true);

  // test vecteur avec valeur == vide
  sv.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv.isComposite() ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv.size() == vq.size() ", sv.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv1 ", sv == sv1, false);
  // tests vect vide == vect avec valeurs
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv1 == sv ", sv1 == sv, false);

  // tests vect taille diff
  vector<double> vq1(3);
  for (int i = 0; i < 3; i++)
  {
    vq1[i] = 3;
  }
  sv1.setValues(vq1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.size() == vq1.size() ", sv1.size() == vq1.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv1 ", sv == sv1, false);

  // test vect meme taille 1er nb différent
  vq[0] = 0;
  sv1.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.size() == vq.size() ", sv1.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv1 ", sv == sv1, false);

  // test vect meme taille dernier nb different
  vq[4] = 5;
  vq[1] = 0;
  sv1.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.size() == vq.size() ", sv1.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv1 ", sv == sv1, false);

  // tests vect taille 1
  vector<double> vq2(1);
  vq2[0] = 100.365;
  sv.setValues(vq2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv.isComposite() ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv.size() == vq2.size() ", sv.size() == vq2.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv1 ", sv == sv1, false);

  sv1.setValues(vq2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp  : sv1.size() == vq2.size() ", sv1.size() == vq2.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : sv == sv1 ", sv == sv1, true);

  cout << "SimpleVectorTest >>> testOperatorComp ............................... OK\n ";
}


void SimpleVectorTest::testOperatorCompDiff()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  //test vecteurs vide
  SimpleVector sv, sv1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv.isComposite() ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv.size() == 0 ", sv.size() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv1 ", sv != sv1, false);

  // test auto comparaison
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv", sv != sv, false);

  // test vecteur avec valeur == vide
  sv.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv.isComposite() ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv.size() == vq.size() ", sv.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv1 ", sv != sv1, true);
  // tests vect vide == vect avec valeurs
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv1 != sv ", sv1 != sv, true);

  // tests vect taille diff
  vector<double> vq1(3);
  for (int i = 0; i < 3; i++)
  {
    vq1[i] = 3;
  }
  sv1.setValues(vq1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.size() == vq1.size() ", sv1.size() == vq1.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv1 ", sv != sv1, true);

  // test vect meme taille 1er nb différent
  vq[0] = 0;
  sv1.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.size() == vq.size() ", sv1.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv1 ", sv != sv1, true);

  // test vect meme taille dernier nb different
  vq[4] = 5;
  vq[1] = 0;
  sv1.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.size() == vq.size() ", sv1.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv1 ", sv != sv1, true);

  // tests vect taille 1
  vector<double> vq2(1);
  vq2[0] = 100.365;
  sv.setValues(vq2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv.isComposite() ", sv.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv.size() == vq2.size() ", sv.size() == vq2.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv1 ", sv != sv1, true);

  sv1.setValues(vq2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.isComposite() ", sv1.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : sv1.size() == vq2.size() ", sv1.size() == vq2.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv != sv1 ", sv != sv1, false);

  SiconosVector *nsv = new SimpleVector(vq2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : nsv->isComposite() ", nsv->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff  : nsv->size() == vq2.size() ", nsv->size() == vq2.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : sv1 != *nsv ", sv1 != *nsv, false);
  delete nsv;

  cout << "SimpleVectorTest >>> testOperatorCompDiff ............................... OK\n ";
}


void SimpleVectorTest::testOperatorEqual()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq);
  SimpleVector vv;
  vv = v;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v.size() == vv.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v == vv, true);

  SimpleVector *pv, *pv1;
  pv = new SimpleVector(vq);
  pv1 =  new SimpleVector();

  *pv1 = *pv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", pv->size() == pv1->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *pv == *pv, true);
  delete pv;
  delete pv1;

  SiconosVector *pv2, *pv3;
  pv2 = new SimpleVector(vq);
  pv3 =  new SimpleVector();

  *pv2 = *pv3;
  pv2->display();
  pv3->display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", pv2->size() == pv3->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *pv2 == *pv3, true);
  delete pv2;
  delete pv3;

  pv = new SimpleVector();
  pv2 =  new SimpleVector(v);

  *pv = *pv2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", pv2->size() == pv->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *pv2 == *pv, true);
  delete pv;
  delete pv2;

  cout << "SimpleVectorTest >>> testOperatorEqual ............................... OK\n ";
}


void SimpleVectorTest::testOperatorAccessRef()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  SimpleVector v(vq);

  v(0) = 10;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", v.isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", v.size() == vq.size() , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", v(0) != vq[0] , true);

  cout << "SimpleVectorTest >>> testOperatorAccessRef ............................... OK\n ";
}

void SimpleVectorTest::testOperatorAccessRefException()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  SimpleVector v(vq);

  v(-1) = 10;
}

void SimpleVectorTest::testOperatorAccessRefException1()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  SimpleVector v(vq);

  v(1000) = 10;
}


void SimpleVectorTest::testOperatorAccessRefException2()
{
  SimpleVector v;
  v(1) = 10;
}

void SimpleVectorTest::testOperatorAccess()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  SimpleVector v(vq);

  const double d = v(3);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", d == v(3), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", d == vq[3] , true);

  cout << "SimpleVectorTest >>> testOperatorAccess ............................... OK\n ";
}

void SimpleVectorTest::testOperatorAccessException()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  SimpleVector v(vq);

  double d = v(-1);
}

void SimpleVectorTest::testOperatorAccessException1()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  SimpleVector v(vq);

  double d =  v(1000);
}


void SimpleVectorTest::testOperatorAccessException2()
{
  SimpleVector v;
  double d =  v(1);
}


//void SimpleVectorTest::testAddException()
//{
//  SimpleVector v;
//  SimpleVector v1;
//  v.add(v1);
//}

void SimpleVectorTest::testDisplay()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  cout << "SimpleVectorTest >>> testDisplay ............................... OK\n ";
}


void SimpleVectorTest::testMiscelleanous()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  cout << "SimpleVectorTest >>> testMiscelleanous ............................... OK\n ";
}


void SimpleVectorTest::testSetValues()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq);

  SimpleVector sv;
  sv.setValues(vq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", v == sv, true);

  SimpleVector v1, v2;
  vector<double> vq2(0);

  v1.setValues(vq2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", v1 == v2, true);

  SimpleVector v3(1);
  v3.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", v3 == sv, true);

  cout << "SimpleVectorTest >>> testSetValues ............................... OK\n ";
}


void SimpleVectorTest::testSize()
{
  SimpleVector v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", v.size() == 0, true);

  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  v.setValues(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", v.size() == vq.size(), true);

  cout << "SimpleVectorTest >>> testBuildSimpleVector ............................... OK\n ";
}

void SimpleVectorTest::testWrite()
{
  SimpleVector v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testWrite : ", v.size() == 0, true);

  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }

  v.setValues(vq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testWrite : ", v.write("testWrite_ascii.dat", "ascii") == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testWrite : ", v.write("testWrite_bin.dat", "binary") == true, true);

  cout << "SimpleVectorTest >>> testWrite ............................... OK\n ";
}


void SimpleVectorTest::testRead()
{
  SimpleVector v, v1;

  v.read("testWrite_ascii.dat", "ascii");
  v1.read("testWrite_bin.dat", "binary");


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v.size() == v1.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", ((v(0) - v1(0) < 0.001) || (v(0) - v1(0) > -0.001)), true);

  cout << "SimpleVectorTest >>> testRead ............................... OK\n ";
}


void SimpleVectorTest::testLapackVSblasMultDouble()
{
  const double d = 17.5;
  const double size = 100000;

  vector<double> vq(size);
  for (int i = 0; i < size; i++)
  {
    vq[i] = d;
  }

  SimpleVector sv(vq), /*sv1(vq),*/ V;
  /*SiconosVector*/
  SimpleVector* sv2 = new SimpleVector(vq);

  clock_t start, end;
  double diff;

  for (int j = 0; j < 10; j++)
  {
    start = clock();

    for (int i = 0; i < 10; i++)
    {
      sv.setValues(vq);
      sv += *sv2;  /* code à timer */
    }
    end = clock();
    diff = diff + (end - start);
  }
  cout << "time = " << (diff / (10.0)) << endl;
  sv.display();

  delete sv2;


  cout << "SimpleVectorTest >>> testLapackVSblasMultDouble ............................... OK\n ";
}

/*******************************************************************************
*         GENERIC INTERNAL OPERATORS                                 *
*******************************************************************************/

void SimpleVectorTest::testOperatorPlusEqualGEN()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq);

  SiconosVector *nsv = new SimpleVector(vq);
  try
  {
    v += *nsv;
  }
  catch (SiconosException e)
  {
    cout << "EXCEPTION testOperatorPlusEqualGEN --- " << e.report() << endl;
    exit(0);
  }
  delete nsv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v(0) == vq[0] + vq[0], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v(4) == vq[4] + vq[4], true);

  cout << "SimpleVectorTest >>> testOperatorPlusEqualGEN ............................... OK\n ";
}

void SimpleVectorTest::testOperatorSubEqualGEN()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq);

  SiconosVector *nsv = new SimpleVector(vq);

  v -= *nsv;
  delete nsv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualGEN : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualGEN : ", v(0) == vq[0] - vq[0], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v(4) == vq[4] - vq[4], true);

  cout << "SimpleVectorTest >>> testOperatorSubEqualGEN ............................... OK\n ";
}

void SimpleVectorTest::testOperatorEqualGEN()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq);

  SiconosVector *nsv = new SimpleVector(vq);

  v = *nsv;
  delete nsv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualGEN : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualGEN : ", v(0) == vq[0], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualGEN : ", v(4) == vq[4], true);

  cout << "SimpleVectorTest >>> testOperatorEqualGEN ............................... OK\n ";
}


/*******************************************************************************
*         SPECIFIC INTERNAL OPERATORS                                *
*******************************************************************************/

void SimpleVectorTest::testOperatorPlusEqualSPC()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq), v1(vq);

  v += v1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualSPC : ", v(0) == vq[0] + vq[0], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualSPC : ", v(4) == vq[4] + vq[4], true);

  cout << "SimpleVectorTest >>> testOperatorPlusEqualSPC ............................... OK\n ";
}


void SimpleVectorTest::testOperatorSubEqualSPC()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1;
  }
  SimpleVector v(vq), v1(vq);

  v -= v1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v(0) == vq[0] - vq[0], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v(4) == vq[4] - vq[4], true);

  cout << "SimpleVectorTest >>> testOperatorSubEqualSPC ............................... OK\n ";
}


void SimpleVectorTest::testOperatorMultEqualSPC()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 1.5;
  }
  SimpleVector v(vq);
  const double d = 2;

  v *= d;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v(0) == vq[0] * d, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v(4) == vq[4] * d, true);

  cout << "SimpleVectorTest >>> testOperatorSubEqualSPC ............................... OK\n ";
}


void SimpleVectorTest::testOperatorDivEqualSPC()
{
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = 3;
  }
  SimpleVector v(vq);
  const double d = 2;
  const double res = vq[0] / d; //cout<<"RES "<<res<<endl;

  v /= d;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", v(0) == res, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", v(4) == res, true);

  SiconosVector *nsv = new SimpleVector(vq);
  *nsv /= d;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", nsv->size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", (*nsv)(0) == res, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", (*nsv)(4) == res, true);

  delete nsv;

  cout << "SimpleVectorTest >>> testOperatorDivEqualSPC ............................... OK\n ";
}


/*******************************************************************************
*         GENERIC EXTERNAL OPERATORS                                 *
*******************************************************************************/

void SimpleVectorTest::testExternalOperatorPlusGEN()
{
  const double d = 1256.36589;
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = d;
  }
  SimpleVector v;

  SiconosVector *nsv = new SimpleVector(vq);
  SiconosVector *nsv1 = new SimpleVector(vq);


  v = *nsv + *nsv1;

  delete nsv;
  delete nsv1;


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v(0) == 2 * d, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v(4) == 2 * d, true);

  cout << "SimpleVectorTest >>> testExternalOperatorPlusGEN ............................... OK\n ";
}

void SimpleVectorTest::testExternalOperatorSubGEN()
{
  const double d = 1256.36589;
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = d;
  }
  SimpleVector v;

  SiconosVector *nsv = new SimpleVector(vq);
  SiconosVector *nsv1 = new SimpleVector(vq);


  v = *nsv - *nsv1;

  delete nsv;
  delete nsv1;


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v(0) == d - d, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v(4) == d - d, true);

  cout << "SimpleVectorTest >>> testExternalOperatorPlusGEN ............................... OK\n ";
}


void SimpleVectorTest::testExternalOperatorMultDoubleSPC()
{
  const double d = 1256.36589;
  const double m = -1546.32658;
  const double res = d * m;
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = d;
  }
  SimpleVector v, v1(vq);

  v = v1 * m;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultDoubleSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultDoubleSPC : ", v(0) == res, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultDoubleSPC : ", v(4) == res, true);

  v.setValues(vq);
  v = m * v1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultDoubleSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultDoubleSPC : ", v(0) == res, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultDoubleSPC : ", v(4) == res, true);

  cout << "SimpleVectorTest >>> testExternalOperatorMultDoubleSPC ............................... OK\n ";
}


void SimpleVectorTest::testExternalOperatorDivDoubleSPC()
{
  const double d = 1256.36589;
  const double m = -1546.32658;
  /*const*/
  double res = d / m;
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = d;
  }
  SimpleVector v, v1(vq);

  v = v1 / m;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorDivDoubleSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorDivDoubleSPC : ", ((v(0) - res < 0.001) || (v(0) - res > -0.001)), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorDivDoubleSPC : ", ((v(4) - res < 0.001) || (v(4) - res > -0.001)), true);

  cout << "SimpleVectorTest >>> testExternalOperatorDivDoubleSPC ............................... OK\n ";
}


void SimpleVectorTest::testExternalOperatorPlusSPC()
{
  const double d = 1256.36589;
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = d;
  }
  SimpleVector v;

  SimpleVector *nsv = new SimpleVector(vq);
  SimpleVector *nsv1 = new SimpleVector(vq);


  v = *nsv + *nsv1;

  delete nsv;
  delete nsv1;


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusSPC : ", v(0) == 2 * d, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusSPC : ", v(4) == 2 * d, true);

  cout << "SimpleVectorTest >>> testExternalOperatorPlusSPC ............................... OK\n ";
}


void SimpleVectorTest::testExternalOperatorSubSPC()
{
  const double d = 1256.36589;
  vector<double> vq(5);
  for (int i = 0; i < 5; i++)
  {
    vq[i] = d;
  }
  SimpleVector v;

  SimpleVector *nsv = new SimpleVector(vq);
  SimpleVector *nsv1 = new SimpleVector(vq);


  v = *nsv - *nsv1;

  delete nsv;
  delete nsv1;


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorSubSPC : ", v.size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorSubSPC : ", v(0) == d - d, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorSubSPC : ", v(4) == d - d, true);

  cout << "SimpleVectorTest >>> testExternalOperatorSubSPC ............................... OK\n ";
}

void SimpleVectorTest::testExternalOperatorMultMat()
{
  SiconosMatrix m(2, 4);
  m(0, 0) = 0;
  m(0, 1) = 1;
  m(0, 2) = -1;
  m(0, 3) = 0;
  m(1, 0) = 2;
  m(1, 1) = 1;
  m(1, 2) = -1;
  m(1, 3) = -2;

  SimpleVector v(4);
  v(0) = 1;
  v(1) = 2;
  v(2) = 3;
  v(3) = 4;

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;


  SimpleVector sv(2);
  sv = m * v;


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);

  cout << "SimpleVectorTest >>> testExternalOperatorMultMat ............................... OK\n ";
}

void SimpleVectorTest::testExternalOperatorMultTransMat()
{
  SiconosMatrix m(4, 2);
  m(0, 0) = 0;
  m(1, 0) = 1;
  m(2, 0) = -1;
  m(3, 0) = 0;
  m(0, 1) = 2;
  m(1, 1) = 1;
  m(2, 1) = -1;
  m(3, 1) = -2;

  SimpleVector v(4);
  v(0) = 1;
  v(1) = 2;
  v(2) = 3;
  v(3) = 4;

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;


  SimpleVector sv(2);
  sv = matTransVecMult(m, v);


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultTransMat : ", sv == res, true);

  cout << "SimpleVectorTest >>> testExternalOperatorMultTransMat ............................... OK\n ";
}

//$Log: SimpleVectorTest.cpp,v $
//Revision 1.12  2005/01/18 11:00:26  jbarbier
//- write wild fixed in simple vector test
//
//Revision 1.11  2005/01/13 09:06:17  charlety
//
//_ documentation is compilable again with : make documentation
//_ xml tests are now available too
//
//Revision 1.10  2004/09/14 13:24:54  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.9  2004/09/09 14:32:50  charlety
//
//_ New tests for operators of multiplication between vectors and matrices.
//
//Revision 1.8  2004/08/23 09:29:01  charlety
//
//_ tests for compositeVector in progress
//
//Revision 1.7  2004/08/20 15:00:35  charlety
//
//_ Tests for operators of SimpleVector
//
//Revision 1.6  2004/08/19 15:21:28  charlety
//
//_ SimpleVector and CompositeVector in progress.
//_ for the operators, we prefer now using directly functions of Blas1++ instead
//  of these of Blas++.h
//
//Revision 1.5  2004/08/16 09:40:01  charlety
//
//_ new tests for the simpleVector
//
//Revision 1.4  2004/08/13 10:36:11  charlety
//
//_ tests for simpleVector in progress
//
//Revision 1.3  2004/08/11 14:16:08  charlety
//
//_ NewSiconosVector in progress...(NewSiconosVector is an abstract class and
//  SimpleVector inherits of NewSiconosVector).
//
//Revision 1.2  2004/07/30 14:21:54  charlety
//
//_ new functions and tests for the new SiconosVector
//
