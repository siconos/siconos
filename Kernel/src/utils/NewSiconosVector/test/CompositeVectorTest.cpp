//$Id: CompositeVectorTest.cpp,v 1.7 2004/09/14 13:24:54 charlety Exp $

#include "CompositeVectorTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(CompositeVectorTest);


void CompositeVectorTest::setUp()
{
  int i;
  vector<double> vq(5);
  vector<double> vdotq(5);

  for (i = 0; i < 5; i++)
  {
    vq.at(i) = 1;
    vdotq.at(i) = 2;
  }

  q.setValues(vq);
  dotq.setValues(vdotq);
}

void CompositeVectorTest::tearDown()
{ }

//______________________________________________________________________________

void CompositeVectorTest::testBuildCompositeVector()
{
  CompositeVector v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector : ", v.size() == 0, true);

  cout << "CompositeVectorTest >>> testBuildCompositeVector ............................... OK\n ";
}

void CompositeVectorTest::testBuildCompositeVector1()
{
  CompositeVector v;
  v.add(q);
  v.add(dotq);
  //  v.display();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v.size() == q.size() + dotq.size(), true);

  CompositeVector v1;
  v1.add(q);
  v1.add(v);
  v1.add(dotq);
  //  v1.display();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v1.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v1.size() == q.size() + dotq.size() + v.size(), true);

  CompositeVector v2(v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v2.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v2 == v, true);

  cout << "CompositeVectorTest >>> testBuildCompositeVector1 ............................... OK\n ";
}


void CompositeVectorTest::testOperatorAccess()
{
  CompositeVector v;
  v.add(q);
  v.add(dotq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", v.size() == q.size() + dotq.size(), true);

  for (int i = 0; i < q.size() + dotq.size(); i++)
  {
    if (i < q.size())
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", v(i) == q(i) , true);
    else
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", v(i) == dotq(i - q.size()) , true);
  }
  cout << "CompositeVectorTest >>> testOperatorAccess ............................... OK\n ";
}


void CompositeVectorTest::testAdd()
{
  CompositeVector v;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", v.size() == 0, true);

  v.add(q);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", v.size() == q.size(), true);

  cout << "CompositeVectorTest >>> testAdd ............................... OK\n ";
}

void CompositeVectorTest::testSize()
{
  CompositeVector v;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", v.size() == 0, true);

  v.add(q);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", v.size() == q.size(), true);

  v.add(dotq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", v.size() == dotq.size() + q.size(), true);

  cout << "CompositeVectorTest >>> testSize ............................... OK\n ";
}


void CompositeVectorTest::testWrite()
{
  CompositeVector v;
  v.add(q);
  v.add(dotq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testWrite : ", v.write("testWriteComposite_ascii.dat", "ascii") == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testWrite : ", v.write("testWriteComposite_bin.dat", "binary") == true, true);

  cout << "CompositeVectorTest >>> testWrite ............................... OK\n ";
}


void CompositeVectorTest::testRead()
{
  SimpleVector sv1(q.size()), sv2(dotq.size()), sv3(q.size()), sv4(dotq.size());
  CompositeVector v, v1;

  v.add(sv1);
  v.add(sv2);
  v1.add(sv3);
  v1.add(sv4);

  v.read("testWriteComposite_ascii.dat", "ascii");
  v1.read("testWriteComposite_bin.dat", "binary");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v.size() == v1.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", ((v(0) - v1(0) < 0.001) || (v(0) - v1(0) > -0.001)), true);

  cout << "CompositeVectorTest >>> testRead ............................... OK\n ";
}


/*******************************************************************************
*         GENERIC INTERNAL OPERATORS                                 *
*******************************************************************************/

void CompositeVectorTest::testOperatorPlusEqualGEN()
{
  vector<double> vq(50);
  int i;
  for (i = 0; i < vq.size(); i++)
  {
    vq.at(i) = 37.5;
  }

  SimpleVector sv(vq), sv1(15), sv2(35);

  for (i = 0; i < sv1.size(); i++)
  {
    sv1(i) = 12.5;
  }

  for (i = 0; i < sv2.size(); i++)
  {
    sv2(i) = 12.5;
  }


  CompositeVector v;
  v.add(sv1);
  v.add(sv2);

  v += sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v.size() == sv.size(), true);
  for (i = 0; i < v.size(); i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v(i) == 50, true);
  }

  cout << "CompositeVectorTest >>> testOperatorPlusEqualGEN ............................... OK\n ";
}


void CompositeVectorTest::testOperatorSubEqualGEN()
{
  vector<double> vq(50);
  int i;
  for (i = 0; i < vq.size(); i++)
  {
    vq.at(i) = 37.5;
  }

  SimpleVector sv(vq), sv1(15), sv2(35);

  for (i = 0; i < sv1.size(); i++)
  {
    sv1(i) = 12.5;
  }

  for (i = 0; i < sv2.size(); i++)
  {
    sv2(i) = 12.5;
  }


  CompositeVector v;
  v.add(sv1);
  v.add(sv2);

  v -= sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualGEN : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualGEN : ", v.size() == sv.size(), true);
  for (i = 0; i < v.size(); i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualGEN : ", v(i) == -25.0, true);
  }

  cout << "CompositeVectorTest >>> testOperatorSubEqualGEN ............................... OK\n ";
}

void CompositeVectorTest::testOperatorEqualGEN()
{
  vector<double> vq(50);
  int i;
  for (i = 0; i < vq.size(); i++)
  {
    vq.at(i) = 37.5;
  }

  SimpleVector sv(vq), sv1(15), sv2(35);

  for (i = 0; i < sv1.size(); i++)
  {
    sv1(i) = 12.5;
  }

  for (i = 0; i < sv2.size(); i++)
  {
    sv2(i) = 12.5;
  }


  CompositeVector v;
  v.add(sv1);
  v.add(sv2);

  v = sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualGEN : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualGEN : ", v.size() == sv.size(), true);
  for (i = 0; i < v.size(); i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualGEN : ", v(i) == 37.5, true);
  }

  cout << "CompositeVectorTest >>> testOperatorEqualGEN ............................... OK\n ";
}


void CompositeVectorTest::testOperatorComp()
{
  CompositeVector v, v1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : ", v == v1, true);

  v.add(q);
  v1.add(q);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : ", v == v1, true);

  v.add(dotq);
  v1.add(dotq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : ", v == v1, true);

  v.add(q);
  v1.add(dotq);
  //  cout<<"V = "<<endl;
  //  v.display();
  //  cout<<"V1 = "<<endl;
  //  v1.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : ", v == v1, false);

  CompositeVector v2(v);
  //  v2.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorComp : ", v2 == v, true);

  cout << "CompositeVectorTest >>> testOperatorComp ............................... OK\n ";
}

void CompositeVectorTest::testOperatorCompDiff()
{
  CompositeVector v, v1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : ", v != v1, false);

  v.add(q);
  v1.add(q);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : ", v != v1, false);

  v.add(dotq);
  v1.add(dotq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : ", v != v1, false);

  v.add(q);
  v1.add(dotq);
  //  cout<<"V = "<<endl;
  //  v.display();
  //  cout<<"V1 = "<<endl;
  //  v1.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : ", v != v1, true);

  CompositeVector v2(v);
  //  v2.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorCompDiff : ", v2 != v, false);

  cout << "CompositeVectorTest >>> testOperatorCompDiff ............................... OK\n ";
}

/*******************************************************************************
*         SPECIFIC INTERNAL OPERATORS                                *
*******************************************************************************/

void CompositeVectorTest::testOperatorPlusEqualSPC()
{
  try
  {
    CompositeVector v, v1;
    SimpleVector sv(q), sv1(dotq);
    const double res = q(0) + dotq(0);

    v.add(q);
    v.add(sv);
    v1.add(dotq);
    v1.add(sv1);

    cout << "ADDED" << endl;

    v += v1;

    cout << "OPERATED" << endl;
    //v.display();
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualSPC : ", v.isComposite(), true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualSPC : ", v.size() == 2 * q.size(), true);
    for (int i = 0; i < v.size(); i++)
    {
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualGEN : ", v(i) == res, true);
    }
  }
  catch (SiconosException e)
  {
    cout << "EXCEPTION testOperatorPlusEqualSPC --- " << e.report() << endl;
    exit(0);
  }
  cout << "CompositeVectorTest >>> testOperatorPlusEqualSPC ............................... OK\n ";
}


void CompositeVectorTest::testOperatorSubEqualSPC()
{
  CompositeVector v, v1;
  SimpleVector sv(q), sv1(dotq);
  const double res = q(0) - dotq(0);
  v.add(q);
  v.add(sv);
  v1.add(dotq);
  v1.add(sv1);

  v -= v1;
  //v.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v.size() == 2 * q.size(), true);
  for (int i = 0; i < v.size(); i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorSubEqualSPC : ", v(i) == res, true);
  }

  cout << "CompositeVectorTest >>> testOperatorSubEqualSPC ............................... OK\n ";
}


void CompositeVectorTest::testOperatorMultEqualSPC()
{
  CompositeVector v;
  const double m = 369.3265;
  SimpleVector sv;

  v.add(q);
  v.add(dotq);

  sv = v;

  v *= m;
  //v.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualSPC : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualSPC : ", v.size() == dotq.size() + q.size(), true);
  for (int i = 0; i < v.size(); i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualSPC : ", v(i) == sv(i) * m, true);
  }

  cout << "CompositeVectorTest >>> testOperatorMultEqualSPC ............................... OK\n ";
}


void CompositeVectorTest::testOperatorDivEqualSPC()
{
  CompositeVector v;
  const double m = 369.3265;
  SimpleVector sv;

  v.add(q);
  v.add(dotq);

  sv = v;
  v /= m;
  //v.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", v.size() == dotq.size() + q.size(), true);
  for (int i = 0; i < v.size(); i++)
  {

    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", ((v(i) - sv(i) / m < 0.001) || (v(i) - sv(i) / m > -0.001)), true);
  }

  cout << "CompositeVectorTest >>> testOperatorDivEqualSPC ............................... OK\n ";
}


void CompositeVectorTest::testOperatorEqualSPC()
{
  CompositeVector v, v1, v2, v3;
  SimpleVector sv(5), sv1(5), sv2(5), sv3(10);

  v1 = v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualSPC : ", v.isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualSPC : ", v.size() == 0, true);

  v.add(q);
  v.add(dotq);
  v1.add(sv);
  v1.add(sv1);

  v1 = v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualSPC : ", v == v1, true);

  v2.add(q);
  v2.add(dotq);
  v2.add(q);
  v3.add(sv);
  v3.add(sv3);

  v3 = v2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqualSPC : ", v3 == v2, true);


  cout << "CompositeVectorTest >>> testOperatorEqualSPC ............................... OK\n ";
}


/*******************************************************************************
*         GENERIC EXTERNAL OPERATORS                                 *
/******************************************************************************/

void CompositeVectorTest::testExternalOperatorPlusGEN()
{
  SimpleVector sv(5), sv1(5); // to compose a composite
  for (int i = 0; i < 5; i++)
  {
    sv(i) = 0.0;
    sv1(i) = 0.0;
  }
  CompositeVector v;
  v.add(sv);
  v.add(sv1);

  CompositeVector v1;
  v1.add(q);
  v1.add(dotq);

  SimpleVector sv3(10);
  for (int i = 0; i < 10; i++)
  {
    sv3(i) = i + 3.250154;
  }
  //cout<<" v1.display "; v1.display();
  v = v1.addition(sv3);
  //cout<<" v1.display "; v1.display(); //getchar();


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v.size() == 10, true);
  for (int i = 0; i < 10; i++)
  {
    //cout<<v1(i) + sv3(i)<<" ";
    //CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : values", v(i) == v1(i) + sv3(i), true);
  }
  cout << endl;
  cout << "CompositeVectorTest >>> testExternalOperatorPlusGEN ............................... OK\n ";
}


void CompositeVectorTest::testExternalOperatorSubGEN()
{

}


void CompositeVectorTest::testExternalOperatorMultDoubleSPC()
{

}


void CompositeVectorTest::testExternalOperatorDivDoubleSPC()
{

}


void CompositeVectorTest::testExternalOperatorPlusSPC()
{

}


void CompositeVectorTest::testExternalOperatorSubSPC()
{

}


void CompositeVectorTest::testExternalOperatorMultMat()
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

  SimpleVector v1(2), v2(2);

  v1(0) = 1;
  v1(1) = 2;
  v2(0) = 3;
  v2(1) = 4;

  CompositeVector cv;
  cv.add(v1);
  cv.add(v2);

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;


  SimpleVector sv(2);
  sv = m * cv;


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);

  cout << "CompositeVectorTest >>> testExternalOperatorMultMat ............................... OK\n ";
}

void CompositeVectorTest::testExternalOperatorMultTransMat()
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

  SimpleVector v1(2), v2(2);

  v1(0) = 1;
  v1(1) = 2;
  v2(0) = 3;
  v2(1) = 4;

  CompositeVector cv;
  cv.add(v1);
  cv.add(v2);

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;


  SimpleVector sv(2);
  sv = matTransVecMult(m, cv);


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultTransMat : ", sv == res, true);

  cout << "CompositeVectorTest >>> testExternalOperatorMultTransMat ............................... OK\n ";
}




//$Log: CompositeVectorTest.cpp,v $
//Revision 1.7  2004/09/14 13:24:54  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.6  2004/09/09 14:32:50  charlety
//
//_ New tests for operators of multiplication between vectors and matrices.
//
//Revision 1.5  2004/08/24 11:29:20  charlety
//
//_ methods replacing generic operators for mixed operations done.
//
//Revision 1.4  2004/08/24 07:35:09  charlety
//
//_ CompositeVector finished at 95%.
//
//TODO :
//_ solve the problem of operators (gcc cannot choice between theoperators of simpleVector and compositeVector in mixed operations).
//_ the behavior of specific operator= of CompositeVector is not ok. Should copy values only.
//
//Revision 1.3  2004/08/23 09:29:01  charlety
//
//_ tests for compositeVector in progress
//
//Revision 1.2  2004/08/19 15:21:28  charlety
//
//_ SimpleVector and CompositeVector in progress.
//_ for the operators, we prefer now using directly functions of Blas1++ instead
//  of these of Blas++.h
//
//Revision 1.1  2004/08/17 15:01:46  charlety
//
//_ composite Vector in progress.
//_ Tests for compositeVector created
//