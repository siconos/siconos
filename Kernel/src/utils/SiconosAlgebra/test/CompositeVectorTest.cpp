
#include "CompositeVectorTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(CompositeVectorTest);


void CompositeVectorTest::setUp()
{
  int i;
  vq.resize(5, 1);   // = [1 1 1 1 1]
  vdotq.resize(3, 1); // = [2 2 2]

  for (i = 0; i < 5; i++)
    vq.at(i) = 1;
  for (i = 0; i < 3; i++)
    vdotq.at(i) = 2;

  simpleVect = new SimpleVector(vdotq);
  q = new SimpleVector(vq);
  r = new SimpleVector(vq);
  CV = new CompositeVector();
  CV->add(*simpleVect);
  CV->add(*r); // CV = [ [2 2 2 ]   [1 1 1 1 1]] tabindex = [3 8]
  SimpleVector * tmp1 = new SimpleVector(3);
  SimpleVector * tmp2 = new SimpleVector(5);
  tmp1->zero();
  tmp2->zero();
  tmp = new CompositeVector();
  tmp->add(*tmp1);
  tmp->add(*tmp2);
  delete tmp2;
  delete tmp1;
}

void CompositeVectorTest::tearDown()
{
  delete tmp;
  delete CV;
  delete r;
  delete q;
  delete simpleVect;
}

// CONSTRUCTORS

// default
void CompositeVectorTest::testBuildCompositeVector()
{
  CompositeVector *v = new CompositeVector();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector : ", v->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector : ", v->size() == 0, true);
  cout << "CompositeVectorTest >>> testBuildCompositeVector ............................... OK\n ";
  delete v;
}

// Copy from SiconosVector and direct copy
void CompositeVectorTest::testBuildCompositeVector1()
{
  // from a SiconosVector which is a simple
  CompositeVector * v = new CompositeVector(*simpleVect);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", v->size() == simpleVect->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", (v->getValues(0)) == (simpleVect->getValues()), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", (v->getTabIndex()).size() == 1 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", (v->getTabIndex())[0] == 3 , true);
  cout << "CompositeVectorTest >>> testBuildCompositeVector1 ............................... OK\n ";
  delete v;

  // copy from composite
  SiconosVector * tmp = new CompositeVector(*CV);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", *tmp == *CV, true);

  // from a SiconosVector which is a composite

  CompositeVector * vv = new CompositeVector(*tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", vv->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", vv->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", vv->getValues(0) == simpleVect->getValues(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", vv->getValues(1) == q->getValues(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", (vv->getTabIndex()).size() == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", (vv->getTabIndex())[0] == 3 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector1 : ", (vv->getTabIndex())[1] == 8 , true);
  cout << "CompositeVectorTest >>> testBuildCompositeVector1 ............................... OK\n ";
  delete vv;
  delete tmp;
}

// display
// getSvref
// getTabIndex
// toString

// ()
void CompositeVectorTest::testOperatorAccess()
{

  (*CV)(2) = 10;
  double d = (*CV)(7);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", CV->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", CV->size() == 8 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", (*(CV->getSvref())[0])(2) == 10 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", d == 1, true);
  cout << "CompositeVectorTest >>> testOperatorAccessRef ............................... OK\n ";
}

// setValue
void CompositeVectorTest::testSetValue()
{
  double a = 4;
  int i = 2;
  CV->setValue(i, a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValue : ", (*CV)(2) == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValue : ", (CV->getValues(0))(2) == 4, true);
  cout << "CompositeVectorTest >>> testSetValue ............................... OK\n ";
}

// getValue
void CompositeVectorTest::testGetValue()
{
  int i = 2;
  CV->setValue(i, 8);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetValue : ", CV->getValue(i) == 8 , true);
  cout << "CompositeVectorTest >>> testGetValue ............................... OK\n ";
}

// setValues
void CompositeVectorTest::testSetValues()
{
  CV->setValues(vq, 0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", CV->size() == 10, true);
  unsigned int size = (CV->getValues(0)).size();
  unsigned int i;
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getValues(0))(i) == 1, true);
  size = (CV->getValues(1)).size();
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getValues(1))(i) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getTabIndex()).size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getTabIndex())[0] == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getTabIndex())[1] == 10, true);
  cout << "CompositeVectorTest >>> testSetValues ............................... OK\n ";
}

// getValues
void CompositeVectorTest::testGetValues()
{
  unsigned int size = (CV->getValues(0)).size();
  unsigned int i;
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getValues(0))(i) == 2, true);
  size = (CV->getValues(1)).size();
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getValues(1))(i) == 1, true);
  cout << "CompositeVectorTest >>> testGetValues ............................... OK\n ";
}

// size
void CompositeVectorTest::testSize()
{
  unsigned int i = CV->size();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", i == 8, true);
  cout << "CompositeVectorTest >>> testSize ............................... OK\n ";
}

// add and addPtr
void CompositeVectorTest::testAdd()
{
  SimpleVector tmp =  *q;
  CV->addPtr(simpleVect);
  CV->add(tmp);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", CV->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", CV->size() == 16, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getTabIndex()).size() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getTabIndex())[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getTabIndex())[1] == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getTabIndex())[2] == 11, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getTabIndex())[3] == 16, true);
  unsigned int size = (CV->getValues(0)).size();
  unsigned int i;
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getValues(0))(i) == 2, true);
  size = (CV->getValues(1)).size();
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getValues(1))(i) == 1, true);
  size = (CV->getValues(2)).size();
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getValues(2))(i) == 2, true);
  size = (CV->getValues(3)).size();
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getValues(3))(i) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getSvref())[2] == simpleVect, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getSvref())[3] != q, true);
  cout << "CompositeVectorTest >>> testAdd ............................... OK\n ";
}

// write, read and read from a file constructor
void CompositeVectorTest::testReadWrite()
{
  // write CV into an ascii file
  bool isok = CV->write("testCompWrite_ascii.dat", "ascii");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWrite : ", isok == true, true);

  // write CV into a binary file
  bool isok2 = CV->write("testCompWrite_bin.dat", "binary");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWrite : ", isok2 == true, true);

  // load v from an ascii file
  string input = "testCompWrite_ascii.dat";
  CompositeVector * v = new CompositeVector(input, true);
  // load v2 from a binary file
  string input_bin = "testCompWrite_bin.dat";
  CompositeVector *v2 = new CompositeVector(input_bin, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size(1) == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size() == v2->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->getTabIndex() == v2->getTabIndex(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (v->getTabIndex())[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (v->getTabIndex())[1] == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (v->getTabIndex()).size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->getValues(0) == simpleVect->getValues() , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->getValues(1) == q->getValues() , true);

  delete v2;
  delete v;
  cout << "CompositeVectorTest >>> testWrite ............................... OK\n ";
}

// OPERATORS

// += -=
void CompositeVectorTest::testOperatorPlusEqual()
{
  SiconosVector *sv = new SimpleVector(*CV);
  *CV += *sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", CV->isComposite(), true);

  *CV -= *sv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 1, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", CV->isComposite(), true);

  delete sv;

  CompositeVector * tmp = new CompositeVector(*CV);

  sv = new CompositeVector(*CV);
  *tmp += *sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", tmp->isComposite(), true);

  *tmp -= *sv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 2, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 1, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", tmp->isComposite(), true);

  delete sv;
  CompositeVector * sv2 = new CompositeVector(*CV);
  *CV += *sv2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", CV->isComposite(), true);

  *CV -= *sv2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 1, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", CV->isComposite(), true);

  delete sv2;
  delete tmp;

  cout << "CompositeVectorTest >>> testOperatorPlusEqualGEN ............................... OK\n ";
}

// =
void CompositeVectorTest::testOperatorEqual()
{
  SiconosVector *v = new SimpleVector(*CV);
  SiconosVector *w = new CompositeVector(*CV);

  *tmp = *v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", tmp->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", tmp->size() == CV->size(), true);
  for (unsigned int i = 0; i < tmp->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == (*CV)(i), true);

  *tmp = *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", tmp->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", tmp->size() == CV->size(), true);
  for (unsigned int i = 0; i < tmp->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == (*CV)(i), true);

  *tmp = *CV;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildCompositeVector3 : ", tmp->isComposite(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", tmp->size() == CV->size(), true);
  for (unsigned int i = 0; i < tmp->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == (*CV)(i), true);

  delete w;
  delete v;
  cout << "CompositeVectorTest >>> testOperatorEqual ............................... OK\n ";
}

// ==, !=

void CompositeVectorTest::testOperatorComp()
{
  SiconosVector *v = new SimpleVector(*CV);
  SiconosVector *w = new CompositeVector(*CV);
  CompositeVector *z = new CompositeVector();
  CompositeVector *z2 = new CompositeVector(*CV);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *CV == *v, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *CV == *w, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *CV == *z2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", !(*CV == *z), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", !(*z == *v), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", !(*z == *w), true);

  delete z2;
  delete z;
  delete w;
  delete v;
  cout << "CompositeVectorTest >>> testOperatorComp ............................... OK\n ";
}


// *= , /=

void CompositeVectorTest::testOperatorMultDivEqual()
{

  double d = 2;

  *CV *= d;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualSPC : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*CV)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*CV)(i) == 2, true);

  *CV /= d;
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*CV)(i) == 2, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*CV)(i) == 1, true);

  cout << "CompositeVectorTest >>> testOperatorMultDivEqual ............................... OK\n ";
}

// addition
void CompositeVectorTest::testAddition()
{
  SiconosVector * sv = new SimpleVector(*CV);

  *tmp = CV->addition(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 2, true);

  delete sv;
  sv = new CompositeVector(*CV);

  *tmp = CV->addition(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 4, true);

  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 2, true);
  delete sv;
  cout << "CompositeVectorTest >>> testAddition ............................... OK\n ";
}

// subtraction
void CompositeVectorTest::testSubtraction()
{
  SiconosVector * sv = new SimpleVector(*CV);

  *tmp = CV->subtraction(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);

  delete sv;
  sv = new CompositeVector(*CV);

  *tmp = CV->subtraction(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);
  delete sv;
  cout << "CompositeVectorTest >>> testSubtraction ............................... OK\n ";
}

// +
void CompositeVectorTest::testExternalOperatorPlusMoins()
{

  *tmp = *CV + *CV;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 2, true);

  *tmp =  *CV - *CV;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);

  cout << "CompositeVectorTest >>> testExternalOperatorPlusMoins ............................... OK\n ";
}

// * /
void CompositeVectorTest::testExternalOperatorMultDiv()
{
  CompositeVector *w = new CompositeVector(*CV);
  double d = 2;

  *w = d * (*CV);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", w->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*w)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*w)(i) == 2, true);

  *w = (*CV) / d;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", w->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*w)(i) == 1, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*w)(i) == 0.5, true);

  delete w;
  cout << "CompositeVectorTest >>> testExternalOperatorMultDiv ............................... OK\n ";
}

// matTransVectMult

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

  CompositeVector *v = new CompositeVector(*simpleVect);
  SimpleVector *v2 = new SimpleVector(1);
  (*v2)(0) = 4;
  (*v)(0) = 1;
  (*v)(1) = 2;
  (*v)(2) = 3;
  v->add(*v2);

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;

  SimpleVector sv(2);
  sv = m * *v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);

  delete v2;
  delete v;
  cout << "CompositeVectorTest >>> testExternalOperatorMultMat ............................... OK\n ";
}

void CompositeVectorTest::testExternalOperatorMatTransMult()
{
  SiconosMatrix m(4, 2);
  m(0, 0) = 0;
  m(0, 1) = 2;
  m(1, 0) = 1;
  m(1, 1) = 1;
  m(2, 0) = -1;
  m(2, 1) = -1;
  m(3, 0) = 0;
  m(3, 1) = -2;

  CompositeVector *v = new CompositeVector(*simpleVect);
  SimpleVector *v2 = new SimpleVector(1);
  (*v2)(0) = 4;
  (*v)(0) = 1;
  (*v)(1) = 2;
  (*v)(2) = 3;
  v->add(*v2);

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;

  SimpleVector sv(2);
  sv = matTransVecMult(m, *v);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMatTransMult : ", sv == res, true);

  delete v2;
  delete v;
  cout << "CompositeVectorTest >>> testExternalOperatorMatTransMult ............................... OK\n ";
}


