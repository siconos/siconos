/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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

#include "BlockVectorTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(BlockVectorTest);


void BlockVectorTest::setUp()
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
  CV = new BlockVector(*simpleVect);
  //  CV->add(*simpleVect);
  CV->add(*r); // CV = [ [2 2 2 ]   [1 1 1 1 1]] tabindex = [3 8]
  SimpleVector * tmp1 = new SimpleVector(3);
  SimpleVector * tmp2 = new SimpleVector(5);
  tmp1->zero();
  tmp2->zero();
  tmp = new BlockVector();
  tmp->add(*tmp1);
  tmp->add(*tmp2);
  delete tmp2;
  delete tmp1;
}

void BlockVectorTest::tearDown()
{
  delete tmp;
  delete CV;
  delete r;
  delete q;
  delete simpleVect;
}

// CONSTRUCTORS

// default
void BlockVectorTest::testBuildBlockVector()
{
  cout << "===================================" << endl;
  cout << "=== Block Vector tests start ...=== " << endl;
  cout << "===================================" << endl;
  cout << "--> Test: buildBlockVector." << endl;
  BlockVector *v = new BlockVector();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector : ", v->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector : ", v->size() == 0, true);
  delete v;
  cout << "--> buildBlockVector test ended with success." << endl;
}

// Copy from SiconosVector and direct copy
void BlockVectorTest::testBuildBlockVector1()
{
  cout << "--> Test: constructor 1." << endl;
  // copy from a block vector
  BlockVector * tmp = new BlockVector(*simpleVect);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", (*tmp) == (*simpleVect), true);

  // from a SiconosVector which is a block vector

  BlockVector * vv = new BlockVector(*CV);

  int a = (vv->getTabIndex())[0];

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", vv->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", vv->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", ((*(vv->getVectorPtr(0)))) == (*simpleVect), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", ((*(vv->getVectorPtr(1)))) == (*q) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", (vv->getTabIndex()).size() == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", a == 3 , true);
  a = (vv->getTabIndex())[1];
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector1 : ", a == 8 , true);
  delete vv;
  delete tmp;
  cout << "--> Constructor 1 test ended with success." << endl;
}

// display
// getSvref
// getTabIndex
// toString

// ()
void BlockVectorTest::testOperatorAccess()
{
  cout << "--> Test: operatorAccess ." << endl;
  (*CV)(2) = 10;
  double d = (*CV)(7);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", CV->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccess : ", CV->size() == 8 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", (*(CV->getSvref())[0])(2) == 10 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", d == 1, true);
  cout << "--> operatorAccess test ended with success." << endl;
}

// setValue
void BlockVectorTest::testSetValue()
{
  cout << "--> Test: setValue." << endl;
  double a = 4;
  int i = 2;
  CV->setValue(i, a);
  double b = (*CV)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValue : ", b == 4, true);
  b = (CV->getValues(0))(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValue : ", b == 4, true);
  cout << "-->  setValue test ended with success." << endl;
}

// getValue
void BlockVectorTest::testGetValue()
{
  cout << "--> Test: getValue." << endl;
  int i = 2;
  CV->setValue(i, 8);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetValue : ", CV->getValue(i) == 8 , true);
  cout << "-->  getValue test ended with success." << endl;
}

// setValues
void BlockVectorTest::testSetValues()
{
  cout << "--> Test: setValues." << endl;
  CV->setValues(vq, 0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", CV->size() == 10, true);
  unsigned int size = (CV->getValues(0)).size();
  unsigned int i;
  LaVectorDouble tmp = CV->getValues(0);
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", tmp(i) == 1, true);
  tmp = CV->getValues(1);
  size = tmp.size();
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", tmp(i) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (CV->getTabIndex()).size() == 2, true);
  double a = (CV->getTabIndex())[0];
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", a == 5, true);
  a = (CV->getTabIndex())[1];
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", a == 10, true);
  cout << "-->  setValues test ended with success." << endl;
}

// getValues
void BlockVectorTest::testGetValues()
{
  cout << "--> Test: getValues." << endl;
  unsigned int size = (CV->getValues(0)).size();
  unsigned int i;
  LaVectorDouble tmp = CV->getValues(0);
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", tmp(i) == 2, true);
  tmp = CV->getValues(1);
  size = tmp.size();
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", tmp(i) == 1, true);
  cout << "-->  getValues test ended with success." << endl;
}

// size
void BlockVectorTest::testSize()
{
  cout << "--> Test: setSize." << endl;
  unsigned int i = CV->size();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", i == 8, true);
  cout << "-->  setSize test ended with success." << endl;
}

// add and addPtr
void BlockVectorTest::testAdd()
{
  cout << "--> Test: add." << endl;
  SimpleVector tmp =  *q;
  CV->addPtr(simpleVect);
  CV->add(tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", CV->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", CV->size() == 16, true);
  double a = (CV->getTabIndex())[0];
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getTabIndex()).size() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", a == 3, true);
  a = (CV->getTabIndex())[1];
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", a == 8, true);
  a = (CV->getTabIndex())[2];
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", a == 11, true);
  a = (CV->getTabIndex())[3];
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", a == 16, true);
  unsigned int size = (CV->getValues(0)).size();
  unsigned int i;
  LaVectorDouble tmp2 = CV->getValues(0);

  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", tmp2(i) == 2, true);
  size = (CV->getValues(1)).size();
  tmp2 = CV->getValues(1);
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", tmp2(i) == 1, true);
  size = (CV->getValues(2)).size();
  tmp2 = CV->getValues(2);
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", tmp2(i) == 2, true);
  size = (CV->getValues(3)).size();
  tmp2 = CV->getValues(3);
  for (i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", tmp2(i) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getSvref())[2] == simpleVect, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAdd : ", (CV->getSvref())[3] != q, true);
  cout << "-->  add test ended with success." << endl;
}

// write, read and read from a file constructor
void BlockVectorTest::testReadWrite()
{
  cout << "--> Test: readWrite." << endl;
  // write CV into an ascii file
  bool isok = CV->write("testCompWrite_ascii.dat", "ascii");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWrite : ", isok == true, true);

  // write CV into a binary file
  bool isok2 = CV->write("testCompWrite_bin.dat", "binary");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWrite : ", isok2 == true, true);

  // load v from an ascii file
  string input = "testCompWrite_ascii.dat";
  BlockVector * v = new BlockVector(input, true);
  // load v2 from a binary file
  string input_bin = "testCompWrite_bin.dat";
  BlockVector *v2 = new BlockVector(input_bin, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size(1) == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size() == v2->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->getTabIndex() == v2->getTabIndex(), true);
  double a = (v->getTabIndex())[0] ;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", a == 3, true);
  a = (v->getTabIndex())[1] ;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", a == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (v->getTabIndex()).size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (*(v->getVectorPtr(0))) == (*simpleVect) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (*(v->getVectorPtr(1))) == (*q) , true);

  delete v2;
  delete v;
  cout << "--> readWrite test ended with success." << endl;
}

// OPERATORS

// += -=
void BlockVectorTest::testOperatorPlusEqual()
{
  cout << "--> Test: operatorPlusEqual." << endl;
  SiconosVector *sv = new SimpleVector(*CV);
  *CV += *sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", CV->isBlock(), true);

  *CV -= *sv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 1, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", CV->isBlock(), true);

  delete sv;

  BlockVector * tmp = new BlockVector(*CV);

  sv = new BlockVector(*CV);
  *tmp += *sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", tmp->isBlock(), true);

  *tmp -= *sv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 2, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == 1, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", tmp->isBlock(), true);

  delete sv;
  BlockVector * sv2 = new BlockVector(*CV);
  *CV += *sv2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", CV->isBlock(), true);

  *CV -= *sv2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", CV->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 2, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*CV)(i) == 1, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", CV->isBlock(), true);

  delete sv2;
  delete tmp;

  cout << "--> operatorPlusEqual test ended with success." << endl;
}

// =
void BlockVectorTest::testOperatorEqual()
{
  cout << "--> Test: operatorEqual." << endl;
  SiconosVector *v = new SimpleVector(*CV);
  SiconosVector *w = new BlockVector(*CV);

  *tmp = *v;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", tmp->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", tmp->size() == CV->size(), true);
  for (unsigned int i = 0; i < tmp->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == (*CV)(i), true);

  *tmp = *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", tmp->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", tmp->size() == CV->size(), true);
  for (unsigned int i = 0; i < tmp->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == (*CV)(i), true);

  *tmp = *CV;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildBlockVector3 : ", tmp->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", tmp->size() == CV->size(), true);
  for (unsigned int i = 0; i < tmp->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*tmp)(i) == (*CV)(i), true);

  delete w;
  delete v;
  cout << "--> operatorEqual test ended with success." << endl;
}


// *= , /=

void BlockVectorTest::testOperatorMultDivEqual()
{

  cout << "--> Test: operatorMultDivEqual." << endl;
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

  cout << "--> operatorMultDivEqual test ended with success." << endl;
}

// addition
void BlockVectorTest::testAddition()
{
  cout << "--> Test: addition." << endl;
  SiconosVector * sv = new SimpleVector(*CV);

  *tmp = CV->addition(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 4, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 2, true);

  delete sv;
  sv = new BlockVector(*CV);

  *tmp = CV->addition(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 4, true);

  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 2, true);
  delete sv;
  cout << "--> addition test ended with success." << endl;
}

// subtraction
void BlockVectorTest::testSubtraction()
{
  cout << "--> Test: subtraction ." << endl;
  SiconosVector * sv = new SimpleVector(*CV);

  *tmp = CV->subtraction(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);

  delete sv;
  sv = new BlockVector(*CV);

  *tmp = CV->subtraction(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", tmp->size() == 8, true);
  for (unsigned int i = 0; i < 3; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);
  for (unsigned int i = 3; i < 8; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*tmp)(i) == 0, true);
  delete sv;
  cout << "--> subtraction test ended with success." << endl;
}

// +
void BlockVectorTest::testExternalOperatorPlusMoins()
{
  cout << "--> Test: externalOperatorPlusMoins ." << endl;

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

  cout << "--> externalOperatorPlusMoins test ended with success." << endl;
}

// * /
void BlockVectorTest::testExternalOperatorMultDiv()
{
  cout << "--> Test: externalOperatorMultDiv ." << endl;
  BlockVector *w = new BlockVector(*CV);
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
  cout << "--> externalOperatorMultDiv test ended with success." << endl;
}

// matTransVectMult

void BlockVectorTest::testExternalOperatorMultMat()
{
  cout << "--> Test: externalOperatorMultMat ." << endl;
  SimpleMatrix m(2, 4);
  m(0, 0) = 0;
  m(0, 1) = 1;
  m(0, 2) = -1;
  m(0, 3) = 0;
  m(1, 0) = 2;
  m(1, 1) = 1;
  m(1, 2) = -1;
  m(1, 3) = -2;

  BlockVector *v = new BlockVector(*simpleVect);
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
  cout << "--> externalOperatorMultMat test ended with success." << endl;
}

void BlockVectorTest::testExternalOperatorMatTransMult()
{
  cout << "--> Test: externalOperatorMatTransMult ." << endl;
  SimpleMatrix m(4, 2);
  m(0, 0) = 0;
  m(0, 1) = 2;
  m(1, 0) = 1;
  m(1, 1) = 1;
  m(2, 0) = -1;
  m(2, 1) = -1;
  m(3, 0) = 0;
  m(3, 1) = -2;

  BlockVector *v = new BlockVector(*simpleVect);
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
  cout << "--> externalOperatorMatTransMult test ended with success." << endl;
}

void BlockVectorTest::End()
{
  cout << "=====================================" << endl;
  cout << " ===== End of BlockVector Tests ===== " << endl;
  cout << "=====================================" << endl;
}



