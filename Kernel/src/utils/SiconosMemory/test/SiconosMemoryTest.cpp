#include "SiconosMemoryTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryTest);


void SiconosMemoryTest::setUp()
{
  sizeMem = 3;
  unsigned int sizeVect = 3;

  vector<double> v(sizeVect);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  vector<double> w(sizeVect);
  w[0] = 4;
  w[1] = 5;
  w[2] = 6;
  vector<double> z(sizeVect);
  z[0] = 7;
  z[1] = 8;
  z[2] = 9;
  SiconosVector * q1 = new SimpleVector(v);
  SiconosVector * q2 = new SimpleVector(w);
  SiconosVector * q3 = new SimpleVector(z);
  SiconosVector *c1 = new Composite();
  SiconosVector *c2 = new Composite();

  c1->add(*q1);
  c1->add(*q2);
  c2->add(*q3);

  V1.push_back(q1);
  V1.push_back(q2);
  V2.push_back(c1);
  V2.push_back(c2);
}

void SiconosMemoryTest::tearDown()
{
  delete c2;
  delete c1;
  delete q3;
  delete q2;
  delete q1;


}

//______________________________________________________________________________

// copy of a std::vector of siconos vectors
void SiconosMemoryTest::testBuildMemory1()
{

  SiconosMemory * tmp1 = new SiconosMemory(V1)
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->gettmp1emorySize() == sizeMem, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->getVectorMemory().size() == 2, true);

  SiconosMemory * tmp2 = new SiconosMemory(V2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp2->getMemorySize() == sizeMem, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp2->getNbVectorsInMemory() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp2->getVectorMemory().size() == 1, true);


  delete tmp2;
  delete tmp1;

  cout << "SiconosMemoryTest >>> testBuildMemory1 .............................. OK\n ";
}


void SiconosMemoryTest::testBuildMemory2()
{
  SiconosMemory * tmp1 = new SiconosMemory(V1, 4)
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->getVectorMemory().size() == 2, true);

  delete tmp1;

  cout << "SiconosMemoryTest >>> testBuildMemory2 .............................. OK\n ";
}

// copy
void SiconosMemoryTest::testBuildMemory3()
{
  SiconosMemory * tmp1 = new SiconosMemory(V1);

  SiconosMemory *tmp2 = new SiconosMemory(*V1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp2->getMemorySize() == sizeMem, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp2->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp2->getVectorMemory().size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", (tmp2->getVectorMemory())[0] == (tmp1->getVectorMemory())[0] , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", (tmp2->getVectorMemory())[1] == (tmp1->getVectorMemory())[1] , true);

  delete tmp2;
  delete tmp1;

  cout << "SiconosMemoryTest >>> testBuildMemory3 .............................. OK\n ";
}

// setVectorMemory
void SiconosMemoryTest::testSetVectorMemory()
{
  SiconosMemory * tmp1 = new SiconosMemory(V1);

  tmp1->setVectorMemory(V2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == sizeMem, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->getVectorMemory().size() == 1, true);

  delete tmp1;

}

// getSiconosVector
void SiconosMemoryTest::testGetSiconosVector()
{
  vector<SiconosVector*> V1(2);
  V1[0] = V[0];
  V1[1] = V[1];
  SiconosMemory M(M_SIZE, V1);
  SiconosVector *v1;

  v1 = M.getSiconosVector(0);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : *v1 size OK", v1->size() == V[0]->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : v1 values OK", *v1 == *V[0], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : v1 adress OK", v1 != V[0], true);

  cout << "SiconosMemoryTest >>> testGetSiconosVector .......................... OK\n ";
}

void SiconosMemoryTest::testGetSiconosVector1()
{
  vector<SiconosVector*> V1(2);
  V1[0] = V[0];
  V1[1] = V[1];
  SiconosMemory M(M_SIZE, V1);
  SiconosVector *v1;

  v1 = M.getSiconosVector(2); // nbVectorsInMemory < 2 < MemorySize : v1 should be NULL

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector1 : v1 is NULL", v1 == NULL, true);

  cout << "SiconosMemoryTest >>> testGetSiconosVector1 ......................... OK\n ";
}

// getVectorMemory
void SiconosMemoryTest::testGetVectorMemory()
{
  vector<SiconosVector*> VGet;

  SiconosMemory M(V);

  VGet = M.getVectorMemory();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetVectorMemory : VGet size OK", VGet.size() == V.size(), true);
  for (int i = 0; i < M_SIZE; i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetVectorMemory : VGet OK", *VGet[i] == *V[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetVectorMemory : VGet adresses OK", VGet[i] != V[i], true);
  }

  cout << "SiconosMemoryTest >>> testGetVectorMemory ........................... OK\n ";
}

// swap

void SiconosMemoryTest::testSwap()
{
  SiconosMemory M(M_SIZE);

  for (int i = 0; i <= M_SIZE; i++)
  {
    M.swap(V[0]);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(M.getSiconosVector(0)) == *V[0], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : nbVectorsInMemory OK", M.getNbVectorsInMemory() <= M.getMemorySize(), true);
  }

  cout << "SiconosMemoryTest >>> testSwap ...................................... OK\n ";
}


void SiconosMemoryTest::testOperatorEqual()
{
  SiconosMemory M(V);
  SiconosMemory M1;

  M1 = M;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : memorysize OK", M1.getMemorySize() ==  M.getMemorySize(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : nbVectorsInMemory OK", M1.getNbVectorsInMemory() == M.getNbVectorsInMemory(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : size vector OK", M1.getVectorMemory().size() == M.getVectorMemory().size(), true);

  cout << "SiconosMemoryTest >>> testOperatorEqual ...................................... OK\n ";
}

SiconosMemory returnSiconosMemory(void)
{
  vector<SiconosVector*> V(M_SIZE);
  for (int i = 0; i < M_SIZE; i++)
  {
    V[i] = new /*SiconosVector*/SimpleVector(SIZE);
    V[i]->zero();
  }

  SiconosMemory M(V);

  for (int i = 0; i < M_SIZE; i++)
  {
    delete V[i];
  }
  V.clear();

  return M;
}

void SiconosMemoryTest::testOperatorEqual1()
{
  SiconosMemory M(V);
  SiconosMemory M1;

  M1 = returnSiconosMemory();

  M1.display();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual1 : memorysize OK", M1.getMemorySize() ==  M.getMemorySize(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual1 : nbVectorsInMemory OK", M1.getNbVectorsInMemory() == M.getNbVectorsInMemory(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual1 : size vector OK", M1.getVectorMemory().size() == M.getVectorMemory().size(), true);

  cout << "SiconosMemoryTest >>> testOperatorEqual1 ...................................... OK\n ";
}




void SiconosMemoryTest::testMemoryException()
{
  /** SiconosMemory::SiconosMemory(int memorySize, vector<SiconosVector*> V)
   * We test the case where memorySize < V.size
   */

  SiconosMemory M(1, V);
}

void SiconosMemoryTest::testMemoryException1()
{
  /** void SiconosMemory::setVectorMemory(int mSize, vector<SiconosVector*> V)
   * We test the case where memorySize < V.size
   */

  SiconosMemory M;
  M.setVectorMemory(1, V);
}

void SiconosMemoryTest::testMemoryException2()
{
  /** SiconosVector* SiconosMemory::getSiconosVector(int index)
   * We test the case where index > memorySize
   */

  SiconosMemory M(V);

  SiconosVector *v = M.getSiconosVector(V.size() + 10);
}

void SiconosMemoryTest::testMemoryException3()
{
  /** void SiconosMemory::swap(SiconosVector* v)
   * We test the case where we try to put a vector in the memory when memorySize <= 0
   */

  SiconosMemory M;

  M.swap(V[0]);
}


