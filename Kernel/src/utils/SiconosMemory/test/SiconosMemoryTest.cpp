//$Id: SiconosMemoryTest.cpp,v 1.4 2004/09/10 11:26:25 charlety Exp $
#include "SiconosMemoryTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryTest);


void SiconosMemoryTest::setUp()
{
  //  A = new SiconosVector(SIZE);
  //  B = new SiconosVector(SIZE);
  //  C = new SiconosVector(SIZE);
  //  A->zero();
  //  B->zero();
  //  C->zero();

  V.resize(M_SIZE);
  for (int i = 0; i < M_SIZE; i++)
  {
    V[i] = new /*SiconosVector*/SimpleVector(SIZE);
    V[i]->zero();
  }
}

void SiconosMemoryTest::tearDown()
{
  //  delete A;
  //  delete B;
  //  delete C;

  for (int i = 0; i < M_SIZE; i++)
  {
    delete V[i];
  }
}

//______________________________________________________________________________


void SiconosMemoryTest::testBuildMemory()
{
  SiconosMemory M;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory : memorysize OK", M.getMemorySize() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory : nbVectorsInMemory OK", M.getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory : size vector OK", M.getVectorMemory().size() == 0, true);

  cout << "SiconosMemoryTest >>> testBuildMemory ............................... OK\n ";
}


void SiconosMemoryTest::testBuildMemory1()
{
  SiconosMemory M(M_SIZE);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", M.getMemorySize() == M_SIZE, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", M.getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", M.getVectorMemory().size() == M_SIZE, true);

  cout << "SiconosMemoryTest >>> testBuildMemory1 .............................. OK\n ";
}


void SiconosMemoryTest::testBuildMemory2()
{
  SiconosMemory M(V);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : memorysize OK", M.getMemorySize() == M_SIZE, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : nbVectorsInMemory OK", M.getNbVectorsInMemory() == M_SIZE, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : size vector OK", M.getVectorMemory().size() == M_SIZE, true);

  cout << "SiconosMemoryTest >>> testBuildMemory2 .............................. OK\n ";
}

void SiconosMemoryTest::testBuildMemory3()
{
  vector<SiconosVector*> V1(2);
  V1[0] = V[0];
  V1[1] = V[1];
  SiconosMemory M(M_SIZE, V1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : memorysize OK", M.getMemorySize() == M_SIZE, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : nbVectorsInMemory OK", M.getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : size vector OK", M.getVectorMemory().size() == M_SIZE, true);

  cout << "SiconosMemoryTest >>> testBuildMemory3 .............................. OK\n ";
}

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

void SiconosMemoryTest::testSetVectorMemory()
{
  SiconosMemory M;

  M.setVectorMemory(V);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory : memorysize OK", M.getMemorySize() == V.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory : nbVectorsInMemory OK", M.getNbVectorsInMemory() == V.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory : size vector OK", M.getVectorMemory().size() == V.size(), true);

  cout << "SiconosMemoryTest >>> testSetVectorMemory ........................... OK\n ";
}

void SiconosMemoryTest::testSetVectorMemory1()
{
  vector<SiconosVector*> V1(2);
  V1[0] = V[0];
  V1[1] = V[1];
  SiconosMemory M;

  M.setVectorMemory(M_SIZE, V1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory1 : memorysize OK", M.getMemorySize() == M_SIZE, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory1 : nbVectorsInMemory OK", M.getNbVectorsInMemory() == V1.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory1 : size vector OK", M.getVectorMemory().size() == M_SIZE, true);

  for (int i = 0; i < 2; i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory1 : vector OK", *(M.getSiconosVector(i)) == *V1[i], true);
  }

  cout << "SiconosMemoryTest >>> testSetVectorMemory1 .......................... OK\n ";
}

void SiconosMemoryTest::testSetVectorMemory2()
{
  SiconosMemory M(M_SIZE + 10);

  M.setVectorMemory(M_SIZE, V); // try to set a memory already initialized. Should be ok.

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory2 : memorysize OK", M.getMemorySize() == M_SIZE, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory2 : nbVectorsInMemory OK", M.getNbVectorsInMemory() == V.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory2 : size vector OK", M.getVectorMemory().size() == M_SIZE, true);

  for (int i = 0; i < V.size(); i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetVectorMemory2 : vector OK", *(M.getSiconosVector(i)) == *V[i], true);
  }

  cout << "SiconosMemoryTest >>> testSetVectorMemory2 .......................... OK\n ";
}

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


//$Log: SiconosMemoryTest.cpp,v $
//Revision 1.4  2004/09/10 11:26:25  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.3  2004/07/27 14:41:15  charlety
//
//_ new tests for operator = of SiconosMemory
//_ modification of the makefile : the variable LIBS was false since the object
//  SiconosMemoryXML is linked with SiconosMemory.
//
//Revision 1.2  2004/07/08 09:13:43  charlety
//
//_ creation of a Siconos exception dedicated to the class SiconosMemory
//_ new tests for SiconosMemory
//
//Revision 1.1  2004/07/07 13:53:13  charlety
//
//_ First version of the memory object
//_ some little bugs corrected otherwhere
//