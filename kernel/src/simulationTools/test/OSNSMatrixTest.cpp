/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "OSNSMatrixTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(OSNSMatrixTest);


void OSNSMatrixTest::setUp()
{
  n = 5;
  tol = 1e-12;
  // Download a Model from Template.xml file
  temp.reset(new Model("Template.xml"));
  SP::TimeStepping s = std11::static_pointer_cast<TimeStepping>(temp->simulation());
  s->initialize();
  // Get a set of Interactions
  indexSet = s->indexSet(0);
  SP::OneStepNSProblem osns = s->getOneStepNSProblems()->begin()->second;
  osns->computeAllBlocks();
  blocks = osns->getBlocks();
}

void OSNSMatrixTest::tearDown()
{}

// Default constructor
void OSNSMatrixTest::testBuildOSNSMatrix0()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== OSNSMatrix tests start ...===== " <<std::endl;
  std::cout << "===========================================" <<std::endl;
  std::cout << "------- Default constructor test -------" <<std::endl;
  SP::OSNSMatrix  M = new OSNSMatrix();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix0 : ", M->size() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix0 : ", M->storagetype() == 0, true);
  std::cout << "------- Default constructor OSNSMatrix ok -------" <<std::endl;
  std::cout <<std::endl <<std::endl;
}

// Basic constructor
void OSNSMatrixTest::testBuildOSNSMatrix1()
{
  std::cout << "------- Constructor with dim. test -------" <<std::endl;
  SP::OSNSMatrix  M(new OSNSMatrix(n));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->size() == n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->storagetype() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->defaultMatrix(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->defaultMatrix()->size(0) == n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->defaultMatrix()->size(1) == n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->defaultMatrix()->normInf() < tol , true);
  std::cout << "------- Constructor(dim) ended with success -------" <<std::endl;
  std::cout <<std::endl <<std::endl;
}

void OSNSMatrixTest::testBuildOSNSMatrix2()
{
  std::cout << "------- Constructor with index set and list of blocks -------" <<std::endl;

  SP::OSNSMatrix  M(new OSNSMatrix(indexSet, blocks));

  unsigned int dim = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->nonSmoothLaw()->size();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (InteractionsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (InteractionsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->nonSmoothLaw()->size();
    }
    col = 0;
    row += (*itRow)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->storagetype() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->defaultMatrix(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->defaultMatrix()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->defaultMatrix()->size(1) == dim, true);
  unsigned int i = 0, pos = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getSizeOfDiagonalBlock(i++) == (*it)->nonSmoothLaw()->size(), true);
    pos += (*it)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", (*M->defaultMatrix() - MRef).normInf() < tol, true);

  std::cout << "------- Constructor(indexSet,blocks) ended with success -------" <<std::endl;
  std::cout <<std::endl <<std::endl;
}

void OSNSMatrixTest::testFill()
{
  std::cout << "------- fill function test -------" <<std::endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  unsigned int dim = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->nonSmoothLaw()->size();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (InteractionsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (InteractionsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->nonSmoothLaw()->size();
    }
    col = 0;
    row += (*itRow)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->storagetype() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->defaultMatrix(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->defaultMatrix()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->defaultMatrix()->size(1) == dim, true);
  unsigned int i = 0, pos = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getSizeOfDiagonalBlock(i++) == (*it)->nonSmoothLaw()->size(), true);
    pos += (*it)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", (*M->defaultMatrix() - MRef).normInf() < tol, true);
  // Start from matrix with maxSize = M and and fill it (with resize)
  M.reset(new OSNSMatrix(30));
  M->fill(indexSet, blocks);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->storagetype() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->defaultMatrix(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->defaultMatrix()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->defaultMatrix()->size(1) == dim, true);
  i = 0;
  pos = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getSizeOfDiagonalBlock(i++) == (*it)->nonSmoothLaw()->size(), true);
    pos += (*it)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", (*M->defaultMatrix() - MRef).normInf() < tol, true);
  std::cout << "------- fill function test ended with success -------" <<std::endl;
  std::cout <<std::endl <<std::endl;
}

void OSNSMatrixTest::testConvert()
{
  std::cout << "------- convert function test -------" <<std::endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  M->convert();
  SP::NumericsMatrix NumMat = M->numericsMatrix();
  unsigned int dim = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->nonSmoothLaw()->size();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (InteractionsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (InteractionsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->nonSmoothLaw()->size();
    }
    col = 0;
    row += (*itRow)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", NumMat->storageType == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", NumMat->size0 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", NumMat->size1 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", !NumMat->matrix1, false);
  double * m1 = NumMat->matrix0;
  double *mRef = MRef.getArray();
  for (unsigned int k = 0; k < dim * dim; k++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", fabs(mRef[k] - m1[k]) < tol, true);

  std::cout << "------- convert function test ended with success -------" <<std::endl;
  std::cout <<std::endl <<std::endl;
}

void OSNSMatrixTest::testFill2()
{
  std::cout << "------- fill2 function test (sparse storage) -------" <<std::endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  unsigned int dim = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->nonSmoothLaw()->size();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (InteractionsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (InteractionsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->nonSmoothLaw()->size();
    }
    col = 0;
    row += (*itRow)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->storagetype() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->defaultMatrix(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->defaultMatrix()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->defaultMatrix()->size(1) == dim, true);
  unsigned int i = 0, pos = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getSizeOfDiagonalBlock(i++) == (*it)->nonSmoothLaw()->size(), true);
    pos += (*it)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", (*M->defaultMatrix() - MRef).normInf() < tol, true);

  // Start from matrix with maxSize = M and and fill it (with resize)
  M.reset(new OSNSMatrix(30));
  M->fill(indexSet, blocks);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->storagetype() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->defaultMatrix(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->defaultMatrix()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->defaultMatrix()->size(1) == dim, true);
  i = 0;
  pos = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getSizeOfDiagonalBlock(i++) == (*it)->nonSmoothLaw()->size(), true);
    pos += (*it)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", (*M->defaultMatrix() - MRef).normInf() < tol, true);

  std::cout << "------- fill2 function test ended with success -------" <<std::endl;
  std::cout <<std::endl <<std::endl;
}

void OSNSMatrixTest::testConvert2()
{
  std::cout << "------- convert2 function test -------" <<std::endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  M->convert();
  SP::NumericsMatrix NumMat = M->numericsMatrix();
  unsigned int dim = 0;
  for (InteractionsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->nonSmoothLaw()->size();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (InteractionsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (InteractionsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->nonSmoothLaw()->size();
    }
    col = 0;
    row += (*itRow)->nonSmoothLaw()->size();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", NumMat->storageType == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", NumMat->size0 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", NumMat->size1 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", !NumMat->matrix1, false);
  double * m1 = NumMat->matrix0;
  double *mRef = MRef.getArray();
  for (unsigned int k = 0; k < dim * dim; k++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", fabs(mRef[k] - m1[k]) < tol, true);

  std::cout << "------- convert2 function test ended with success -------" <<std::endl;
  std::cout <<std::endl <<std::endl;
}

void OSNSMatrixTest::End()
{
  std::cout << "==========================================" <<std::endl;
  std::cout << " ===== End of OSNSMatrix tests ===== " <<std::endl;
  std::cout << "==========================================" <<std::endl;
}
