/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "OSNSMatrixTest.h"
using namespace std;

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
  SP::TimeStepping s = boost::static_pointer_cast<TimeStepping>(temp->getSimulationPtr());
  s->initialize();
  // Get a set of Unitary Relations
  indexSet = s->getIndexSetPtr(0);
  SP::OneStepNSProblem osns = s->getOneStepNSProblems()->begin()->second;
  osns->computeAllBlocks();
  blocks = osns->getBlocks();
}

void OSNSMatrixTest::tearDown()
{}

// Default constructor
void OSNSMatrixTest::testBuildOSNSMatrix0()
{
  cout << "===========================================" << endl;
  cout << " ===== OSNSMatrix tests start ...===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Default constructor test -------" << endl;
  SP::OSNSMatrix  M = new OSNSMatrix();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix0 : ", M->size() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix0 : ", M->getStorageType() == 0, true);
  cout << "------- Default constructor OSNSMatrix ok -------" << endl;
  cout << endl << endl;
}

// Basic constructor
void OSNSMatrixTest::testBuildOSNSMatrix1()
{
  cout << "------- Constructor with dim. test -------" << endl;
  SP::OSNSMatrix  M(new OSNSMatrix(n));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->size() == n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->getStorageType() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->getDefaultMatrixPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->getDefaultMatrixPtr()->size(0) == n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->getDefaultMatrixPtr()->size(1) == n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix1 : ", M->getDefaultMatrixPtr()->normInf() < tol , true);
  cout << "------- Constructor(dim) ended with success -------" << endl;
  cout << endl << endl;
}

void OSNSMatrixTest::testBuildOSNSMatrix2()
{
  cout << "------- Constructor with index set and list of blocks -------" << endl;

  SP::OSNSMatrix  M(new OSNSMatrix(indexSet, blocks));

  unsigned int dim = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->getNonSmoothLawSize();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (UnitaryRelationsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (UnitaryRelationsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->getNonSmoothLawSize();
    }
    col = 0;
    row += (*itRow)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getStorageType() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getDefaultMatrixPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getDefaultMatrixPtr()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getDefaultMatrixPtr()->size(1) == dim, true);
  unsigned int i = 0, pos = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", M->getSizeOfDiagonalBlock(i++) == (*it)->getNonSmoothLawSize(), true);
    pos += (*it)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildOSNSMatrix2 : ", (*M->getDefaultMatrixPtr() - MRef).normInf() < tol, true);

  cout << "------- Constructor(indexSet,blocks) ended with success -------" << endl;
  cout << endl << endl;
}

void OSNSMatrixTest::testFill()
{
  cout << "------- fill function test -------" << endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  unsigned int dim = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->getNonSmoothLawSize();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (UnitaryRelationsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (UnitaryRelationsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->getNonSmoothLawSize();
    }
    col = 0;
    row += (*itRow)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getStorageType() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getDefaultMatrixPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getDefaultMatrixPtr()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getDefaultMatrixPtr()->size(1) == dim, true);
  unsigned int i = 0, pos = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getSizeOfDiagonalBlock(i++) == (*it)->getNonSmoothLawSize(), true);
    pos += (*it)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", (*M->getDefaultMatrixPtr() - MRef).normInf() < tol, true);
  // Start from matrix with maxSize = M and and fill it (with resize)
  M.reset(new OSNSMatrix(30));
  M->fill(indexSet, blocks);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getStorageType() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getDefaultMatrixPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getDefaultMatrixPtr()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getDefaultMatrixPtr()->size(1) == dim, true);
  i = 0;
  pos = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", M->getSizeOfDiagonalBlock(i++) == (*it)->getNonSmoothLawSize(), true);
    pos += (*it)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", (*M->getDefaultMatrixPtr() - MRef).normInf() < tol, true);
  cout << "------- fill function test ended with success -------" << endl;
  cout << endl << endl;
}

void OSNSMatrixTest::testConvert()
{
  cout << "------- convert function test -------" << endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  M->convert();
  SP::NumericsMatrix NumMat = M->getNumericsMatrix();
  unsigned int dim = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->getNonSmoothLawSize();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (UnitaryRelationsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (UnitaryRelationsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->getNonSmoothLawSize();
    }
    col = 0;
    row += (*itRow)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", NumMat->storageType == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", NumMat->size0 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", NumMat->size1 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", !NumMat->matrix1, false);
  double * m1 = NumMat->matrix0;
  double *mRef = MRef.getArray();
  for (unsigned int k = 0; k < dim * dim; k++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert : ", fabs(mRef[k] - m1[k]) < tol, true);

  cout << "------- convert function test ended with success -------" << endl;
  cout << endl << endl;
}

void OSNSMatrixTest::testFill2()
{
  cout << "------- fill2 function test (sparse storage) -------" << endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  unsigned int dim = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->getNonSmoothLawSize();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (UnitaryRelationsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (UnitaryRelationsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->getNonSmoothLawSize();
    }
    col = 0;
    row += (*itRow)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getStorageType() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getDefaultMatrixPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getDefaultMatrixPtr()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getDefaultMatrixPtr()->size(1) == dim, true);
  unsigned int i = 0, pos = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getSizeOfDiagonalBlock(i++) == (*it)->getNonSmoothLawSize(), true);
    pos += (*it)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", (*M->getDefaultMatrixPtr() - MRef).normInf() < tol, true);

  // Start from matrix with maxSize = M and and fill it (with resize)
  M.reset(new OSNSMatrix(30));
  M->fill(indexSet, blocks);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->size() == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getStorageType() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getDefaultMatrixPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getDefaultMatrixPtr()->size(0) == dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getDefaultMatrixPtr()->size(1) == dim, true);
  i = 0;
  pos = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getPositionOfBlock(*it) == pos, true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", M->getSizeOfDiagonalBlock(i++) == (*it)->getNonSmoothLawSize(), true);
    pos += (*it)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill2 : ", (*M->getDefaultMatrixPtr() - MRef).normInf() < tol, true);

  cout << "------- fill2 function test ended with success -------" << endl;
  cout << endl << endl;
}

void OSNSMatrixTest::testConvert2()
{
  cout << "------- convert2 function test -------" << endl;
  // Start from empty matrix and fill it
  SP::OSNSMatrix  M(new OSNSMatrix());
  M->fill(indexSet, blocks);

  M->convert();
  SP::NumericsMatrix NumMat = M->getNumericsMatrix();
  unsigned int dim = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
    dim += (*it)->getNonSmoothLawSize();
  SimpleMatrix MRef(dim, dim);
  int row = 0, col = 0;
  for (UnitaryRelationsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (UnitaryRelationsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      if (blocks[*itRow].find(*itCol) == blocks[*itRow].end())
      {}
      else
        MRef.setBlock(row, col, *(blocks[*itRow][*itCol]));
      col += (*itCol)->getNonSmoothLawSize();
    }
    col = 0;
    row += (*itRow)->getNonSmoothLawSize();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", NumMat->storageType == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", NumMat->size0 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", NumMat->size1 == (int)dim, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", !NumMat->matrix1, false);
  double * m1 = NumMat->matrix0;
  double *mRef = MRef.getArray();
  for (unsigned int k = 0; k < dim * dim; k++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConvert2 : ", fabs(mRef[k] - m1[k]) < tol, true);

  cout << "------- convert2 function test ended with success -------" << endl;
  cout << endl << endl;
}

void OSNSMatrixTest::End()
{
  cout << "==========================================" << endl;
  cout << " ===== End of OSNSMatrix tests ===== " << endl;
  cout << "==========================================" << endl;
}
