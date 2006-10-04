/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
#include "BlockMatrix.h"
using namespace std;


// --- CONSTRUCTORS ---

// Default (private)
BlockMatrix::BlockMatrix(): SiconosMatrix(true)
{}

void BlockMatrix::addInTab(const unsigned int& i, const unsigned int& j)
{
  tabRow.push_back(i);
  tabCol.push_back(j);
  unsigned int dim = tabRow.size();
  if (tabRow.size() > 1)
  {
    tabRow[dim - 1] += tabRow[dim - 2];
    tabCol[dim - 1] += tabCol[dim - 2];
  }
}

void BlockMatrix::makeTab(const unsigned int& row, const unsigned int& col)
{
  //blocksMap::iterator it;
  tabRow.clear();
  tabCol.clear();
  vector<unsigned int> index(2);
  index[0] = 0;
  index[1] = 0;
  unsigned int dim;

  // col sweeping
  for (index[1] = 0; index[1] < col ; ++index[1])
  {
    tabCol.push_back(matrixOfBlocks[index]->size(1));
    dim = tabCol.size();
    if (dim > 1)
      tabCol[dim - 1] += tabCol[dim - 2];
  }

  index[1] = 0;

  // row sweeping
  for (index[0] = 0; index[0] < row ; ++index[0])
  {
    tabRow.push_back(matrixOfBlocks[index]->size(0));
    dim = tabRow.size();
    if (dim > 1)
      tabRow[dim - 1] += tabRow[dim - 2];
  }
}

// copy
BlockMatrix::BlockMatrix(const SiconosMatrix& m): SiconosMatrix(true)
{
  if (!m.isBlock())
    SiconosMatrixException::selfThrow("BlockMatrix copy constructor from a SimpleMatrix: forbidden operation.");

  isBlockAllocatedIn.resize(m.size(0)*m.size(1));

  blocksMap::iterator it;
  blocksMap mBlocks = (static_cast<BlockMatrix>(m)).getListOfBlocks();

  tabRow = (static_cast<BlockMatrix>(m)).getTabRow();
  tabCol = (static_cast<BlockMatrix>(m)).getTabCol();

  unsigned int i, j, nCol = tabCol.size();

  for (it = mBlocks.begin(); it != mBlocks.end(); ++it)
  {
    matrixOfBlocks[it->first] = new SimpleMatrix(*(it->second));
    i = (it->first)[0];
    j = (it->first)[1];
    isBlockAllocatedIn[i * nCol + j] = true;
  }
}

BlockMatrix::BlockMatrix(const BlockMatrix& m): SiconosMatrix(true)
{

  isBlockAllocatedIn.resize(m.size(0)*m.size(1));

  blocksMap::iterator it;
  blocksMap mBlocks = m.getListOfBlocks();

  tabRow = m.getTabRow();
  tabCol = m.getTabCol();

  unsigned int i, j, nCol = tabCol.size();

  for (it = mBlocks.begin(); it != mBlocks.end(); ++it)
  {
    matrixOfBlocks[it->first] = new SimpleMatrix(*(it->second));
    i = (it->first)[0];
    j = (it->first)[1];
    isBlockAllocatedIn[i * nCol + j] = true;
  }
}

// copy a list of LaGenMatDouble
BlockMatrix::BlockMatrix(const std::vector<LaGenMatDouble>& m, const unsigned int& row, const unsigned int& col):
  SiconosMatrix(true)
{
  if (m.size() != (row * col))
    SiconosMatrixException::selfThrow("BlockMatrix constructor from a vector<LaGenMatDouble>, number of blocks inconsistent with provided dimensions.");

  isBlockAllocatedIn.resize(row * col);
  tabRow.reserve(row);
  tabCol.reserve(col);

  vector<unsigned int> index(2);
  index[0] = 0;
  index[1] = 0;

  for (unsigned int i = 0; i < row; ++i)
  {
    for (unsigned int j = 0 ; j < col ; ++j)
    {
      matrixOfBlocks[index] = new SimpleMatrix(m[i * col + j]);
      isBlockAllocatedIn.push_back(true);

      index[1]++;
    }
    index[0]++;
    index[1] = 0;
  }

  makeTab(row, col);

}

BlockMatrix::BlockMatrix(std::vector<SiconosMatrix*> m, const unsigned int& row, const unsigned int& col):
  SiconosMatrix(true)
{
  if (m.size() != (row * col))
    SiconosMatrixException::selfThrow("BlockMatrix constructor from a vector<SiconosMatrix*>, number of blocks inconsistent with provided dimensions.");

  isBlockAllocatedIn.resize(row * col);

  vector<unsigned int> index(2);
  index[0] = 0;
  index[1] = 0;

  for (unsigned int i = 0; i < row; ++i)
  {
    for (unsigned int j = 0 ; j < col ; ++j)
    {
      matrixOfBlocks[index] = m[i * col + j];
      isBlockAllocatedIn.push_back(false);

      index[1]++;
    }
    index[0]++;
    index[1] = 0;
  }
  makeTab(row, col);
}

BlockMatrix::BlockMatrix(SiconosMatrix* m1, SiconosMatrix* m2, SiconosMatrix* m3, SiconosMatrix* m4):
  SiconosMatrix(true)
{
  isBlockAllocatedIn.resize(4, false);

  vector<unsigned int> index(2);
  index[0] = 0;
  index[1] = 0;
  matrixOfBlocks[index] = m1;
  index[0] = 0;
  index[1] = 1;
  matrixOfBlocks[index] = m2;
  index[0] = 1;
  index[1] = 0;
  matrixOfBlocks[index] = m3;
  index[0] = 1;
  index[1] = 1;
  matrixOfBlocks[index] = m4;
  makeTab(2, 2);
}

// --- DESTRUCTOR ---

BlockMatrix::~BlockMatrix()
{
  blocksMap::iterator it;
  unsigned int col = tabCol.size();
  vector<unsigned int> index(2);
  unsigned int pos ;
  for (it = matrixOfBlocks.begin(); it != matrixOfBlocks.end(); ++it)
  {
    index = it->first;
    pos = index[0] * col + index[1];
    if (isBlockAllocatedIn[pos]) delete it->second;
    it->second = NULL;
  }
}

// --- FUNCTIONS TO GET INFO ABOUT THE MATRIX ---

unsigned int BlockMatrix::size(const unsigned int& d) const
{
  if ((d != 0) && (d != 1))
    SiconosMatrixException::selfThrow("function size() : Index out of range");

  if (!d)
    return tabRow[tabRow.size() - 1];
  else
    return tabCol[tabCol.size() - 1];
}

bool BlockMatrix::isSquare() const
{
  return (size(0) == size(1));
}

bool BlockMatrix::isInversed() const
{
  SiconosMatrixException::selfThrow("BlockMatrix isInversed: not implemented.");
  return false;
}

bool BlockMatrix::isFactorized() const
{
  SiconosMatrixException::selfThrow("BlockMatrix isFactorized: not implemented.");
  return false;
}

// --- GETTERS/SETTERS ---

const LaGenMatDouble BlockMatrix::getLaGenMatDouble(const unsigned int& row, const unsigned int& col) const
{
  vector<unsigned int> index(2);
  index[0] = row;
  index[1] = col;
  return (getListOfBlocks()[index]->getLaGenMatDouble());
}

const LaGenMatDouble* BlockMatrix::getLaGenMatDoubleRef(const unsigned int& row, const unsigned int& col) const
{
  vector<unsigned int> index(2);
  index[0] = row;
  index[1] = col;
  return (getListOfBlocks()[index]->getLaGenMatDoubleRef());
}

void BlockMatrix::setValue(const LaGenMatDouble& newMat, const unsigned int& row, const unsigned int& col)
{
  vector<unsigned int> index(2);
  index[0] = row;
  index[1] = col;
  matrixOfBlocks[index]->setValue(newMat);
}

void BlockMatrix::setRow(const unsigned int& row, const SiconosVector &v)
{
  SiconosMatrixException::selfThrow("BlockMatrix setRow: not implemented.");
}

void BlockMatrix::getRow(const unsigned int& index, const SimpleVector& vOut) const
{
  SiconosMatrixException::selfThrow("BlockMatrix getRow: not implemented.");
}

void BlockMatrix::setCol(const unsigned int& col, const SiconosVector &v)
{
  SiconosMatrixException::selfThrow("BlockMatrix setCol: not implemented.");
}

void BlockMatrix::getCol(const unsigned int& index, const SimpleVector& vOut) const
{
  SiconosMatrixException::selfThrow("BlockMatrix getCol: not implemented.");
}

void BlockMatrix::getBlock(const vector<unsigned int>& index_list, SiconosMatrix& block) const
{
  SiconosMatrixException::selfThrow("BlockMatrix getBlock: not implemented.");
}

void BlockMatrix::getBlock(const vector<unsigned int>& indexRow, const vector<unsigned int>& indexCol, SiconosMatrix& block) const
{
  SiconosMatrixException::selfThrow("BlockMatrix getBlock: not implemented.");
}

SiconosMatrix* BlockMatrix::getBlockPtr(const unsigned int& row, const unsigned int& col)
{
  vector<unsigned int> index(2);
  index[0] = row;
  index[1] = col;
  return matrixOfBlocks[index];
}

double* BlockMatrix::getArray(const unsigned int& row, const unsigned int& col)
{
  vector<unsigned int> index(2);
  index[0] = row;
  index[1] = col;
  return (matrixOfBlocks[index])->getArray();
}

// --- READ, WRITE ... ---

bool BlockMatrix::read(const string& fileName, const string& mode)
{
  SiconosMatrixException::selfThrow("BlockMatrix read: not implemented.");
  return false;
}

bool BlockMatrix::write(const string& fileName, const string& mode) const
{
  if ((mode != "binary") && (mode != "ascii"))
    SiconosMatrixException::selfThrow("BlockMatrix write: incorrect mode for writing");

  // open the file
  ofstream outFile(fileName.c_str());           // checks that it's opened
  if (!outFile.is_open())
    SiconosMatrixException::selfThrow("function write error : Fail to open \"" + fileName + "\"");

  unsigned int m = size(0);
  unsigned int n = size(1);
  if (mode == "binary")
  {
    outFile.write((char*)&m, sizeof(unsigned int));
    outFile.write((char*)&n, sizeof(unsigned int));
  }
  else if (mode == "ascii")
    outFile << m << ' ' << n;

  for (unsigned int i = 0; i < m; i++)
  {
    if (mode == "ascii")
    {
      outFile << endl;
    }
    for (unsigned int j = 0; j < n; j++)
    {
      if (mode == "binary")
      {
        outFile.write((char*) & (*this)(i, j), sizeof(double));
      }
      else if (mode == "ascii")
      {
        char buffer[30];
        sprintf(buffer, "%1.17e ", (*this)(i, j)); // /!\ depends on machine precision
        outFile << buffer;
      }
    }
  }
  outFile.close();
  return true;
}

bool BlockMatrix::rawWrite(const string& fileName, const string& mode) const
{
  //  if( (size(0) == 0)||(size(1) == 0) ) SiconosMatrixException::selfThrow("write impossible - SiconosMatrix empty");

  if ((mode != "binary") && (mode != "ascii"))
    SiconosMatrixException::selfThrow("Incorrect mode for writing");

  // open the file
  ofstream outFile(fileName.c_str());           // checks that it's opened
  if (!outFile.is_open())
    SiconosMatrixException::selfThrow("function write error : Fail to open \"" + fileName + "\"");

  unsigned int m = size(0);
  unsigned int n = size(1);

  for (unsigned int i = 0; i < m; i++)
  {
    if (mode == "ascii")
    {
      outFile << endl;
    }
    for (unsigned int j = 0; j < n; j++)
    {
      if (mode == "binary")
      {
        outFile.write((char*) & (*this)(i, j), sizeof(double));
      }
      else if (mode == "ascii")
      {
        char buffer[30];
        sprintf(buffer, "%1.17e ", (*this)(i, j)); // /!\ depends on machine precision
        outFile << buffer;
      }
    }
  }
  outFile.close();
  return true;
}

void BlockMatrix::zero()
{
  blocksMap::iterator it;
  for (it = matrixOfBlocks.begin(); it != matrixOfBlocks.end(); ++it)
    (it->second)->zero();
}

void BlockMatrix::eye()
{
  blocksMap::iterator it;
  for (it = matrixOfBlocks.begin(); it != matrixOfBlocks.end(); ++it)
    (it->second)->eye();
}

void BlockMatrix::display() const
{
  cout << "=======> BlockMatrix, with " << size(0) << " lines and " << size(1) << " columns." << endl;
  cout << " There are " << tabRow.size() << " blocks in each row, and " << tabCol.size() << " blocks in each column." << endl;
  cout << "For details, call display() of the specific block you want." << endl;
  cout << "===== End of BlockMatrix display ====" << endl;
}

// --- MATRICES HANDLING AND OPERATORS ---

SimpleMatrix BlockMatrix::multTranspose(const SiconosMatrix & B)
{
  SiconosMatrixException::selfThrow("BlockMatrix multTranspose: not implemented.");
  return B;
}

void BlockMatrix::blockMatrixCopy(const SiconosMatrix &blockMat, const unsigned int& xPos, const unsigned int& yPos)
{
  SiconosMatrixException::selfThrow("BlockMatrix blockMatrixCopy: not implemented.");
}

/****************** () ******************/
// subscript operator to get/set individual elements
double& BlockMatrix::operator()(const int& row, const int& col)
{
  if ((row >= (int)size(0)) || (col >= (int)size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  unsigned int numBlockRow = 0, numBlockCol = 0;
  while ((unsigned int)row >= tabRow[numBlockRow] && numBlockRow < tabRow.size()) numBlockRow++;
  while ((unsigned int)col >= tabCol[numBlockCol] && numBlockCol < tabCol.size()) numBlockCol++;

  vector<unsigned int> index(2);
  index[0] = numBlockRow;
  index[1] = numBlockCol;

  unsigned int posRow = row , posCol = col;
  if (numBlockRow != 0)
    posRow -= tabRow[numBlockRow - 1];
  if (numBlockCol != 0)
    posCol -= tabCol[numBlockCol - 1];

  return (*matrixOfBlocks[index])(posRow, posCol);
}

// subscript operator to get/set individual elements
double& BlockMatrix::operator()(const unsigned int& row, const unsigned int& col)
{
  if ((row >= size(0)) || (col >= size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  unsigned int numBlockRow = 0, numBlockCol = 0;
  while (row >= tabRow[numBlockRow] && numBlockRow < tabRow.size()) numBlockRow++;
  while (col >= tabCol[numBlockCol] && numBlockCol < tabCol.size()) numBlockCol++;

  vector<unsigned int> index(2);
  index[0] = numBlockRow;
  index[1] = numBlockCol;

  unsigned int posRow = row, posCol = col;
  if (numBlockRow != 0)
    posRow -= tabRow[numBlockRow - 1];
  if (numBlockCol != 0)
    posCol -= tabCol[numBlockCol - 1];

  return (*matrixOfBlocks[index])(posRow, posCol);
}

// subscript operator to get/set individual elements
double& BlockMatrix::operator()(const unsigned int& row, const unsigned int& col) const
{
  if ((row >= size(0)) || (col >= size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  unsigned int numBlockRow = 0, numBlockCol = 0;
  while (row >= tabRow[numBlockRow] && numBlockRow < tabRow.size()) numBlockRow++;
  while (col >= tabCol[numBlockCol] && numBlockCol < tabCol.size()) numBlockCol++;

  vector<unsigned int> index(2);
  index[0] = numBlockRow;
  index[1] = numBlockCol;

  unsigned int posRow = row, posCol = col;
  if (numBlockRow != 0)
    posRow -= tabRow[numBlockRow - 1];
  if (numBlockCol != 0)
    posCol -= tabCol[numBlockCol - 1];
  return (getListOfBlocks()[index]->getLaGenMatDouble())(posRow, posCol);
}

/*************************************************/
BlockMatrix& BlockMatrix::operator = (const SiconosMatrix& m)
{

  if (&m == this) return *this;

  if (m.size(0) != size(0) || m.size(1) != size(1))
    SiconosMatrixException::selfThrow("operator = : left and right value have inconsistent sizes.");

  if (m.isBlock())
  {
    blocksMap::iterator it;
    blocksMap mBlocks = (static_cast<BlockMatrix>(m)).getListOfBlocks();
    for (it = mBlocks.begin(); it != mBlocks.end(); ++it)
      matrixOfBlocks[it->first] = new SimpleMatrix(*(it->second));
    tabRow = (static_cast<BlockMatrix>(m)).getTabRow();
    tabCol = (static_cast<BlockMatrix>(m)).getTabCol();
  }
  else
    SiconosMatrixException::selfThrow("operator = : BlockMatrix = SimpleMatrix is a forbidden operation.");

  return *this;
}

BlockMatrix& BlockMatrix::operator = (const BlockMatrix& m)
{
  if (&m == this) return *this;

  if (m.size(0) != size(0) || m.size(1) != size(1))
    SiconosMatrixException::selfThrow("operator = : left and right value have inconsistent sizes.");

  if (!m.isBlock())
    SiconosMatrixException::selfThrow("BlockMatrix copy constructor from a SimpleMatrix: forbidden operation.");

  isBlockAllocatedIn.clear();
  isBlockAllocatedIn.resize(m.size(0)*m.size(1));
  tabRow = m.getTabRow();
  tabCol = m.getTabCol();

  blocksMap::iterator it;
  blocksMap mBlocks = m.getListOfBlocks();

  unsigned int i, j, nCol = tabCol.size();

  for (it = mBlocks.begin(); it != mBlocks.end(); ++it)
  {
    matrixOfBlocks[it->first] = new SimpleMatrix(*(it->second));
    i = (it->first)[0];
    j = (it->first)[1];
    isBlockAllocatedIn[i * nCol + j] = true;
  }

  return *this;
}

BlockMatrix& BlockMatrix::operator+=(const SiconosMatrix &m)
{
  // Check global sizes
  if (m.size(0) != size(0) || m.size(1) != size(1))
    SiconosMatrixException::selfThrow("operator = : left and right value have inconsistent sizes.");

  if (m.isBlock())
  {
    // get list of blocks of m and iterate over them
    blocksMap mBlocks = (static_cast<BlockMatrix>(m)).getListOfBlocks();
    blocksMap::iterator it;
    for (it = mBlocks.begin(); it != mBlocks.end(); ++it)
      *matrixOfBlocks[it->first] += *(it->second) ;
    // Note that size consistency is checked in += operator of SimpleMatrix.
  }
  else
  {
    // sum using () operator
    for (unsigned int i = 0; i < size(0); ++i)
      for (unsigned int j = 0; j < size(1); ++j)
        (*this)(i, j) += m(i, j);
  }

  return *this;
}

BlockMatrix& BlockMatrix::operator-=(const SiconosMatrix &m)
{
  // Check global sizes
  if (m.size(0) != size(0) || m.size(1) != size(1))
    SiconosMatrixException::selfThrow("operator = : left and right value have inconsistent sizes.");

  if (m.isBlock())
  {
    // get list of blocks of m and iterate over them
    blocksMap mBlocks = (static_cast<BlockMatrix>(m)).getListOfBlocks();
    blocksMap::iterator it;
    for (it = mBlocks.begin(); it != mBlocks.end(); ++it)
      *matrixOfBlocks[it->first] -= *(it->second) ;
    // Note that size consistency is checked in += operator of SimpleMatrix.
  }
  else
  {
    // sum using () operator
    for (unsigned int i = 0; i < size(0); ++i)
      for (unsigned int j = 0; j < size(1); ++j)
        (*this)(i, j) -= m(i, j);
  }

  return *this;
}

BlockMatrix& BlockMatrix::operator*=(const double& d)
{
  blocksMap::iterator it;
  for (it = matrixOfBlocks.begin(); it != matrixOfBlocks.end(); ++it)
    *(it->second) *= d;

  return *this;
}

const double BlockMatrix::normInf() const
{
  double sum = 0, norm = 0;
  for (unsigned int i = 0; i < size(0); ++i)
  {
    for (unsigned int j = 0; j < size(0); ++j)
      sum += (*this)(i, j);
    if (sum > norm) norm = sum;
    sum = 0;
  }
  return norm;
}

// --- COMPUTING WITH MATRICES  ---

void BlockMatrix::linearSolve(const SiconosMatrix &B, SiconosMatrix &X)
{
  SiconosMatrixException::selfThrow("BlockMatrix linearSolve: not implemented.");
};

void  BlockMatrix::PLUFactorizationInPlace()
{
  SiconosMatrixException::selfThrow("BlockMatrix PLUFactorizationInPlace: not implemented.");
}

void  BlockMatrix::PLUInverseInPlace()
{
  SiconosMatrixException::selfThrow("BlockMatrix PLUInverseInPlace: not implemented.");
}

void  BlockMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  SiconosMatrixException::selfThrow("BlockMatrix PLUForwardBackwardInPlace: not implemented.");
}

void BlockMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  SiconosMatrixException::selfThrow("BlockMatrix PLUForwardBackwardInPlace: not implemented.");
}

