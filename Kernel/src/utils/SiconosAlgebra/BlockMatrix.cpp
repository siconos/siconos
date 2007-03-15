#include "BlockMatrix.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"
#include "SiconosMatrixException.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

void BlockMatrix::addInTab(unsigned int row, unsigned int col)
{

  tabRow.push_back(row);
  tabCol.push_back(col);

  unsigned int nbRow = tabRow.size();
  unsigned int nbCol = tabCol.size();

  if (nbRow > 1)
    tabRow[nbRow - 1] += tabRow[nbRow - 2];
  if (nbCol > 1)
    tabCol[nbCol - 1] += tabCol[nbCol - 2];
}

void BlockMatrix::makeTab(unsigned int row, unsigned int col)
{

  unsigned int dim;
  tabRow.clear();
  tabCol.clear();

  if (STDMAP == 1)
  {
    // Col sweeping
    for (unsigned int j = 0; j < map.size2(); j++)
    {
      tabCol.push_back((*map(0, j)).size(1));
      dim = tabCol.size();
      if (dim > 1)
        tabCol[dim - 1] += tabCol[dim - 2];
    }
    // Row sweeping
    for (unsigned int i = 0; i < map.size1(); i++)
    {
      tabRow.push_back((*map(i, 0)).size(0));
      dim = tabRow.size();
      if (dim > 1)
        tabRow[dim - 1] += tabRow[dim - 2];
    }
  }
  else
  {
    //BlocksMat Mmap = map;
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;

    // Col sweeping
    for (it2 = map.begin2(); it2 != map.end2(); ++it2)
    {
      tabCol.push_back((**it2).size(1));
      dim = tabCol.size();
      if (dim > 1)
        tabCol[dim - 1] += tabCol[dim - 2];
    }

    // Row sweeping
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      tabRow.push_back((**it).size(0));
      dim = tabRow.size();
      if (dim > 1)
        tabRow[dim - 1] += tabRow[dim - 2];
    }
  }
  // Update matrix dim.
  computeDim();
}

// Default (private)
BlockMatrix::BlockMatrix(): SiconosMatrix(true)
{}

BlockMatrix::BlockMatrix(const SiconosMatrix &m): SiconosMatrix(true)
{
  if (!m.isBlock())
    SiconosMatrixException::selfThrow("BlockMatrix copy constructor from a SimpleMatrix: forbidden operation.");

  // BlocksMat Mmap = (dynamic_cast<const BlockMatrix&>(m)).map;
  BlocksMat Mmap = m.getAllBlocks();
  unsigned int col = Mmap.size2();
  unsigned int row = Mmap.size1();
  isBlockAllocatedIn.resize(col * row);

  map.resize(row, col, false);
  unsigned int i = 0, j = 0;
  // STDMAP = 1 means that we use the map::iterator of standard library, else we use the iterator of boost map

  if (STDMAP == 1)
  {
    //       BlocksMat::array_type::iterator it;

    //       for(it=(Mmap.data ()).begin (); it!=(Mmap.data ()).end (); ++it){
    //  // Computation of index row and index column corresponding
    //  i = (it->first)/col;
    //  j = (it->first) - i * col;
    //  map(i, j) = new SimpleMatrix( *(it->second) );
    //  isBlockAllocatedIn[ (it->first) ] = true;
    //       }
  }

  else
  {
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;

    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2 ++)
      {
        // Computation of index row and index column corresponding
        i = it2.index1();
        j = it2.index2();
        map(i, j) = new SimpleMatrix(**it2);
        isBlockAllocatedIn[i * col + j] = true;
      }
    }
  }
  makeTab(row, col);
}

BlockMatrix::BlockMatrix(const BlockMatrix &m): SiconosMatrix(true)
{
  // BlocksMat Mmap = (dynamic_cast<const BlockMatrix&>(m)).map;
  BlocksMat Mmap = m.getAllBlocks();
  unsigned int col = Mmap.size2();
  unsigned int row = Mmap.size1();
  isBlockAllocatedIn.resize(col * row);

  map.resize(row, col, false);
  unsigned int i = 0, j = 0;
  // STDMAP = 1 means that we use the map::iterator of standard library, else we use the iterator of boost map

  if (STDMAP == 1)
  {
    //       BlocksMat::array_type::iterator it;

    //       for(it=(Mmap.data ()).begin (); it!=(Mmap.data ()).end (); ++it){
    //  // Computation of index row and index column corresponding
    //  i = (it->first)/col;
    //  j = (it->first) - i * col;
    //  map(i, j) = new SimpleMatrix( *(it->second) );
    //  isBlockAllocatedIn[ (it->first) ] = true;
    //       }
  }

  else
  {
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;

    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2 ++)
      {
        // Computation of index row and index column corresponding
        i = it2.index1();
        j = it2.index2();
        map(i, j) = new SimpleMatrix(**it2);
        isBlockAllocatedIn[i * col + j] = true;
      }
    }
  }
  makeTab(row, col);
}

BlockMatrix::BlockMatrix(BlocksMat& Mmap): SiconosMatrix(true)
{

  unsigned int i = 0, j = 0;
  unsigned int row = Mmap.size1();
  unsigned int col = Mmap.size2();
  isBlockAllocatedIn.resize(row * col);
  map.resize(row, col, false);

  // STDMAP = 1 means that we use the map::iterator of standard library, else we use the iterator of boost map
  if (STDMAP == 1)
  {
    //     BlocksMat::array_type::iterator it1;

    //     for(it1=(Mmap.data ()).begin (); it1!=(Mmap.data ()).end (); ++it1){
    //       // Computation of index row and index column corresponding
    //       i = (it1->first)/col;
    //       j = (it1->first) - i * col;
    //       map(i, j) = new SimpleMatrix( *(it1->second) );
    //       isBlockAllocatedIn[ (it1->first) ] = true;
    //     }
  }
  else
  {

    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;

    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2 ++)
      {
        // Computation of index row and index column corresponding
        i = it2.index1();
        j = it2.index2();
        map(i, j) = new SimpleMatrix(**it2);
        isBlockAllocatedIn[i * col + j] = true;
      }
    }
  }
  makeTab(row, col);
}

BlockMatrix::BlockMatrix(const std::vector<SiconosMatrix* > &m, unsigned int row, unsigned int col): SiconosMatrix(true)
{
  if (m.size() != (row * col))
    SiconosMatrixException::selfThrow("BlockMatrix constructor from a vector<SiconosMatrix*>, number of blocks inconsistent with provided dimensions.");

  isBlockAllocatedIn.resize(row * col);
  map.resize(row, col, false);

  for (unsigned int i = 0; i < row; i++)
  {
    for (unsigned int j = 0; j < col; j++)
    {
      map(i, j) = m [i * col + j];
      isBlockAllocatedIn[i * col + j] = false;
    }
  }
  makeTab(row, col);
}

BlockMatrix::BlockMatrix(SiconosMatrix* A, SiconosMatrix* B, SiconosMatrix* C, SiconosMatrix* D): SiconosMatrix(true)
{
  if (A->size(0) != B->size(0) || C->size(0) != D->size(0) ||  A->size(1) != C->size(1) ||  B->size(1) != D->size(1))
    SiconosMatrixException::selfThrow("BlockMatrix constructor(A,B,C,D), inconsistent sizes between A, B, C or D SiconosMatrices.");

  isBlockAllocatedIn.resize(4, false);
  map.resize(2, 2, false);
  tabRow.reserve(2);
  tabCol.reserve(2);
  map(0, 0) = A;
  map(0, 1) = B;
  map(1, 0) = C;
  map(1, 1) = D;
  makeTab(2, 2);
}

BlockMatrix::~BlockMatrix(void)
{
  unsigned row = map.size1();
  unsigned col = map.size2();
  // we call delete on the map only if the pointer has been allocated in the class
  for (unsigned i = 0; i < row; ++i)
  {
    for (unsigned j = 0; j < col; ++j)
    {
      if (isBlockAllocatedIn[i * col + j] == true)
        delete(map(i, j));
    }
  }
}

// return the number of rows of blocks
void BlockMatrix::computeDim()
{
  dim[0] = tabRow[ tabRow.size() - 1 ];
  dim[1] = tabCol[ tabCol.size() - 1 ];
}

void BlockMatrix::resize(unsigned int row, unsigned int col, unsigned int lower, unsigned int upper, bool preserve)
{
  //   if(lower != 0 || upper != 0)
  //     SiconosMatrixException::selfThrow("BlockMatrix::resize : lower or upper not equal to 0.0");
  //   map.resize(row, col, preserve);
  SiconosMatrixException::selfThrow("BlockMatrix::resize : forbidden for block matrices.");
}

double* BlockMatrix::getArray(unsigned int i, unsigned int j) const
{
  return (map(i, j))->getArray();
}

BlockIterator1 BlockMatrix::begin()
{
  return map.begin1();
}

BlockIterator1 BlockMatrix::end()
{
  return map.end1();
}

ConstBlockIterator1 BlockMatrix::begin() const
{
  return map.begin1();
}

ConstBlockIterator1 BlockMatrix::end() const
{
  return map.end1();
}

const double BlockMatrix::normInf(void)const
{
  double sum = 0, norm = 0;
  for (unsigned int i = 0; i < size(0); i++)
  {
    for (unsigned int j = 0; j < size(1); j++)
    {
      sum += (*this)(i, j);
    }
    if (sum > norm) norm = sum;
    sum = 0;
  }
  return norm;
}

void BlockMatrix::trans()
{
  SiconosMatrixException::selfThrow("BlockMatrix::trans(): not yet implemented.");
}

void BlockMatrix::trans(const SiconosMatrix &m)
{
  SiconosMatrixException::selfThrow("BlockMatrix::trans(M): not yet implemented.");
}

void BlockMatrix::zero(void)
{
  if (STDMAP == 1)
  {
    //     BlocksMat::array_type::iterator it;
    //     for(it=(map.data ()).begin (); it!=(map.data ()).end (); ++it){
    //       (it->second)->zero ();
    //}
  }
  else
  {
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        (**it2).zero();
      }
    }
  }
}

void BlockMatrix::eye(void)
{
  if (STDMAP == 1)
  {
    //     unsigned int col = map.size2 ();
    //     unsigned int p =0, q =0;
    //     BlocksMat::array_type::iterator it;
    //     for(it=(map.data ()).begin (); it!=(map.data ()).end (); ++it){
    //       p = (it->first)/col;
    //       q = (it->first) - p * col;
    //       if(p == q){
    //  (it->second)->eye ();
    //       }
    //       else{
    //  (it->second)->zero ();
    //       }
    //     }
  }
  else
  {
    BlocksMat Mmap = map;
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;

    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        if ((it2.index1() == 0 && it2.index2() == 0) || (it2.index1() == it2.index2()))
        {
          (**it2).eye();
        }
        else
        {
          (**it2).zero();
        }
      }
    }
  }
}

void BlockMatrix::display(void)const
{
  std::cout << "==========> BlockMatrix (" << getNumberOfBlocks(0) << " blocks X " << getNumberOfBlocks(1) << " blocks): " << std::endl;
  if (STDMAP == 1)
  {
    //     BlocksMat::array_type Mmap = map.data ();
    //     BlocksMat::array_type::iterator it;
    //     for(it=Mmap.begin (); it!=Mmap.end (); ++it){
    //       (it->second)->display ();
    //     }
  }
  else
  {
    BlocksMat Mmap = map;
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        (**it2).display();
      }
    }
  }

}

void BlockMatrix::getBlock(unsigned int row, unsigned int col, SiconosMatrix &m)const
{
  SiconosMatrixException::selfThrow("BlockMatrix::getBlock, forbidden for BlockMatrices.");
  //m = dynamic_cast<SiconosMatrix&> ( *(map (row, col)) );
}

SiconosMatrix* BlockMatrix::getBlockPtr(unsigned int row, unsigned int col)
{
  return map(row, col);
}

const std::deque<bool> BlockMatrix::getBlockAllocated(void)const
{
  return isBlockAllocatedIn;
}

unsigned int BlockMatrix::getNum(void)const
{
  SiconosMatrixException::selfThrow("BlockMatrix::getNum of a block is forbidden.");
  return 1;
}

unsigned int BlockMatrix::getNumberOfBlocks(unsigned int dim) const
{
  if (dim == 1)
    return tabRow.size();
  else
    return tabCol.size();
}

void BlockMatrix::getRow(unsigned int r, SimpleVector &v) const
{
  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  if (r > size(0))
    SiconosMatrixException::selfThrow("BlockMatrix:getRow : row number is out of range");

  // Verification of the size of the result vector
  if (v.size() != size(1))
    SiconosMatrixException::selfThrow("BlockMatrix:getRow : inconsistent sizes");

  // Find the row-block number where "r" is
  while (r >= tabRow[numRow] && numRow < tabRow.size())
    numRow ++;

  // Computation of the value of the index row into this block
  if (numRow != 0)
    posRow -= tabRow[numRow - 1];

  for (unsigned int j = 0; j < tabCol.size(); j++)
  {
    start = stop;
    stop += (*map(numRow, j)).size(1);
    ublas::subrange(*(v.getDensePtr()), start, stop) = ublas::row((*map(numRow, j)).getDense(), posRow);
  }
}

void BlockMatrix::getCol(unsigned int c, SimpleVector &v) const
{
  unsigned int numCol = 0, posCol = c, start = 0, stop = 0;

  if (c > size(1))
    SiconosMatrixException::selfThrow("BlockMatrix:getCol : column number is out of range");

  // Verification of the size of the result vector
  if (v.size() != size(0))
    SiconosMatrixException::selfThrow("BlockMatrix:getcol : inconsistent sizes");

  // Find the column-block number where "c" is
  while (c >= tabCol[numCol] && numCol < tabCol.size())
    numCol ++;

  // Computation of the value of the index column into this block
  if (numCol != 0)
    posCol -= tabCol[numCol - 1];

  for (unsigned int i = 0; i < tabRow.size(); i++)
  {
    start = stop;
    stop += (*map(i, numCol)).size(0);
    ublas::subrange(*(v.getDensePtr()), start, stop) = ublas::column((*map(i, numCol)).getDense(), posCol);
  }
}

void BlockMatrix::setRow(unsigned int r, const SimpleVector &v)
{

  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  if (v.size() != size(1))
    SiconosMatrixException::selfThrow("BlockMatrix:setRow : inconsistent sizes");

  while (r >= tabRow[numRow] && numRow < tabRow.size())
    numRow ++;

  if (numRow != 0)
    posRow -= tabRow[numRow - 1];

  for (unsigned int j = 0; j < tabCol.size(); j++)
  {
    start = stop;
    stop += (*map(numRow, j)).size(1);
    ublas::row(*((*map(numRow, j)).getDensePtr()), posRow) = ublas::subrange(*(v.getDensePtr()), start, stop);
  }
}

void BlockMatrix::setCol(unsigned int col, const SimpleVector &v)
{

  unsigned int numCol = 0, posCol = col, start = 0, stop = 0;

  if (v.size() != size(0))
    SiconosMatrixException::selfThrow("BlockMatrix:setCol : inconsistent sizes");

  while (col >= tabCol[numCol] && numCol < tabCol.size())
    numCol ++;

  if (numCol != 0)
    posCol -= tabCol[numCol - 1];

  for (unsigned int i = 0; i < tabRow.size(); i++)
  {
    start = stop;
    stop += (*map(i, numCol)).size(0);
    ublas::column(*((*map(i, numCol)).getDensePtr()), posCol) = ublas::subrange(*(v.getDensePtr()), start, stop);
  }
}

// return the boost dense matrix of the block (i, j)
const DenseMat  BlockMatrix::getDense(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 1);
  SiconosMatrixException::selfThrow("DenseMat getDense(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Dense matrix");

  return (*map(row, col)).getDense();
}

// return the boost triangular matrix of the block (i, j)
const TriangMat BlockMatrix::getTriang(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 2);
  SiconosMatrixException::selfThrow("TriangMat getTriang(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Triangular matrix");
  return (*map(row, col)).getTriang();
}

// return the boost symmetric matrix of the block (i, j)
const SymMat BlockMatrix::getSym(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 3);
  SiconosMatrixException::selfThrow("SymMat getSym(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Symmmetric matrix");
  return (*map(row, col)).getSym();
}

// return the boost sparse matrix of the block (i, j)
const SparseMat  BlockMatrix::getSparse(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 4);
  SiconosMatrixException::selfThrow("SparseMat getSparse(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Sparse matrix");

  return (*map(row, col)).getSparse();
}

// return the boost banded matrix of the block (i, j)
const BandedMat  BlockMatrix::getBanded(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("BandedMat getBanded(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Banded matrix");

  return (*map(row, col)).getBanded();
}

// return the boost zero matrix of the block (i, j)
const ZeroMat  BlockMatrix::getZero(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("ZeroMat getZero(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Zero matrix");

  return (*map(row, col)).getZero();
}

// return the boost identity matrix of the block (i, j)
const IdentityMat  BlockMatrix::getIdentity(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("IdentityMat getIdentity(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Identity matrix");

  return (*map(row, col)).getIdentity();
}

// The following functions return the corresponding pointers
DenseMat*  BlockMatrix::getDensePtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 1);
  SiconosMatrixException::selfThrow("DenseMat* getDensePtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Dense matrix");

  return (*map(row, col)).getDensePtr();
}

TriangMat* BlockMatrix::getTriangPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 2);
  SiconosMatrixException::selfThrow("TriangMat* getTriangPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Triangular matrix");

  return (*map(row, col)).getTriangPtr();
}
SymMat* BlockMatrix::getSymPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 3);
  SiconosMatrixException::selfThrow("SymMat* getSymPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Symmmetric matrix");
  return (*map(row, col)).getSymPtr();
}

SparseMat*  BlockMatrix::getSparsePtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 4);
  SiconosMatrixException::selfThrow("SparseMat* getSparsePtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Sparse matrix");

  return (*map(row, col)).getSparsePtr();
}

BandedMat*  BlockMatrix::getBandedPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("BandedMat* getBandedPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Banded matrix");

  return (*map(row, col)).getBandedPtr();
}

ZeroMat*  BlockMatrix::getZeroPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("ZeroMat* getZeroPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Zero matrix");

  return (*map(row, col)).getZeroPtr();
}

IdentityMat*  BlockMatrix::getIdentityPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("IdentityMat* getIdentityPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Identity matrix");

  return (*map(row, col)).getIdentityPtr();
}

// return the boost BlocksMat matrix of a block matrix
const BlocksMat BlockMatrix::getAllBlocks(void)const
{
  return map;
}

void BlockMatrix::matrixCopy(const SiconosMatrix &m, unsigned int row, unsigned int col)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("BlockMatrix::matrixCopy of a block into an other block is forbidden.");

  if (row > tabRow.size() || col > tabCol.size())
    SiconosMatrixException::selfThrow("BlockMatrix::matrixCopy(m,x,y), x or y is out of range.");

  // Check dim
  if ((*map(row, col)).size(0) != m.size(0) || (*map(row, col)).size(1) != m.size(1))
    SiconosMatrixException::selfThrow("BlockMatrix::matrixCopy(m,x,y), block(x,y) of current matrix and m have inconsistent sizes.");
  *(map(row, col)) = m; // copy

  // Warning: this suppose that no block of the matrix can be a NULL pointer.
}

// OPERATORS

double BlockMatrix::getValue(unsigned int row, unsigned int col)
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= tabRow[nbRow] && nbRow < tabRow.size())
    nbRow ++;

  while (col >= tabCol[nbCol] && nbCol < tabCol.size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= tabRow[nbRow - 1];
  if (nbCol != 0)
    posCol -= tabCol[nbCol - 1];

  return (*map(nbRow, nbCol))(posRow, posCol);
}

void BlockMatrix::setValue(unsigned int row, unsigned int col, double value)
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= tabRow[nbRow] && nbRow < tabRow.size())
    nbRow ++;

  while (col >= tabCol[nbCol] && nbCol < tabCol.size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= tabRow[nbRow - 1];
  if (nbCol != 0)
    posCol -= tabCol[nbCol - 1];

  (*map(nbRow, nbCol))(posRow, posCol) = value;
}

double& BlockMatrix::operator()(unsigned int row, unsigned int col)
{
  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= tabRow[nbRow] && nbRow < tabRow.size())
    nbRow ++;

  while (col >= tabCol[nbCol] && nbCol < tabCol.size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= tabRow[nbRow - 1];
  if (nbCol != 0)
    posCol -= tabCol[nbCol - 1];

  return (*map(nbRow, nbCol))(posRow, posCol);
}

double BlockMatrix::operator()(unsigned int row, unsigned int col)const
{

  unsigned int nbRow = 0;
  unsigned int nbCol = 0;

  while (row >= tabRow[nbRow] && nbRow < tabRow.size())
    nbRow ++;

  while (col >= tabCol[nbCol] && nbCol < tabCol.size())
    nbCol ++;

  unsigned int posRow = row;
  unsigned int posCol = col;

  if (nbRow != 0)
    posRow -= tabRow[nbRow - 1];
  if (nbCol != 0)
    posCol -= tabCol[nbCol - 1];

  return (*map(nbRow, nbCol))(posRow, posCol);
}

BlockMatrix& BlockMatrix::operator = (const SiconosMatrix &m)
{
  if (&m == this) return *this;

  if (m.size(0) != size(0) || m.size(1) != size(1))
  {
    SiconosMatrixException::selfThrow("operator = (const SiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  // Warning: we do not reallocate the blocks, but only copy the values. This means that
  // all blocks are already allocated and that dim of m and mat are to be consistent.
  // Thus, tabRow and tabCol remains unchanged.
  if (m.isBlock())
  {
    unsigned int n1 = (static_cast<BlockMatrix>(m)).numberOfBlockInARow();
    unsigned int n2 = (static_cast<BlockMatrix>(m)).numberOfBlockInACol();
    if (tabRow.size() != n1 || tabCol.size() != n2)
      SiconosMatrixException::selfThrow("BlockMatrix operator = Left and Right blocks have inconsistent sizes.");
    unsigned int i = 0, j = 0;
    BlocksMat Mmap = m.getAllBlocks();
    if (STDMAP == 1)
    {
      //    BlocksMat::array_type::iterator it;
      //      unsigned int col = Mmap.size2();

      //    for(it=(Mmap.data ()).begin (); it!=(Mmap.data ()).end (); ++it)
      //      {
      //        i = (it->first)/col;
      //        j = (it->first) - i * col;
      //        *(map(i, j)) = *(it->second);
      //      }
    }
    else
    {
      BlocksMat::iterator1 it;
      BlocksMat::iterator2 it2;
      for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
      {
        for (it2 = it.begin(); it2 != it.end(); it2++)
        {
          i = it2.index1();
          j = it2.index2();
          *(map(i, j)) = **it2;
        }
      }
    }
  }
  else // if m is a SimpleMatrix
  {
    for (unsigned int i = 0; i < size(0); ++i)
      for (unsigned int j = 0; j < size(1); ++j)
        (*this)(i, j) = m(i, j);
  }

  return *this;
}

// BlockMatrix& BlockMatrix::operator = (const BlockMatrix &m){
//   if(&m == this) return *this;

//   if(m.size(0) != size(0)||m.size(1) != size(1)){
//     SiconosMatrixException::selfThrow("operator = (const BlockMatrix): Left and Right values have inconsistent sizes.");
//   }

//   return *this;
// }

BlockMatrix& BlockMatrix::operator *= (double d)
{
  //  unsigned int col = tabCol.size ();
  if (STDMAP == 1)
  {
    //     BlocksMat::array_type::iterator it;

    //     for(it=(map.data ()).begin (); it!=(map.data ()).end (); ++it)
    //       if( (it->second) != NULL)
    //  *(it->second) *= d;
    //       else
    //  SiconosMatrixException::selfThrow("BlockMatrix operator *=, a block is a NULL pointer.");
  }
  else
  {
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    unsigned int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if ((*it2) != NULL)
          (**it2) *= d;
        else
          SiconosMatrixException::selfThrow("BlockMatrix operator *=, a block is a NULL pointer.");
      }
    }
  }
  return *this;
}

BlockMatrix& BlockMatrix::operator *= (int d)
{
  if (STDMAP == 1)
  {
    //     BlocksMat::array_type::iterator it;

    //     for(it=(map.data ()).begin (); it!=(map.data ()).end (); ++it)
    //       {
    //  if( (it->second) != NULL)
    //    *(it->second) *= d;
    //  else
    //    SiconosMatrixException::selfThrow("BlockMatrix operator *=, a block is a NULL pointer.");
    //       }
  }
  else
  {
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    unsigned int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if ((*it2) != NULL)
          (**it2) *= d;
        else
          SiconosMatrixException::selfThrow("BlockMatrix operator *=, a block is a NULL pointer.");
      }
    }
  }
  return *this;
}

BlockMatrix& BlockMatrix::operator /= (double d)
{
  if (STDMAP == 1)
  {
    //     BlocksMat::array_type::iterator it;

    //     for(it=(map.data ()).begin (); it!=(map.data ()).end (); ++it){
    //       if( (it->second) != NULL)
    //  *(it->second) /= d;
    //       else
    //  SiconosMatrixException::selfThrow("BlockMatrix operator /=, a block is a NULL pointer.");
    //     }
  }
  else
  {
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    unsigned int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if ((*it2) != NULL)
          (**it2) /= d;
        else
          SiconosMatrixException::selfThrow("BlockMatrix operator /=, a block is a NULL pointer.");
      }
    }
  }
  return *this;
}

BlockMatrix& BlockMatrix::operator /= (int d)
{
  if (STDMAP == 1)
  {
    //     BlocksMat::array_type::iterator it;

    //     for(it=(map.data ()).begin (); it!=(map.data ()).end (); ++it){
    //       if( (it->second) != NULL)
    //  *(it->second) /= d;
    //       else
    //  SiconosMatrixException::selfThrow("BlockMatrix operator /=, a block is a NULL pointer.");
    //     }
  }
  else
  {
    BlocksMat::iterator1 it;
    BlocksMat::iterator2 it2;
    unsigned int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if ((*it2) != NULL)
          (**it2) /= d;
        else
          SiconosMatrixException::selfThrow("BlockMatrix operator /=, a block is a NULL pointer.");
      }
    }
  }
  return *this;
}

BlockMatrix& BlockMatrix::operator += (const SiconosMatrix &m)
{
  if (m.size(0) != size(0) || m.size(1) != size(1))
  {
    SiconosMatrixException::selfThrow("operator += (const SiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  if (m.isBlock())
  {
    BlocksMat Mmap = (dynamic_cast<const BlockMatrix&>(m)).map;

    if (STDMAP == 1)
    {
      //    unsigned int col = tabCol.size ();
      //       BlocksMat::array_type::iterator it;
      //       unsigned int i=0, j=0;


      //       for(it=(Mmap.data ()).begin (); it!=(Mmap.data ()).end (); ++it){
      //  i = (it->first)/col;
      //  j = (it->first) - i * col;
      //  *map(i, j) +=  *(it->second);
      //       }
    }
    else
    {
      BlocksMat::iterator1 it;
      BlocksMat::iterator2 it2;
      unsigned int i = 0, j = 0;
      for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
      {
        for (it2 = it.begin(); it2 != it.end(); it2++)
        {
          i = it2.index1();
          j = it2.index2();
          *map(i, j) += (**it2);
        }
      }
    }
  }

  else
  {
    //SUM USING OPERATOR ()
    for (unsigned int i = 0; i < size(0); i++)
    {

      for (unsigned int j = 0; j < size(1); j++)
      {
        (*this)(i, j) += m(i, j);
      }
    }
  }
  return *this;
}

BlockMatrix& BlockMatrix::operator -= (const SiconosMatrix &m)
{
  if (m.size(0) != size(0) || m.size(1) != size(1))
  {
    SiconosMatrixException::selfThrow("operator -= (const SiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  if (m.isBlock() == true)
  {
    BlocksMat Mmap = (dynamic_cast<const BlockMatrix&>(m)).map;

    if (STDMAP == 1)
    {
      //    unsigned int col = tabCol.size ();
      //       BlocksMat::array_type::iterator it;
      //       unsigned int i=0, j=0;

      //       for(it=(Mmap.data ()).begin (); it!=(Mmap.data ()).end (); ++it){
      //  i = (it->first)/col;
      //  j = (it->first) - i * col;
      //    *map(i, j) -=  *(it->second);
      //       }
    }
    else
    {
      BlocksMat::iterator1 it;
      BlocksMat::iterator2 it2;
      unsigned int i = 0, j = 0;
      for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
      {
        for (it2 = it.begin(); it2 != it.end(); it2++)
        {
          i = it2.index1();
          j = it2.index2();
          *map(i, j) -= (**it2);
        }
      }
    }
  }

  else
  {
    //SUBSTRACTION USING OPERATOR ()
    for (unsigned int i = 0; i < size(0); i++)
    {
      for (unsigned int j = 0; j < size(1); j++)
      {
        (*this)(i, j) -= m(i, j);
      }
    }
  }
  return *this;
}

void BlockMatrix::PLUFactorizationInPlace()
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUFactorizationInPlace: not yet implemented for Block Matrices.");
}

void BlockMatrix::PLUInverseInPlace()
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUInverseInPlace: not yet implemented for Block Matrices.");
}

void BlockMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUForwardBackwardInPlace: not yet implemented for Block Matrices.");
}

void BlockMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  SiconosMatrixException::selfThrow(" BlockMatrix::PLUForwardBackwardInPlace: not yet implemented for Block Matrices.");
}
