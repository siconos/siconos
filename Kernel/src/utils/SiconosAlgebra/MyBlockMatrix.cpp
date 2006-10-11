#include "MyBlockMatrix.h"

//#define STDMAP 0

void MyBlockMatrix::addInTab(unsigned int row, unsigned int col)
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

void MyBlockMatrix::makeTab(unsigned int row, unsigned int col)
{

  unsigned int dim;
  tabRow.clear();
  tabCol.clear();

  if (STDMAP == 1)
  {
    // Col sweeping
    for (unsigned int j = 0; j < map.size2(); j++)
    {
      tabCol.push_back((*map(0, j)).size2());
      dim = tabCol.size();
      if (dim > 1)
        tabCol[dim - 1] += tabCol[dim - 2];
    }
    // Row sweeping
    for (unsigned int i = 0; i < map.size1(); i++)
    {
      tabRow.push_back((*map(i, 0)).size1());
      dim = tabRow.size();
      if (dim > 1)
        tabRow[dim - 1] += tabRow[dim - 2];
    }
  }
  else
  {
    //mapped Mmap = map;
    mapped::iterator1 it;
    mapped::iterator2 it2;

    // Col sweeping
    for (it2 = map.begin2(); it2 != map.end2(); ++it2)
    {
      tabCol.push_back((**it2).size2());
      dim = tabCol.size();
      if (dim > 1)
        tabCol[dim - 1] += tabCol[dim - 2];
    }

    // Row sweeping
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      tabRow.push_back((**it).size1());
      dim = tabRow.size();
      if (dim > 1)
        tabRow[dim - 1] += tabRow[dim - 2];
    }
  }
}

MyBlockMatrix::MyBlockMatrix(const MySiconosMatrix &m): MySiconosMatrix(true)
{
  if (!m.isBlock())
    SiconosMatrixException::selfThrow("BlockMatrix copy constructor from a SimpleMatrix: forbidden operation.");

  // mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;
  mapped Mmap = m.getMap();
  unsigned int col = Mmap.size2();
  unsigned int row = Mmap.size1();
  isBlockAllocatedIn.resize(col * row);

  map.resize(row, col, false);
  unsigned int i = 0, j = 0;
  // STDMAP = 1 means that we use the map::iterator of standard library, else we use the iterator of boost map

  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
    {
      // Computation of index row and index column corresponding
      i = (it->first) / col;
      j = (it->first) - i * col;
      map(i, j) = new MySimpleMatrix(*(it->second));
      isBlockAllocatedIn[(it->first) ] = true;
    }
  }

  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;

    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2 ++)
      {
        // Computation of index row and index column corresponding
        i = it2.index1();
        j = it2.index2();
        map(i, j) = new MySimpleMatrix(**it2);
        isBlockAllocatedIn[i * col + j] = true;
      }
    }
  }
  makeTab(row, col);
}

MyBlockMatrix::MyBlockMatrix(const MyBlockMatrix &m): MySiconosMatrix(true)
{
  // mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;
  mapped Mmap = m.getMap();
  unsigned int col = Mmap.size2();
  unsigned int row = Mmap.size1();
  isBlockAllocatedIn.resize(col * row);

  map.resize(row, col, false);
  unsigned int i = 0, j = 0;
  // STDMAP = 1 means that we use the map::iterator of standard library, else we use the iterator of boost map

  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
    {
      // Computation of index row and index column corresponding
      i = (it->first) / col;
      j = (it->first) - i * col;
      map(i, j) = new MySimpleMatrix(*(it->second));
      isBlockAllocatedIn[(it->first) ] = true;
    }
  }

  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;

    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2 ++)
      {
        // Computation of index row and index column corresponding
        i = it2.index1();
        j = it2.index2();
        map(i, j) = new MySimpleMatrix(**it2);
        isBlockAllocatedIn[i * col + j] = true;
      }
    }
  }
  makeTab(row, col);
}

MyBlockMatrix::MyBlockMatrix(mapped& Mmap): MySiconosMatrix(true)
{

  unsigned int i = 0, j = 0;
  unsigned int row = Mmap.size1();
  unsigned int col = Mmap.size2();
  isBlockAllocatedIn.resize(row * col);
  map.resize(row, col, false);

  // STDMAP = 1 means that we use the map::iterator of standard library, else we use the iterator of boost map
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it1;

    for (it1 = (Mmap.data()).begin(); it1 != (Mmap.data()).end(); ++it1)
    {
      // Computation of index row and index column corresponding
      i = (it1->first) / col;
      j = (it1->first) - i * col;
      map(i, j) = new MySimpleMatrix(*(it1->second));
      isBlockAllocatedIn[(it1->first) ] = true;
    }
  }
  else
  {

    mapped::iterator1 it;
    mapped::iterator2 it2;

    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2 ++)
      {
        // Computation of index row and index column corresponding
        i = it2.index1();
        j = it2.index2();
        map(i, j) = new MySimpleMatrix(**it2);
        isBlockAllocatedIn[i * col + j] = true;
      }
    }
  }
  makeTab(row, col);
}

MyBlockMatrix::MyBlockMatrix(const std::vector<MySiconosMatrix* > &m, unsigned int row, unsigned int col): MySiconosMatrix(true)
{
  if (m.size() != (row * col))
    SiconosMatrixException::selfThrow("BlockMatrix constructor from a vector<MySiconosMatrix*>, number of blocks inconsistent with provided dimensions.");

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

MyBlockMatrix::MyBlockMatrix(MySiconosMatrix* A, MySiconosMatrix* B, MySiconosMatrix* C, MySiconosMatrix* D): MySiconosMatrix(true)
{
  if (A->size1() != B->size1() || C->size1() != D->size1() ||  A->size2() != C->size2() ||  B->size2() != D->size2())
    SiconosMatrixException::selfThrow("BlockMatrix constructor(A,B,C,D), inconsistent sizes between A, B, C or D MySiconosMatrices.");

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

MyBlockMatrix::~MyBlockMatrix(void)
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
unsigned int MyBlockMatrix::size1(void)const
{
  return tabRow[ tabRow.size() - 1 ];
}

// return the number of columns of blocks
unsigned int MyBlockMatrix::size2(void)const
{
  return tabCol[ tabCol.size() - 1 ];
}

void MyBlockMatrix::resize(unsigned int row, unsigned int col, unsigned int lower, unsigned int upper, bool preserve)
{
  //   if(lower != 0 || upper != 0)
  //     SiconosMatrixException::selfThrow("MyBlockMatrix::resize : lower or upper not equal to 0.0");
  //   map.resize(row, col, preserve);
  SiconosMatrixException::selfThrow("MyBlockMatrix::resize : forbidden for block matrices.");
}

const double MyBlockMatrix::normInf(void)const
{
  double sum = 0, norm = 0;
  for (unsigned int i = 0; i < size1(); i++)
  {
    for (unsigned int j = 0; j < size2(); j++)
    {
      sum += (*this)(i, j);
    }
    if (sum > norm) norm = sum;
    sum = 0;
  }
  return norm;
}

void MyBlockMatrix::zero(void)
{
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;
    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      (it->second)->zero();
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        (**it2).zero();
      }
    }
  }
}

void MyBlockMatrix::eye(void)
{
  if (STDMAP == 1)
  {
    unsigned int col = map.size2();
    unsigned int p = 0, q = 0;
    mapped::array_type::iterator it;
    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      p = (it->first) / col;
      q = (it->first) - p * col;
      if (p == q)
      {
        (it->second)->eye();
      }
      else
      {
        (it->second)->zero();
      }
    }
  }
  else
  {
    mapped Mmap = map;
    mapped::iterator1 it;
    mapped::iterator2 it2;

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

void MyBlockMatrix::display(void)const
{
  if (STDMAP == 1)
  {
    mapped::array_type Mmap = map.data();
    mapped::array_type::iterator it;
    for (it = Mmap.begin(); it != Mmap.end(); ++it)
    {
      (it->second)->display();
    }
  }
  else
  {
    mapped Mmap = map;
    mapped::iterator1 it;
    mapped::iterator2 it2;
    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        (**it2).display();
      }
    }
  }

}

void MyBlockMatrix::getBlock(unsigned int row, unsigned int col, MySiconosMatrix &m)const
{
  m = dynamic_cast<MySiconosMatrix&>(*(map(row, col)));
}

MySiconosMatrix* MyBlockMatrix::getBlockPtr(unsigned int row, unsigned int col)
{
  return map(row, col);
}

const std::deque<bool> MyBlockMatrix::getBlockAllocated(void)const
{
  return isBlockAllocatedIn;
}

unsigned int MyBlockMatrix::getNum(void)const
{
  SiconosMatrixException::selfThrow("getNum of a block is forbidden.");
  return 0;
}

void MyBlockMatrix::getRow(unsigned int r, MySimpleVector &v) const
{
  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  if (r > size1())
    SiconosMatrixException::selfThrow("MyBlockMatrix:getRow : row number is out of range");

  // Verification of the size of the result vector
  if (v.size() != size2())
    SiconosMatrixException::selfThrow("MyBlockMatrix:getRow : inconsistent sizes");

  // Find the row-block number where "r" is
  while (r >= tabRow[numRow] && numRow < tabRow.size())
    numRow ++;

  // Computation of the value of the index row into this block
  if (numRow != 0)
    posRow -= tabRow[numRow - 1];

  for (unsigned int j = 0; j < tabCol.size(); j++)
  {
    start = stop;
    stop += (*map(numRow, j)).size2();
    subrange(*(v.getDensePtr()), start, stop) = row((*map(numRow, j)).getDense(), posRow);
  }
}

void MyBlockMatrix::getCol(unsigned int c, MySimpleVector &v) const
{
  unsigned int numCol = 0, posCol = c, start = 0, stop = 0;

  if (c > size2())
    SiconosMatrixException::selfThrow("MyBlockMatrix:getCol : column number is out of range");

  // Verification of the size of the result vector
  if (v.size() != size1())
    SiconosMatrixException::selfThrow("MyBlockMatrix:getcol : inconsistent sizes");

  // Find the column-block number where "c" is
  while (c >= tabCol[numCol] && numCol < tabCol.size())
    numCol ++;

  // Computation of the value of the index column into this block
  if (numCol != 0)
    posCol -= tabCol[numCol - 1];

  for (unsigned int i = 0; i < tabRow.size(); i++)
  {
    start = stop;
    stop += (*map(i, numCol)).size1();
    subrange(*(v.getDensePtr()), start, stop) = column((*map(i, numCol)).getDense(), posCol);
  }
}

void MyBlockMatrix::setRow(unsigned int r, const MySimpleVector &v)
{

  unsigned int numRow = 0, posRow = r, start = 0, stop = 0;

  unsigned int nbCol = size2();
  assert(v.size() == nbCol);

  while (r >= tabRow[numRow] && numRow < tabRow.size())
    numRow ++;

  if (numRow != 0)
    posRow -= tabRow[numRow - 1];

  for (unsigned int j = 0; j < tabCol.size(); j++)
  {
    start = stop;
    stop += (*map(numRow, j)).size2();
    row(*((*map(numRow, j)).getDensePtr()), posRow) = subrange(*(v.getDensePtr()), start, stop);
  }
}

void MyBlockMatrix::setCol(unsigned int col, const MySimpleVector &v)
{

  unsigned int numCol = 0, posCol = col, start = 0, stop = 0;

  unsigned int nbRow = size1();
  assert(v.size() == nbRow);

  while (col >= tabCol[numCol] && numCol < tabCol.size())
    numCol ++;

  if (numCol != 0)
    posCol -= tabCol[numCol - 1];

  for (unsigned int i = 0; i < tabRow.size(); i++)
  {
    start = stop;
    stop += (*map(i, numCol)).size1();
    column(*((*map(i, numCol)).getDensePtr()), posCol) = subrange(*(v.getDensePtr()), start, stop);
  }
}

// return the boost dense matrix of the block (i, j)
const DenseMat  MyBlockMatrix::getDense(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 1);
  SiconosMatrixException::selfThrow("DenseMat getDense(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Dense matrix");

  return (*map(row, col)).getDense();
}

// return the boost triangular matrix of the block (i, j)
const TriangMat MyBlockMatrix::getTriang(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 2);
  SiconosMatrixException::selfThrow("TriangMat getTriang(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Triangular matrix");
  return (*map(row, col)).getTriang();
}

// return the boost symmetric matrix of the block (i, j)
const SymMat MyBlockMatrix::getSym(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 3);
  SiconosMatrixException::selfThrow("SymMat getSym(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Symmmetric matrix");
  return (*map(row, col)).getSym();
}

// return the boost sparse matrix of the block (i, j)
const SparseMat  MyBlockMatrix::getSparse(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 4);
  SiconosMatrixException::selfThrow("SparseMat getSparse(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Sparse matrix");

  return (*map(row, col)).getSparse();
}

// return the boost banded matrix of the block (i, j)
const BandedMat  MyBlockMatrix::getBanded(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("BandedMat getBanded(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Banded matrix");

  return (*map(row, col)).getBanded();
}

// The following functions return the corresponding pointers
DenseMat*  MyBlockMatrix::getDensePtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 1);
  SiconosMatrixException::selfThrow("DenseMat* getDensePtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Dense matrix");

  return (*map(row, col)).getDensePtr();
}

TriangMat* MyBlockMatrix::getTriangPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 2);
  SiconosMatrixException::selfThrow("TriangMat* getTriangPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Triangular matrix");

  return (*map(row, col)).getTriangPtr();
}
SymMat* MyBlockMatrix::getSymPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 3);
  SiconosMatrixException::selfThrow("SymMat* getSymPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Symmmetric matrix");
  return (*map(row, col)).getSymPtr();
}

SparseMat*  MyBlockMatrix::getSparsePtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 4);
  SiconosMatrixException::selfThrow("SparseMat* getSparsePtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Sparse matrix");

  return (*map(row, col)).getSparsePtr();
}

BandedMat*  MyBlockMatrix::getBandedPtr(unsigned int row, unsigned int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("BandedMat* getBandedPtr(unsigned int row, unsigned int col) : the matrix at (row, col) is not a Banded matrix");

  return (*map(row, col)).getBandedPtr();
}

// return the boost mapped matrix of a block matrix
const mapped MyBlockMatrix::getMap(void)const
{
  return map;
}

void MyBlockMatrix::matrixCopy(const MySiconosMatrix &m, unsigned int row, unsigned int col)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("MyBlockMatrix::matrixCopy of a block into an other block is forbidden.");

  if (row > tabRow.size() || col > tabCol.size())
    SiconosMatrixException::selfThrow("MyBlockMatrix::matrixCopy(m,x,y), x or y is out of range.");

  // Check dim
  if ((*map(row, col)).size1() != m.size1() || (*map(row, col)).size2() != m.size2())
    SiconosMatrixException::selfThrow("MyBlockMatrix::matrixCopy(m,x,y), block(x,y) of current matrix and m have inconsistent sizes.");
  *(map(row, col)) = m; // copy

  // Warning: this suppose that no block of the matrix can be a NULL pointer.
}

// OPERATORS

double& MyBlockMatrix::operator()(unsigned int row, unsigned int col)
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

double MyBlockMatrix::operator()(unsigned int row, unsigned int col)const
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

MyBlockMatrix& MyBlockMatrix::operator = (const MySiconosMatrix &m)
{
  if (&m == this) return *this;

  if (m.size1() != size1() || m.size2() != size2())
  {
    SiconosMatrixException::selfThrow("operator = (const MySiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  // Warning: we do not reallocate the blocks, but only copy the values. This means that
  // all blocks are already allocated and that dim of m and mat are to be consistent.
  // Thus, tabRow and tabCol remains unchanged.
  if (m.isBlock())
  {
    unsigned int n1 = (static_cast<MyBlockMatrix>(m)).numberOfBlockInARow();
    unsigned int n2 = (static_cast<MyBlockMatrix>(m)).numberOfBlockInACol();
    if (tabRow.size() != n1 || tabCol.size() != n2)
      SiconosMatrixException::selfThrow("BlockMatrix operator = Left and Right blocks have inconsistent sizes.");
    unsigned int i = 0, j = 0;
    mapped Mmap = m.getMap();
    unsigned int col = Mmap.size2();
    if (STDMAP == 1)
    {
      mapped::array_type::iterator it;

      for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
      {
        i = (it->first) / col;
        j = (it->first) - i * col;
        *(map(i, j)) = *(it->second);
      }
    }
    else
    {
      mapped::iterator1 it;
      mapped::iterator2 it2;
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
  else // if m is a MySimpleMatrix
  {
    for (unsigned int i = 0; i < size1(); ++i)
      for (unsigned int j = 0; j < size2(); ++j)
        (*this)(i, j) = m(i, j);
  }

  return *this;
}

// MyBlockMatrix& MyBlockMatrix::operator = (const MyBlockMatrix &m){
//   if(&m == this) return *this;

//   if(m.size1 () != size1 ()||m.size2 () != size2 ()){
//     SiconosMatrixException::selfThrow("operator = (const MyBlockMatrix): Left and Right values have inconsistent sizes.");
//   }

//   return *this;
// }

MyBlockMatrix& MyBlockMatrix::operator *= (double d)
{
  //  unsigned int col = tabCol.size ();
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
      if ((it->second) != NULL)
        *(it->second) *= d;
      else
        SiconosMatrixException::selfThrow("BlockMatrix operator *=, a block is a NULL pointer.");
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
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

MyBlockMatrix& MyBlockMatrix::operator *= (int d)
{
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      if ((it->second) != NULL)
        *(it->second) *= d;
      else
        SiconosMatrixException::selfThrow("BlockMatrix operator *=, a block is a NULL pointer.");
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
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

MyBlockMatrix& MyBlockMatrix::operator /= (double d)
{
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      if ((it->second) != NULL)
        *(it->second) /= d;
      else
        SiconosMatrixException::selfThrow("BlockMatrix operator /=, a block is a NULL pointer.");
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
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

MyBlockMatrix& MyBlockMatrix::operator /= (int d)
{
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      if ((it->second) != NULL)
        *(it->second) /= d;
      else
        SiconosMatrixException::selfThrow("BlockMatrix operator /=, a block is a NULL pointer.");
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
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

MyBlockMatrix& MyBlockMatrix::operator += (const MySiconosMatrix &m)
{
  if (m.size1() != size1() || m.size2() != size2())
  {
    SiconosMatrixException::selfThrow("operator += (const MySiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  if (m.isBlock())
  {
    mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;
    unsigned int col = tabCol.size();

    if (STDMAP == 1)
    {
      mapped::array_type::iterator it;
      unsigned int i = 0, j = 0;


      for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
      {
        i = (it->first) / col;
        j = (it->first) - i * col;
        *map(i, j) +=  *(it->second);
      }
    }
    else
    {
      mapped::iterator1 it;
      mapped::iterator2 it2;
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
    for (unsigned int i = 0; i < size1(); i++)
    {

      for (unsigned int j = 0; j < size2(); j++)
      {
        (*this)(i, j) += m(i, j);
      }
    }
  }
  return *this;
}

MyBlockMatrix& MyBlockMatrix::operator -= (const MySiconosMatrix &m)
{
  if (m.size1() != size1() || m.size2() != size2())
  {
    SiconosMatrixException::selfThrow("operator -= (const MySiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  if (m.isBlock() == true)
  {
    mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;
    unsigned int col = tabCol.size();

    if (STDMAP == 1)
    {
      mapped::array_type::iterator it;
      unsigned int i = 0, j = 0;

      for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
      {
        i = (it->first) / col;
        j = (it->first) - i * col;
        *map(i, j) -=  *(it->second);
      }
    }
    else
    {
      mapped::iterator1 it;
      mapped::iterator2 it2;
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
    for (unsigned int i = 0; i < size1(); i++)
    {
      for (unsigned int j = 0; j < size2(); j++)
      {
        (*this)(i, j) -= m(i, j);
      }
    }
  }
  return *this;
}

void MyBlockMatrix::PLUFactorizationInPlace()
{
  SiconosMatrixException::selfThrow(" MyBlockMatrix::PLUFactorizationInPlace: not yet implemented for Block Matrices.");
}

void MyBlockMatrix::PLUInverseInPlace()
{
  SiconosMatrixException::selfThrow(" MyBlockMatrix::PLUInverseInPlace: not yet implemented for Block Matrices.");
}

void MyBlockMatrix::PLUForwardBackwardInPlace(MySiconosMatrix &B)
{
  SiconosMatrixException::selfThrow(" MyBlockMatrix::PLUForwardBackwardInPlace: not yet implemented for Block Matrices.");
}

void MyBlockMatrix::PLUForwardBackwardInPlace(MySiconosVector &B)
{
  SiconosMatrixException::selfThrow(" MyBlockMatrix::PLUForwardBackwardInPlace: not yet implemented for Block Matrices.");
}
