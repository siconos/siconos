#include "MyBlockMatrix.h"

//#define STDMAP 0

void MyBlockMatrix::addInTab(int row, int col)
{

  tabRow.push_back(row);
  tabCol.push_back(col);

  int nbRow = tabRow.size();
  int nbCol = tabCol.size();

  if (nbRow > 1)
    tabRow[nbRow - 1] += tabRow[nbRow - 2];
  if (nbCol > 1)
    tabCol[nbCol - 1] += tabCol[nbCol - 2];
}

void MyBlockMatrix::makeTab(int row, int col)
{

  int dim;
  tabRow.clear();
  tabCol.clear();

  if (STDMAP == 1)
  {
    // Col sweeping
    for (int j = 0; j < map.size2(); j++)
    {
      tabCol.push_back((*map(0, j)).size2());
      dim = tabCol.size();
      if (dim > 1)
        tabCol[dim - 1] += tabCol[dim - 2];
    }
    // Row sweeping
    for (int i = 0; i < map.size1(); i++)
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



MyBlockMatrix::MyBlockMatrix(const MySiconosMatrix &m)
{
  if (m.isBlock() == false)
    SiconosMatrixException::selfThrow("BlockMatrix copy constructor from a SimpleMatrix: forbidden operation.");

  setIsBlock(true);
  isBlockAllocatedIn.resize(m.size1() * m.size2());
  mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;

  tabRow = (dynamic_cast<const MyBlockMatrix&>(m)).tabRow;
  tabCol = (dynamic_cast<const MyBlockMatrix&>(m)).tabCol;

  int col = Mmap.size2();
  int row = Mmap.size1();
  map.resize(row, col, false);
  int i = 0, j = 0;
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
}

MyBlockMatrix::MyBlockMatrix(const MyBlockMatrix &m)
{

  setIsBlock(true);
  isBlockAllocatedIn.resize(m.size1() * m.size2());
  mapped Mmap = m.map;

  tabRow = (dynamic_cast<const MyBlockMatrix&>(m)).tabRow;
  tabCol = (dynamic_cast<const MyBlockMatrix&>(m)).tabCol;

  int col = Mmap.size2();
  int row = Mmap.size1();
  map.resize(row, col, false);
  int i = 0, j = 0;

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
}

MyBlockMatrix::MyBlockMatrix(mapped& Mmap)
{

  int i = 0, j = 0;
  int row = Mmap.size1();
  int col = Mmap.size2();
  setIsBlock(true);
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

MyBlockMatrix::MyBlockMatrix(const std::vector<MySiconosMatrix* > &m, int row, int col)
{
  if (m.size() != (row * col))
    SiconosMatrixException::selfThrow("BlockMatrix constructor from a vector<LaGenMatDouble>, number of blocks inconsistent with provided dimensions.");

  setIsBlock(true);
  isBlockAllocatedIn.resize(row * col);
  map.resize(row, col, false);
  tabRow.reserve(row);
  tabCol.reserve(col);

  for (int i = 0; i < row; i++)
  {
    for (int j = 0; j < col; j++)
    {
      map(i, j) = new MySimpleMatrix(*m [i * col + j]);
      isBlockAllocatedIn[i * col + j] == true;
      //isBlockAllocatedIn.push_back (false);
    }
  }
  makeTab(row, col);
}

MyBlockMatrix::~MyBlockMatrix(void)
{
  int row = map.size1();
  int col = map.size2();
  // we call delete on the map only if the pointer was allocated
  for (unsigned i = 0; i < row; ++i)
  {
    for (unsigned j = 0; j < col; ++j)
    {
      if (isBlockAllocatedIn[i * col + j] == true)
        delete(map(i, j));
    }
  }
}


void MyBlockMatrix::getBlock(int row, int col, MySiconosMatrix &m)const
{
  m = dynamic_cast<MySiconosMatrix&>(*(map(row, col)));
}

const std::deque<bool> MyBlockMatrix::getBlockAllocated(void)const
{
  return isBlockAllocatedIn;
}

const double MyBlockMatrix::normInf(void)const
{
  double sum = 0;
  double norm = 0;
  for (int i = 0; i < size1(); i++)
  {
    for (int j = 0; j < size2(); j++)
    {
      sum += (*this)(i, j);
    }
    if (sum > norm) norm = sum;
    sum = 0;
  }
  return norm;
}

void MyBlockMatrix::blockMatrixCopy(const MySiconosMatrix &m, int row, int col)
{

  if (m.isBlock() == true)
    SiconosMatrixException::selfThrow("blockMatrixCopy of a block to an other block is forbidden.");


  int nbRow = 0;
  int nbCol = 0;

  // Computation of the index row and the index column of the block corresponding
  while (row >= tabRow[nbRow] && nbRow < tabRow.size())
    nbRow ++;

  while (col >= tabCol[nbCol] && nbCol < tabCol.size())
    nbCol ++;
  // we delete the block if it was allocated to avoid leak memory
  if (isBlockAllocatedIn[nbRow * tabCol.size() + nbCol] == true)
    delete map(nbRow, nbCol);
  map(nbRow, nbCol) = new MySimpleMatrix(m);
  isBlockAllocatedIn[nbRow * tabCol.size() + nbCol] = true;
}

int MyBlockMatrix::getNum(void)const
{
  SiconosMatrixException::selfThrow("getNum of a block is forbidden.");
}

void MyBlockMatrix::setNum(int)
{
  SiconosMatrixException::selfThrow("setNum of a block is forbidden.");
}


void MyBlockMatrix::getRow(int row, MySimpleVector &v)const
{

  int numRow = 0, posRow = row, start = 0, stop = 0, step = 0;

  if (row > size1() || row < 0)
    SiconosMatrixException::selfThrow("getRow : row is out of range");

  // Verification of the size of the result vector
  if (v.size() != size2())
    SiconosMatrixException::selfThrow("getRow : inconsistent sizes");

  int nbCol = size2();

  // Determination of the row of the block corresponding

  while (row >= tabRow[numRow] && numRow < tabRow.size())
    numRow ++;
  // Computation of the value of the index row into this block
  if (numRow != 0)
    posRow -= tabRow[numRow - 1];

  // loop for copying the values into the vector
  DenseVect vect(size2());

  for (int j = 0; j < tabCol.size(); j++)
  {
    step = (*map(numRow, j)).size2();
    start = stop;
    stop += step;
    MySimpleVector tmp(DENSE, step);
    (*map(numRow, j)).getRow(posRow, tmp);

    subrange(vect, start, stop) = tmp.getDense();
  }
  MySimpleVector p(vect);
  v = p;
}

void MyBlockMatrix::getCol(int col, MySimpleVector &v)const
{

  int numCol = 0, posCol = col, start = 0, stop = 0, step = 0;

  if (col > size2() || col < 0)
    SiconosMatrixException::selfThrow("getCol : col is out of range");

  // Verification of the size of the result vector
  if (v.size() != size1())
    SiconosMatrixException::selfThrow("getcol : inconsistent sizes");

  int nbRow = size1();

  // Determination of the column of the block corresponding
  while (col >= tabCol[numCol] && numCol < tabCol.size())
    numCol ++;
  // Computation of the value of the index column into this block
  if (numCol != 0)
    posCol -= tabCol[numCol - 1];

  DenseVect vect(size1());

  // loop for copying the values into the vector
  for (int i = 0; i < tabRow.size(); i++)
  {
    step = (*map(i, numCol)).size1();
    start = stop;
    stop += step;
    MySimpleVector tmp(DENSE, step);
    (*map(i, numCol)).getCol(posCol, tmp);

    subrange(vect, start, stop) = tmp.getDense();
  }
  MySimpleVector p(vect);
  v = p;
}

void MyBlockMatrix::setRow(int row, const MySimpleVector &v)
{

  int numRow = 0, posRow = row, start = 0, stop = 0, step = 0;

  int nbCol = size2();
  assert(v.size() == nbCol);

  while (row >= tabRow[numRow] && numRow < tabRow.size())
    numRow ++;

  if (numRow != 0)
    posRow -= tabRow[numRow - 1];

  DenseVect vect(size2());
  vect = v.getDense();

  for (int j = 0; j < tabCol.size(); j++)
  {
    step = (*map(numRow, j)).size2();
    start = stop;
    stop += step;
    DenseVect tmp(step);
    tmp = subrange(vect, start, stop);
    MySimpleVector p(tmp);
    (*map(numRow, j)).setRow(posRow, p);

  }
}

void MyBlockMatrix::setCol(int col, const MySimpleVector &v)
{

  int numCol = 0, posCol = col, start = 0, stop = 0, step = 0;

  int nbRow = size1();
  assert(v.size() == nbRow);

  while (col >= tabCol[numCol] && numCol < tabCol.size())
    numCol ++;

  if (numCol != 0)
    posCol -= tabCol[numCol - 1];

  DenseVect vect(size1());
  vect = v.getDense();
  for (int i = 0; i < tabRow.size(); i++)
  {
    step = (*map(i, numCol)).size1();
    start = stop;
    stop += step;
    DenseVect tmp(step);
    tmp = subrange(vect, start, stop);
    MySimpleVector p(tmp);
    (*map(i, numCol)).setCol(posCol, p);

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
    int col = map.size2();
    int p = 0, q = 0;
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

// return the boost dense matrix of the block (i, j)
const DenseMat  MyBlockMatrix::getDense(int row, int col)const
{

  if ((*map(row, col)).getNum() != 1);
  SiconosMatrixException::selfThrow("DenseMat getDense(int row, int col) : the matrix at (row, col) is not a Dense matrix");

  return (*map(row, col)).getDense();
}

// return the boost triangular matrix of the block (i, j)
const TriangMat MyBlockMatrix::getTriang(int row, int col)const
{

  if ((*map(row, col)).getNum() != 2);
  SiconosMatrixException::selfThrow("TriangMat getTriang(int row, int col) : the matrix at (row, col) is not a Triangular matrix");
  return (*map(row, col)).getTriang();
}

// return the boost symmetric matrix of the block (i, j)
const SymMat MyBlockMatrix::getSym(int row, int col)const
{

  if ((*map(row, col)).getNum() != 3);
  SiconosMatrixException::selfThrow("SymMat getSym(int row, int col) : the matrix at (row, col) is not a Symmmetric matrix");
  return (*map(row, col)).getSym();
}

// return the boost sparse matrix of the block (i, j)
const SparseMat  MyBlockMatrix::getSparse(int row, int col)const
{

  if ((*map(row, col)).getNum() != 4);
  SiconosMatrixException::selfThrow("SparseMat getSparse(int row, int col) : the matrix at (row, col) is not a Sparse matrix");

  return (*map(row, col)).getSparse();
}

// return the boost banded matrix of the block (i, j)
const BandedMat  MyBlockMatrix::getBanded(int row, int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("BandedMat getBanded(int row, int col) : the matrix at (row, col) is not a Banded matrix");

  return (*map(row, col)).getBanded();
}

// The following functions return the corresponding pointers
const DenseMat*  MyBlockMatrix::getDensePtr(int row, int col)const
{

  if ((*map(row, col)).getNum() != 1);
  SiconosMatrixException::selfThrow("DenseMat* getDensePtr(int row, int col) : the matrix at (row, col) is not a Dense matrix");

  return (*map(row, col)).getDensePtr();
}

const TriangMat* MyBlockMatrix::getTriangPtr(int row, int col)const
{

  if ((*map(row, col)).getNum() != 2);
  SiconosMatrixException::selfThrow("TriangMat* getTriangPtr(int row, int col) : the matrix at (row, col) is not a Triangular matrix");

  return (*map(row, col)).getTriangPtr();
}
const SymMat* MyBlockMatrix::getSymPtr(int row, int col)const
{

  if ((*map(row, col)).getNum() != 3);
  SiconosMatrixException::selfThrow("SymMat* getSymPtr(int row, int col) : the matrix at (row, col) is not a Symmmetric matrix");
  return (*map(row, col)).getSymPtr();
}

const SparseMat*  MyBlockMatrix::getSparsePtr(int row, int col)const
{

  if ((*map(row, col)).getNum() != 4);
  SiconosMatrixException::selfThrow("SparseMat* getSparsePtr(int row, int col) : the matrix at (row, col) is not a Sparse matrix");

  return (*map(row, col)).getSparsePtr();
}

const BandedMat*  MyBlockMatrix::getBandedPtr(int row, int col)const
{

  if ((*map(row, col)).getNum() != 5);
  SiconosMatrixException::selfThrow("BandedMat* getBandedPtr(int row, int col) : the matrix at (row, col) is not a Banded matrix");

  return (*map(row, col)).getBandedPtr();
}

// return the boost mapped matrix of a block matrix
const mapped MyBlockMatrix::getMap(void)const
{
  return map;
}

// return the number of rows of blocks
int MyBlockMatrix::size1(void)const
{
  return tabRow[ tabRow.size() - 1 ];
}

// return the number of columns of blocks
int MyBlockMatrix::size2(void)const
{
  return tabCol[ tabCol.size() - 1 ];
}


void MyBlockMatrix::resize(int row, int col, int lower, int upper, bool preserve)
{
  if (lower != 0 || upper != 0)
    SiconosMatrixException::selfThrow("MyBlockMatrix::resize : lower or upper not equal to 0.0");
  map.resize(row, col, preserve);
}

// OPERATORS

double& MyBlockMatrix::operator()(int row, int col)
{

  int nbRow = 0;
  int nbCol = 0;

  while (row >= tabRow[nbRow] && nbRow < tabRow.size())
    nbRow ++;

  while (col >= tabCol[nbCol] && nbCol < tabCol.size())
    nbCol ++;

  int posRow = row;
  int posCol = col;

  if (nbRow != 0)
    posRow -= tabRow[nbRow - 1];
  if (nbCol != 0)
    posCol -= tabCol[nbCol - 1];

  return (*map(nbRow, nbCol))(posRow, posCol);
}

double MyBlockMatrix::operator()(int row, int col)const
{

  int nbRow = 0;
  int nbCol = 0;

  while (row >= tabRow[nbRow] && nbRow < tabRow.size())
    nbRow ++;

  while (col >= tabCol[nbCol] && nbCol < tabCol.size())
    nbCol ++;

  int posRow = row;
  int posCol = col;

  if (nbRow != 0)
    posRow -= tabRow[nbRow - 1];
  if (nbCol != 0)
    posCol -= tabCol[nbCol - 1];

  return (*map(nbRow, nbCol))(posRow, posCol);
}

const MyBlockMatrix& MyBlockMatrix::operator = (const MySiconosMatrix &m)
{
  if (&m == this) return *this;

  if (m.size1() != size1() || m.size2() != size2())
  {
    SiconosMatrixException::selfThrow("operator = (const MySiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  if (m.isBlock() == true)
  {
    mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;
    tabRow = (dynamic_cast<const MyBlockMatrix&>(m)).tabRow;
    tabCol = (dynamic_cast<const MyBlockMatrix&>(m)).tabCol;

    int col = tabCol.size();
    int i = 0, j = 0;
    if (STDMAP == 1)
    {
      mapped::array_type::iterator it;

      for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
      {
        i = (it->first) / col;
        j = (it->first) - i * col;
        if (isBlockAllocatedIn[i * col + j] == true)
        {
          delete(map(i, j));
        }
        map(i, j) = new MySimpleMatrix(*(it->second));
        isBlockAllocatedIn[i * col + j] = true;
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
          if (isBlockAllocatedIn[i * col + j] == true)
          {
            delete(map(i, j));
          }
          map(i, j) = new MySimpleMatrix(**it2);
          isBlockAllocatedIn[i * col + j] = true;
        }
      }
    }
    return *this;
  }

  else
  {
    SiconosMatrixException::selfThrow("operator = : BlockMatrix=MySimpleMatrix is a forbidden operation");
  }
}

const MyBlockMatrix& MyBlockMatrix::operator = (const MyBlockMatrix &m)
{
  if (&m == this) return *this;

  if (m.size1() != size1() || m.size2() != size2())
  {
    SiconosMatrixException::selfThrow("operator = (const MyBlockMatrix): Left and Right values have inconsistent sizes.");
  }

  mapped Mmap = m.map;
  tabRow = (dynamic_cast<const MyBlockMatrix&>(m)).tabRow;
  tabCol = (dynamic_cast<const MyBlockMatrix&>(m)).tabCol;

  int col = tabCol.size();
  int i = 0, j = 0;
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
    {
      i = (it->first) / col;
      j = (it->first) - i * col;
      if (isBlockAllocatedIn[i * col + j] == true)
      {
        delete(map(i, j));
      }
      map(i, j) = new MySimpleMatrix(*(it->second));
      isBlockAllocatedIn[i * col + j] = true;
    }
  }
  else
  {
    //m.isBlock () = true because it's a BlockMatrix
    mapped::iterator1 it;
    mapped::iterator2 it2;
    for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if (isBlockAllocatedIn[i * col + j] == true)
        {
          delete(map(i, j));
        }
        map(i, j) = new MySimpleMatrix(**it2);
        isBlockAllocatedIn[i * col + j] = true;
      }
    }
  }
  return *this;
}

const MyBlockMatrix& MyBlockMatrix::operator *= (double d)
{


  int col = tabCol.size();
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      if (isBlockAllocatedIn[(it->first)] == true)
      {
        *(it->second) *= d;
      }
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
    int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if (isBlockAllocatedIn[i * col + j] == true)
        {
          (**it2) *= d;
        }
      }
    }
  }
  return *this;
}

const MyBlockMatrix& MyBlockMatrix::operator *= (int d)
{


  int col = tabCol.size();
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      if (isBlockAllocatedIn[(it->first)] == true)
      {
        *(it->second) *= d;
      }
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
    int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if (isBlockAllocatedIn[i * col + j] == true)
        {
          (**it2) *= d;
        }
      }
    }
  }
  return *this;
}

const MyBlockMatrix& MyBlockMatrix::operator /= (double d)
{
  int col = tabCol.size();
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      if (isBlockAllocatedIn[(it->first)] == true)
      {
        *(it->second) /= d;
      }
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
    int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if (isBlockAllocatedIn[i * col + j] == true)
        {
          (**it2) /= d;
        }
      }
    }
  }
  return *this;
}

const MyBlockMatrix& MyBlockMatrix::operator /= (int d)
{
  int col = tabCol.size();
  if (STDMAP == 1)
  {
    mapped::array_type::iterator it;

    for (it = (map.data()).begin(); it != (map.data()).end(); ++it)
    {
      if (isBlockAllocatedIn[(it->first)] == true)
      {
        *(it->second) /= d;
      }
    }
  }
  else
  {
    mapped::iterator1 it;
    mapped::iterator2 it2;
    int i = 0, j = 0;
    for (it = map.begin1(); it != map.end1(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); it2++)
      {
        i = it2.index1();
        j = it2.index2();
        if (isBlockAllocatedIn[i * col + j] == true)
        {
          (**it2) /= d;
        }
      }
    }
  }
  return *this;
}

const MyBlockMatrix& MyBlockMatrix::operator += (const MySiconosMatrix &m)
{
  if (m.size1() != size1() || m.size2() != size2())
  {
    SiconosMatrixException::selfThrow("operator += (const MySiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  if (m.isBlock() == true)
  {
    mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;
    int col = tabCol.size();

    if (STDMAP == 1)
    {
      mapped::array_type::iterator it;
      int i = 0, j = 0;


      for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
      {
        i = (it->first) / col;
        j = (it->first) - i * col;
        if (isBlockAllocatedIn[(it->first)] == true)
        {
          *map(i, j) +=  *(it->second);
        }
      }
    }
    else
    {
      mapped::iterator1 it;
      mapped::iterator2 it2;
      int i = 0, j = 0;
      for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
      {
        for (it2 = it.begin(); it2 != it.end(); it2++)
        {
          i = it2.index1();
          j = it2.index2();
          if (isBlockAllocatedIn[i * col + j] == true)
          {
            *map(i, j) += (**it2);
          }
        }
      }
    }
  }

  else
  {
    //SUM USING OPERATOR ()
    for (int i = 0; i < size1(); i++)
    {

      for (int j = 0; j < size2(); j++)
      {
        (*this)(i, j) += m(i, j);
      }
    }
  }
  return *this;
}

const MyBlockMatrix& MyBlockMatrix::operator -= (const MySiconosMatrix &m)
{
  if (m.size1() != size1() || m.size2() != size2())
  {
    SiconosMatrixException::selfThrow("operator -= (const MySiconosMatrix&): Left and Right values have inconsistent sizes.");
  }

  if (m.isBlock() == true)
  {
    mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).map;
    int col = tabCol.size();

    if (STDMAP == 1)
    {
      mapped::array_type::iterator it;
      int i = 0, j = 0;

      for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
      {
        i = (it->first) / col;
        j = (it->first) - i * col;
        if (isBlockAllocatedIn[(it->first)] == true)
        {
          *map(i, j) -=  *(it->second);
        }
      }
    }
    else
    {
      mapped::iterator1 it;
      mapped::iterator2 it2;
      int i = 0, j = 0;
      for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
      {
        for (it2 = it.begin(); it2 != it.end(); it2++)
        {
          i = it2.index1();
          j = it2.index2();
          if (isBlockAllocatedIn[i * col + j] == true)
          {
            *map(i, j) -= (**it2);
          }
        }
      }
    }
  }

  else
  {
    //SUBSTRACTION USING OPERATOR ()
    for (int i = 0; i < size1(); i++)
    {
      for (int j = 0; j < size2(); j++)
      {
        (*this)(i, j) -= m(i, j);
      }
    }
  }
  return *this;
}



