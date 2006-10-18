
#include "MyBlockVector.h"
#include "MySimpleMatrix.h"
#include "ioVector.h"

using namespace boost::numeric::ublas;

/***************************** CONSTRUCTORS ****************************/
// Default (private)
MyBlockVector::MyBlockVector(): MySiconosVector(true)
{}

MyBlockVector::MyBlockVector(const std::string & file, bool ascii): MySiconosVector(true)
{
  SiconosVectorException::selfThrow(" MyBlockVector::constructor from a file : read BlockVector is not implemented");
}

MyBlockVector::MyBlockVector(const MyBlockVector &v): MySiconosVector(true)
{
  unsigned int numberOfBlocks;

  tabIndex = v.getTabIndex();
  numberOfBlocks = tabIndex.size();
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
  {
    // Call copy-constructor of MySimpleVector
    vect.push_back(new MySimpleVector(*v[i])) ;
    isBlockAllocatedIn.push_back(true);
  }
}

MyBlockVector::MyBlockVector(const MySiconosVector &v): MySiconosVector(true)
{
  if (v.isBlock())
  {
    //       unsigned int numberOfBlocks;

    //       MySiconosVector* tmp1 = const_cast<MySiconosVector*>(&v);
    //       MyBlockVector * tmp = static_cast<MyBlockVector*>(tmp1);
    //       tabIndex = tmp->getTabIndex();
    //       numberOfBlocks = tabIndex.size();
    //       for (unsigned int i=0;i < numberOfBlocks; ++i)
    //  {
    //    // Call copy-constructor of MySimpleVector
    //    vect.push_back(new MySimpleVector(*(tmp->getVectorPtr(i)))) ;
    //    isBlockAllocatedIn.push_back(true);
    //  }
    unsigned int numberOfBlocks = v.getNumberOfBlocks();;
    tabIndex.resize(numberOfBlocks);
    for (unsigned int i = 0; i < numberOfBlocks; ++i)
    {
      // Call copy-constructor of MySimpleVector
      vect.push_back(new MySimpleVector(*v[i])) ;
      isBlockAllocatedIn.push_back(true);
      tabIndex[i] = v[i]->size();
      if (i > 0) tabIndex[i] += tabIndex[i - 1];
    }
  }
  else
  {
    // Call copy-constructor of MySimpleVector
    vect.push_back(new MySimpleVector(v)) ;
    isBlockAllocatedIn.push_back(true);
    tabIndex.push_back(v.size());
  }
}

MyBlockVector::MyBlockVector(MySiconosVector* v1, MySiconosVector* v2): MySiconosVector(true)
{
  // Add the two vectors in the container
  // NO COPY !!
  if (v1->isBlock() || v2->isBlock())
    SiconosVectorException::selfThrow("MyBlockVector::constructor(v1,v2) : v1 or v2 is a block vector and cannot be added into another block vector.");

  vect.push_back(v1);
  vect.push_back(v2);
  isBlockAllocatedIn.resize(2, false);
  tabIndex.push_back(v1->size());
  tabIndex.push_back(v1->size() + v2->size());
}

MyBlockVector::MyBlockVector(unsigned int numberOfBlocks, unsigned int dim): MySiconosVector(true)
{
  vect.resize(numberOfBlocks, NULL);
  BlockVectIterator it;
  unsigned int i = 1;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    *it = new MySimpleVector(dim);
    tabIndex.push_back(dim * i);
    i++;
    isBlockAllocatedIn.push_back(true);
  }
}

/****************************** DESTRUCTOR  ****************************/
MyBlockVector::~MyBlockVector()
{
  unsigned int numberOfBlocks = tabIndex.size();
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
    if (isBlockAllocatedIn[i]) delete vect[i];
  vect.clear();
  tabIndex.clear();
  isBlockAllocatedIn.clear();
}

/******************************** METHODS ******************************/
unsigned int MyBlockVector::getNum() const
{
  return 0;
}

const DenseVect MyBlockVector::getDense(unsigned int i) const
{
  if (vect[i]->getNum() != 1)
    SiconosVectorException::selfThrow("MyBlockVector::getDense(unsigned int num) : the vector[num] is not a Dense vector");
  return (vect[i])->getDense();
}

const SparseVect MyBlockVector::getSparse(unsigned int i)const
{
  if (vect[i]->getNum() != 4)
    SiconosVectorException::selfThrow("MyBlockVector::getSparse(unsigned int num) : the vector[num] is not a Sparse vector");
  return (vect[i])->getSparse();
}

DenseVect* MyBlockVector::getDensePtr(unsigned int i) const
{
  if (vect[i]->getNum() != 1)
    SiconosVectorException::selfThrow("MyBlockVector::getDensePtr(unsigned int num) : the vector[num] is not a Dense vector");
  return (vect[i])->getDensePtr();
}

SparseVect* MyBlockVector::getSparsePtr(unsigned int i)const
{

  if (vect[i]->getNum() != 4)
    SiconosVectorException::selfThrow("MyBlockVector::getSparsePtr(unsigned int num) : the vector[num] is not a Sparse vector");
  return (vect[i])->getSparsePtr();
}

void MyBlockVector::zero()
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->zero();
}

bool MyBlockVector::check() const
{
  bool res = true;

  for (unsigned int i = 0; i < tabIndex.size(); ++i)
  {
    if (vect[i] == NULL)
    {
      res = false;
      std::cout << "BlockVector warning: one of the block is a NULL pointer, at position:" << i << std::endl;
    }
  }
  return res;
}

void MyBlockVector::resize(unsigned int, bool)
{
  SiconosVectorException::selfThrow("MyBlockVector::resize, not allowed for block vectors.");
}

const double MyBlockVector::normInf() const
{
  double d = 0, tmp;
  ConstBlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    tmp = (*it)->normInf();
    if (tmp > d) d = tmp;
  }
  return d;
}


const double MyBlockVector::norm() const
{
  double d = 0;
  ConstBlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    if ((*it) != NULL)
      d += pow((*it)->norm(), 2);
    else
      SiconosVectorException::selfThrow("MyBlockVector::norm, one of the blocks is equal to NULL pointer.");
  }
  return sqrt(d);
}

void MyBlockVector::display() const
{
  ConstBlockVectIterator it;
  std::cout << "=======> Block Vector Display (" << tabIndex.size() << " block(s)): " << std::endl;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->display();
}

MySimpleVector MyBlockVector::getVector(unsigned int num) const
{
  if (vect[num] == NULL)
    SiconosVectorException::selfThrow("MyBlockVector::getVector(num), vector[num] == NULL pointer.");

  return *(vect[num]);
}

MySiconosVector * MyBlockVector::getVectorPtr(unsigned int num)
{
  return vect[num];
}

void MyBlockVector::fill(double value)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    if ((*it) != NULL)(*it)->fill(value);
}

/***************************** OPERATORS ******************************/

double& MyBlockVector::operator()(unsigned int pos)
{
  unsigned int blockNum = 0;

  while (pos >= tabIndex[blockNum] && blockNum < tabIndex.size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= tabIndex[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

double MyBlockVector::operator()(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= tabIndex[blockNum] && blockNum < tabIndex.size())
    blockNum ++;
  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= tabIndex[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

MySiconosVector* MyBlockVector::operator [](unsigned int pos)
{
  return  vect[pos];
}

const MySiconosVector* MyBlockVector::operator [](unsigned int pos) const
{
  return  vect[pos];
}

MyBlockVector& MyBlockVector::operator = (const MySiconosVector& m)
{
  if (&m == this) return *this;
  if (!check())
    SiconosVectorException::selfThrow("MyBlockVector:operator = m, one of the block is a NULL pointer.");

  if (m.size() != size())
    SiconosVectorException::selfThrow("MyBlockVector:operator = m, inconsistent size between this and m.");

  if (m.isBlock())
  {
    double numberOfBlocks = tabIndex.size();
    for (unsigned int i = 0; i < numberOfBlocks; ++i)
    {
      // Call = of MySimpleVector
      *(vect[i]) = *(m[i]) ;
    }
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < size(); ++i)
      (*this)(i) = m(i);
  }

  return *this;
}

MyBlockVector& MyBlockVector::operator += (const MySiconosVector& m)
{
  if (m.size() != size())
    SiconosVectorException::selfThrow("MyBlockVector:operator += m, inconsistent size between this and m.");
  if (m.isBlock())
  {
    unsigned int n = getNumberOfBlocks();
    for (unsigned int i = 0; i < n; ++i)
      *vect[i] += *(m[i]);
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < size(); i++)
      (*this)(i) += m(i);
  }
  return *this;
}

MyBlockVector& MyBlockVector::operator -= (const MySiconosVector& m)
{
  if (m.size() != size())
    SiconosVectorException::selfThrow("MyBlockVector:operator -= m, inconsistent size between this and m.");
  if (m.isBlock())
  {
    unsigned int n = getNumberOfBlocks();
    for (unsigned int i = 0; i < n; ++i)
      *vect[i] -= *(m[i]);
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < size(); i++)
      (*this)(i) -= m(i);
  }
  return *this;
}

MyBlockVector& MyBlockVector::operator *= (double m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) *= m;

  return *this;
}

MyBlockVector& MyBlockVector::operator *= (int m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) *= m;

  return *this;
}

MyBlockVector& MyBlockVector::operator /= (double m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) /= m;

  return *this;
}

MyBlockVector& MyBlockVector::operator /= (int m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) /= m;

  return *this;
}

void MyBlockVector::add(const  MySiconosVector& v)
{
  if (v.isBlock())
    SiconosVectorException::selfThrow("MyBlockVector::add(v) : v is a block vector and cannot be added into another block vector.");

  unsigned int tabVal = size() + v.size();
  vect.push_back(new MySimpleVector(v)); // Copy
  isBlockAllocatedIn.push_back(true);
  tabIndex.push_back(tabVal);
}

void MyBlockVector::addPtr(MySiconosVector* v)
{
  if (v->isBlock())
    SiconosVectorException::selfThrow("MyBlockVector::addPtr(v) : v is a block vector and cannot be added into another block vector.");

  unsigned int tabVal = size() + v->size();
  vect.push_back(v);
  isBlockAllocatedIn.push_back(false);
  tabIndex.push_back(tabVal);
}
