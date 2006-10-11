
#include "MySimpleVector.h"
#include "MySimpleMatrix.h"
#include "ioVector.h"

using namespace boost::numeric::ublas;

/***************************** CONSTRUCTORS ****************************/
// Default (private)
MySimpleVector::MySimpleVector(TYP typ): MySiconosVector(false), num(1)
{
  if (typ == DENSE)
  {
    vect.Dense = new DenseVect();
    num = 1;
  }
  else if (typ == SPARSE)
  {
    vect.Sparse = new SparseVect();
    num = 4;
  }
  else
    SiconosVectorException::selfThrow("MySimpleVector:constructor(TYP) : invalid type given");
  zero();
}

MySimpleVector::MySimpleVector(unsigned int row, TYP typ): MySiconosVector(false), num(1)
{
  if (typ == SPARSE)
  {
    vect.Sparse = new SparseVect(row);
    num = 4;
  }
  else if (typ == DENSE)
  {
    vect.Dense = new DenseVect(row);
    num = 1;
  }
  else
  {
    SiconosVectorException::selfThrow("MySimpleVector::constructor(TYP, unsigned int) : invalid type given");
  }
  zero();
}

MySimpleVector::MySimpleVector(const std::vector<double> v, TYP typ): MySiconosVector(false), num(1)
{
  if (typ == DENSE)
  {
    vect.Dense = new DenseVect(v.size(), v);
    num = 1;
  }
  else
  {
    SiconosVectorException::selfThrow("MySimpleVector::constructor(TYP, std::vector<double>, unsigned int) : invalid type given");
  }
}

MySimpleVector::MySimpleVector(const MySimpleVector &svect): MySiconosVector(false), num(svect.getNum())
{
  if (num == 1)
    vect.Dense = new DenseVect(svect.getDense());

  else if (num == 4)
    vect.Sparse = new SparseVect(svect.getSparse());

  else
    SiconosVectorException::selfThrow("MySimpleVector:constructor(const MySimpleVector) : invalid type given");
}

MySimpleVector::MySimpleVector(const MySiconosVector &svect): MySiconosVector(false), num(svect.getNum())
{
  assert(svect.isBlock() == false);
  if (num == 1)
    vect.Dense = new DenseVect(svect.getDense());

  else if (num == 4)
    vect.Sparse = new SparseVect(svect.getSparse());

  else
    SiconosVectorException::selfThrow("MySimpleVector:constructor(const MySiconosVector) : invalid type given");
}

MySimpleVector::MySimpleVector(const DenseVect& m): MySiconosVector(false), num(1)
{
  vect.Dense = new DenseVect(m);
}

MySimpleVector::MySimpleVector(const SparseVect& m): MySiconosVector(false), num(4)
{
  vect.Sparse = new SparseVect(m);
}

MySimpleVector::MySimpleVector(const std::string &file, bool ascii): MySiconosVector(false), num(1)
{
  vect.Dense = new DenseVect();
  if (ascii)
  {
    ioVector io(file, "ascii");
    io.read(*this);
  }
  else
  {
    ioVector io(file, "binary");
    io.read(*this);
  }
}

/****************************** DESTRUCTOR  ****************************/
MySimpleVector::~MySimpleVector(void)
{
  if (num == 1)
    delete(vect.Dense);
  else if (num == 4)
    delete(vect.Sparse);
}

/******************************** METHODS ******************************/
unsigned int  MySimpleVector::getNum(void)const
{
  return num;
}

void MySimpleVector::zero(void)
{
  unsigned int size = (*this).size();
  if (num == 1)
  {
    zero_vector<double> p(size);
    *vect.Dense = p;
  }
  else //if(num==4)
  {
    zero_vector<double> p(size);
    *vect.Sparse = p;
  }
}


const double MySimpleVector::normInf(void)const
{
  double d;
  if (num == 1)
    d = norm_inf(*vect.Dense);
  else if (num == 4)
    d = norm_inf(*vect.Sparse);
  return d;
}


const double MySimpleVector::norm(void)const
{
  double d;
  if (num == 1)
    d = norm_2(*vect.Dense);
  else if (num == 4)
    d = norm_2(*vect.Sparse);
  return d;
}

unsigned int MySimpleVector::size(void)const
{
  unsigned int n;
  if (num == 1)
    n = (vect.Dense)->size();
  if (num == 4)
    n = (vect.Sparse)->size();

  return n;
}

void MySimpleVector::resize(unsigned int n, bool preserve)
{

  if (num == 1)
    (vect.Dense)->resize(n, preserve);
  if (num == 4)
    (vect.Sparse)->resize(n, preserve);
}

const DenseVect MySimpleVector::getDense(void)const
{

  if (num != 1)
    SiconosVectorException::selfThrow("MySimpleVector::getDense(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return *vect.Dense;
}

const SparseVect MySimpleVector::getSparse(void)const
{

  if (num != 4)
    SiconosVectorException::selfThrow("MySimpleVector::getSparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return *vect.Sparse;
}

DenseVect* MySimpleVector::getDensePtr(void) const
{
  if (num != 1)
    SiconosVectorException::selfThrow("MySimpleVector::getDensePtr(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return vect.Dense;
}

SparseVect* MySimpleVector::getSparsePtr(void)const
{

  if (num != 4)
    SiconosVectorException::selfThrow("MySimpleVector::getSparsePtr(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return vect.Sparse;
}

void MySimpleVector::display(void)const
{
  std::cout << "vect: " ;
  if (num == 1)
    std::cout << *vect.Dense << std::endl;
  else if (num == 4)
    std::cout << *vect.Sparse << std::endl;
}


/***************************** OPERATORS ******************************/

double& MySimpleVector::operator()(unsigned int row)
{
  if (num == 1)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

double MySimpleVector::operator()(unsigned int row)const
{
  double d;
  switch (num)
  {
  case 1:
    d = (*vect.Dense)(row);
    break;
  case 2:
    d = ((*vect.Sparse)(row)).ref();
    break;
  default:
    SiconosVectorException::selfThrow("operator() (unsigned int) : invalid type given");
    break;
  }
  return d;
}

MySimpleVector& MySimpleVector::operator = (const MySiconosVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense = m.getDense();
      break;
    case 4:
      *vect.Dense = m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("MySimpleVector::operator = : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse = m.getSparse();
    else
      SiconosVectorException::selfThrow("MySimpleVector::operator = : can not set sparse = dense.");
    break;
  default:
    SiconosVectorException::selfThrow("MySimpleVector::operator = : invalid type given");
    break;
  }
  return *this;
}

MySimpleVector& MySimpleVector::operator = (const MySimpleVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense = m.getDense();
      break;
    case 4:
      *vect.Dense = m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("MySimpleVector::operator = : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse = m.getSparse();
    else
      SiconosVectorException::selfThrow("MySimpleVector::operator = : can not set sparse = dense.");
    break;
  }
  return *this;
}

bool operator == (const MySiconosVector &m, const MySiconosVector &x)
{
  assert(m.isBlock() == false && x.isBlock() == false);
  return ((m - x).norm() < tolerance);
}

MySimpleVector& MySimpleVector::operator += (const MySiconosVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense += m.getDense();
      break;
    case 4:
      *vect.Dense += m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("MySimpleVector::operator += : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse += m.getSparse();
    else SiconosVectorException::selfThrow("MySimpleVector::operator += : can not add a dense to a sparse.");
    break;
  default:
    SiconosVectorException::selfThrow("MySimpleVector::operator += : invalid type given");
    break;
  }
  return *this;
}

MySimpleVector& MySimpleVector::operator -= (const MySiconosVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense -= m.getDense();
      break;
    case 4:
      *vect.Dense -= m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("MySimpleVector::operator -= : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse -= m.getSparse();
    else SiconosVectorException::selfThrow("MySimpleVector::operator -= : can not sub a dense to a sparse.");
    break;
  default:
    SiconosVectorException::selfThrow("MySimpleVector::operator -= : invalid type given");
    break;
  }
  return *this;
}

MySimpleVector& MySimpleVector::operator *= (double m)
{
  switch (num)
  {
  case 1:
    *vect.Dense *= m;
    break;
  case 4:
    *vect.Sparse *= m;
    break;
  default:
    SiconosVectorException::selfThrow("MySimpleVector::operator *=  : invalid type given");
    break;
  }
  return *this;
}

MySimpleVector& MySimpleVector::operator *= (int m)
{
  switch (num)
  {
  case 1:
    *vect.Dense *= m;
    break;
  case 4:
    *vect.Sparse *= m;
    break;
  default:
    SiconosVectorException::selfThrow("MySimpleVector::operator *=  : invalid type given");
    break;
  }
  return *this;
}

MySimpleVector& MySimpleVector::operator /= (double m)
{
  switch (num)
  {
  case 1:
    *vect.Dense /= m;
    break;
  case 4:
    *vect.Sparse /= m;
    break;
  default:
    SiconosVectorException::selfThrow("MySimpleVector::operator /=  : invalid type given");
    break;
  }
  return *this;
}

MySimpleVector& MySimpleVector::operator /= (int m)
{
  switch (num)
  {
  case 1:
    *vect.Dense /= m;
    break;
  case 4:
    *vect.Sparse /= m;
    break;
  default:
    SiconosVectorException::selfThrow("MySimpleVector::operator /=  : invalid type given");
    break;
  }
  return *this;
}

const double inner_prod(const MySiconosVector &x, const MySiconosVector &m)
{
  double p;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("MySimpleVector::operator * (const MySiconosVector&): inconsistent sizes");

  unsigned int num = m.getNum();
  if (x.getNum() == num)
  {
    if (num == 1)
    {
      p = inner_prod(x.getDense(), m.getDense());
    }
    else if (num == 4)
    {
      p = inner_prod(x.getSparse(), m.getSparse());
    }
  }
  else
  {
    if (num == 1)
    {
      p = inner_prod(x.getSparse(), m.getDense());
    }
    else if (num == 4)
    {
      p = inner_prod(x.getDense(), m.getSparse());
    }
  }
  return p;
}

MySimpleVector operator + (const MySimpleVector &x, const MySimpleVector &m)
{
  if (x.size() != m.size())
    SiconosVectorException::selfThrow("MySimpleVector::Vector addition: inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numM = m.getNum();

  if (numX == numM)
  {
    if (numX == 1)
    {
      DenseVect p = x.getDense() + m.getDense();
      return p;
    }
    else
    {
      SparseVect s = x.getSparse() + m.getSparse();
      return s;
    }
  }
  else
  {
    DenseVect p;
    if (numX == 1)
    {
      p = x.getDense() + m.getSparse();
      return p;
    }
    else
    {
      p = x.getSparse() + m.getDense();
      return p;
    }
  }
}

MySimpleVector operator - (const MySimpleVector &x, const MySimpleVector &m)
{
  if (x.size() != m.size())
    SiconosVectorException::selfThrow("MySimpleVector::Vector subtraction: inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numM = m.getNum();
  if (numX == numM)
  {
    if (numX == 1)
    {
      DenseVect p = x.getDense() - m.getDense();
      return p;
    }
    else
    {
      SparseVect s = x.getSparse() - m.getSparse();
      return s;
    }
  }
  else
  {
    DenseVect p;
    if (numX == 1)
    {
      p = x.getDense() - m.getSparse();
      return p;
    }
    else
    {
      p = x.getSparse() - m.getDense();
      return p;
    }
  }
}

// outer_prod(v,w) = trans(v)*w
MySimpleMatrix outer_prod(const MySiconosVector &x, const MySiconosVector& m)
{
  DenseMat p;
  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();
  if (numM == 1)
  {
    if (numX == 1)
      p = outer_prod(x.getDense(), m.getDense());
    else if (numX == 4)
      p = outer_prod(x.getSparse(), m.getDense());
    else
      SiconosVectorException::selfThrow("vector function outer_prod : invalid type of vector");
  }
  else if (numM == 4)
  {
    if (numX == 1)
      p = outer_prod(x.getDense(), m.getSparse());
    else if (numX == 4)
      p = outer_prod(x.getSparse(), m.getSparse());
    else
      SiconosVectorException::selfThrow("vector function outer_prod : invalid type of vector");
  }
  else
    SiconosVectorException::selfThrow("vector function outer_prod : invalid type of vector");

  return p;
}

MySimpleVector operator * (const MySimpleVector &m, double d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (const MySimpleVector&, double) : invalid type of vector");

  if (m.getNum() == 1)
  {
    DenseVect p = m.getDense() * d;
    return p;
  }
  else// if(m.getNum()==4)
  {
    SparseVect s = m.getSparse() * d;
    return s;
  }
}

MySimpleVector operator * (const MySimpleVector &m, int d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (const MySimpleVector&, int) : invalid type of vector");

  if (m.getNum() == 1)
  {
    DenseVect p = m.getDense() * d;
    return p;
  }
  else // if(m.getNum()==4){
  {
    SparseVect s = m.getSparse() * d;
    return s;
  }
}

MySimpleVector operator * (double d, const MySimpleVector &m)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (double, const MySimpleVector&) : invalid type of vector");

  if (m.getNum() == 1)
  {
    DenseVect p = d * m.getDense();
    return p;
  }
  else // if(m.getNum()==4){
  {
    SparseVect s = d * m.getSparse();
    return s;
  }
}

MySimpleVector operator * (int d, const MySimpleVector &m)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (int, const MySimpleVector&) : invalid type of vector");

  if (m.getNum() == 1)
  {
    DenseVect p = d * m.getDense();
    return p;
  }
  else
  {
    // if(m.getNum()==4){
    SparseVect s = d * m.getSparse();
    return s;
  }
}

MySimpleVector operator / (const MySimpleVector &m, double d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor / (const MySimpleVector&, double) : invalid type of vector");

  if (m.getNum() == 1)
  {
    DenseVect p = m.getDense() / d;
    return p;
  }
  else // if(m.getNum()==4){
  {
    SparseVect s = m.getSparse() / d;
    return s;
  }
}

MySimpleVector operator / (const MySimpleVector &m, int d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor / (const MySimpleVector&, int) : invalid type of vector");

  if (m.getNum() == 1)
  {
    DenseVect p = m.getDense() / d;
    return p;
  }
  else // if(m.getNum()==4){
  {
    SparseVect s = m.getSparse() / d;
    return s;
  }
}



