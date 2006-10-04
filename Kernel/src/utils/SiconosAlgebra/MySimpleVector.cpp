
#include "MySimpleVector.h"
#include "MySimpleMatrix.h"
#include "ioVector.h"

using namespace boost::numeric::ublas;

/***************************** CONSTRUCTORS ****************************/
// Default (private)
MySimpleVector::MySimpleVector(TYP typ)
{
  setIsBlock(false);
  if (typ == DENSE)
  {
    vect.Dense = new DenseVect();
    num = 1;
  }
  else if (typ == SPARSE)
  {
    vect.Sparse = new SparseVect();
    num = 2;
  }
  else
    SiconosVectorException::selfThrow("constructor(TYP) : invalid type given");

}

MySimpleVector::MySimpleVector(const MySimpleVector &svect)
{
  setIsBlock(false);
  if (svect.getNum() == 1)
  {
    vect.Dense = new DenseVect(svect.getDense());
    num = 1;
  }
  else if (svect.getNum() == 2)
  {
    vect.Sparse = new SparseVect(svect.getSparse());
    num = 2;
  }
  else
    SiconosVectorException::selfThrow("constructor(const MySimpleVector) : invalid type given");
}

MySimpleVector::MySimpleVector(const MySiconosVector &svect)
{
  setIsBlock(false);
  assert(svect.isBlock() == false);
  if (svect.getNum() == 1)
  {
    vect.Dense = new DenseVect(svect.getDense());
    num = 1;
  }
  else if (svect.getNum() == 2)
  {
    vect.Sparse = new SparseVect(svect.getSparse());
    num = 2;
  }
  else
    SiconosVectorException::selfThrow("constructor(const MySiconosVector) : invalid type given");
}

MySimpleVector::MySimpleVector(const DenseVect& m)
{
  setIsBlock(false);
  vect.Dense = new DenseVect(m);
  num = 1;
}

MySimpleVector::MySimpleVector(const SparseVect& m)
{
  setIsBlock(false);
  vect.Sparse = new SparseVect(m);
  num = 2;
}

MySimpleVector::MySimpleVector(unsigned int row, TYP typ)
{
  setIsBlock(false);
  if (typ == SPARSE)
  {
    vect.Sparse = new SparseVect(row);
    num = 2;
  }
  else if (typ == DENSE)
  {
    vect.Dense = new DenseVect(row);
    num = 1;
  }
  else
  {
    SiconosVectorException::selfThrow("constructor(TYP, unsigned int) : invalid type given");
  }
}

MySimpleVector::MySimpleVector(const std::vector<double> v, unsigned int row, TYP typ)
{
  setIsBlock(false);
  if (typ == DENSE)
  {
    vect.Dense = new DenseVect(row, v);
    num = 1;
  }
  else
  {
    SiconosVectorException::selfThrow("constructor(TYP, std::vector<double>, unsigned int) : invalid type given");
  }
}

MySimpleVector::MySimpleVector(const std::string &file, bool ascii)
{

  setIsBlock(false);
  num = 1;
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
  else if (num == 2)
    delete(vect.Sparse);
}

/******************************** METHODS ******************************/
unsigned int  MySimpleVector::getNum(void)const
{
  return num;
}
void MySimpleVector::setNum(unsigned int n)
{
  num = n;
}

void MySimpleVector::zero(void)
{
  unsigned int size = (*this).size();
  if (num == 1)
  {
    zero_vector<double> p(size);
    *vect.Dense = p;
  }
  else if (num == 2)
  {
    zero_vector<double> p(size);
    *vect.Sparse = p;
  }
}

unsigned int MySimpleVector::size(void)const
{
  unsigned int n;
  if (num == 1)
    n = (vect.Dense)->size();
  if (num == 2)
    n = (vect.Sparse)->size();

  return n;
}

void MySimpleVector::resize(unsigned int n, bool preserve)
{

  if (num == 1)
    (vect.Dense)->resize(n, preserve);
  if (num == 2)
    (vect.Sparse)->resize(n, preserve);
}

const DenseVect MySimpleVector::getDense(void)const
{

  if (num != 1)
    SiconosVectorException::selfThrow("DenseVect getDense(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return *vect.Dense;
}

const SparseVect MySimpleVector::getSparse(void)const
{

  if (num != 2)
    SiconosVectorException::selfThrow("SparseVect getSparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return *vect.Sparse;
}

DenseVect* MySimpleVector::getDensePtr(void) const
{
  if (num != 1)
    SiconosVectorException::selfThrow("DenseVect* getDensePtr(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return vect.Dense;
}

SparseVect* MySimpleVector::getSparsePtr(void)const
{

  if (num != 2)
    SiconosVectorException::selfThrow("SparseVect* getSparsePtr(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return vect.Sparse;
}

const double MySimpleVector::normInf(void)const
{
  double d;
  if (num == 1)
    d = norm_inf(*vect.Dense);
  else if (num == 2)
    d = norm_inf(*vect.Sparse);
  return d;
}


void MySimpleVector::display(void)const
{
  std::cout << "vect: " ;
  if (num == 1)
    std::cout << *vect.Dense << std::endl;
  else if (num == 2)
    std::cout << *vect.Sparse << std::endl;
}


/***************************** OPERATORS ******************************/

double& MySimpleVector::operator()(unsigned int row)
{
  double d;  // \Warning : warning at compile due to reference return?
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
    break;
  default:
    SiconosVectorException::selfThrow("operator() (unsigned int) : invalid type given");

    break;
  }
  return d;
}

const MySimpleVector& MySimpleVector::operator = (const MySiconosVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense = m.getDense();
      break;
    case 2:
      *vect.Dense = m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator = : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 1:
      *vect.Sparse = m.getDense();
      break;
    case 2:
      *vect.Sparse = m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator = : invalid type given");
      break;
    }
    break;
  default:
    SiconosVectorException::selfThrow("operator = : invalid type given");
    break;
  }
  return *this;
}

const MySimpleVector& MySimpleVector::operator = (const MySimpleVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense = m.getDense();
      break;
    case 2:
      *vect.Dense = m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator = : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 1:
      *vect.Sparse = m.getDense();
      break;
    case 2:
      *vect.Sparse = m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator = : invalid type given");
      break;
    }
    break;
  default:
    SiconosVectorException::selfThrow("operator = : invalid type given");
    break;
  }
  return *this;
}

bool operator == (const MySiconosVector &m, const MySiconosVector &x)
{
  assert(m.isBlock() == false && x.isBlock() == false);
  if (m.isBlock() ==  true || m.isBlock() == true)
    SiconosVectorException::selfThrow("operator == : SimpleVector == BlockVector is forbidden");

  double norm = (m - x).normInf();
  return (norm < tolerance);
}

const MySimpleVector& MySimpleVector::operator += (const MySiconosVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense += m.getDense();
      break;
    case 2:
      *vect.Dense += m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator += : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 1:
      *vect.Sparse += m.getDense();
      break;
    case 2:
      *vect.Sparse += m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator += : invalid type given");
      break;
    }
    break;
  default:
    SiconosVectorException::selfThrow("operator += : invalid type given");
    break;
  }
  return *this;
}

const MySimpleVector& MySimpleVector::operator -= (const MySiconosVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *vect.Dense -= m.getDense();
      break;
    case 2:
      *vect.Dense -= m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator -= : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 1:
      *vect.Sparse -= m.getDense();
      break;
    case 2:
      *vect.Sparse -= m.getSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator -= : invalid type given");
      break;
    }
    break;
  default:
    SiconosVectorException::selfThrow("operator -= : invalid type given");
    break;
  }
  return *this;
}

const MySimpleVector& MySimpleVector::operator *= (double m)
{
  switch (num)
  {
  case 1:
    *vect.Dense *= m;
    break;
  case 2:
    *vect.Sparse *= m;
    break;
  default:
    SiconosVectorException::selfThrow("operator *=  : invalid type given");
    break;
  }
  return *this;
}

const MySimpleVector& MySimpleVector::operator *= (int m)
{
  switch (num)
  {
  case 1:
    *vect.Dense *= m;
    break;
  case 2:
    *vect.Sparse *= m;
    break;
  default:
    SiconosVectorException::selfThrow("operator *=  : invalid type given");
    break;
  }
  return *this;
}

const MySimpleVector& MySimpleVector::operator /= (double m)
{
  switch (num)
  {
  case 1:
    *vect.Dense /= m;
    break;
  case 2:
    *vect.Sparse /= m;
    break;
  default:
    SiconosVectorException::selfThrow("operator /=  : invalid type given");
    break;
  }
  return *this;
}

const MySimpleVector& MySimpleVector::operator /= (int m)
{
  switch (num)
  {
  case 1:
    *vect.Dense /= m;
    break;
  case 2:
    *vect.Sparse /= m;
    break;
  default:
    SiconosVectorException::selfThrow("operator /=  : invalid type given");
    break;
  }
  return *this;
}

MySimpleVector operator + (const MySimpleVector &x, const MySimpleVector &m)
{
  DenseVect p;
  SparseVect s;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("Vector addition: inconsistent sizes");

  if (x.getNum() != m.getNum())
    SiconosVectorException::selfThrow("Vector addition: use function add in order to add vectors of different type");

  if (x.getNum() == 1)
  {
    p = x.getDense() + m.getDense();
    return p;
  }
  else if (x.getNum() == 2)
  {
    s = x.getSparse() + m.getSparse();
    return s;
  }
}

MySimpleVector operator - (const MySimpleVector &x, const MySimpleVector &m)
{
  DenseVect p;
  SparseVect s;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("Vector subtraction: inconsistent sizes");

  if (x.getNum() != m.getNum())
    SiconosVectorException::selfThrow("Vector subtraction: use function sub in order to subtract vectors of different type");

  if (x.getNum() == 1)
  {
    p = x.getDense() - m.getDense();
    return p;
  }
  else if (x.getNum() == 2)
  {
    s = x.getSparse() - m.getSparse();
    return s;
  }
}

double MySimpleVector::operator * (const MySiconosVector &x)
{
  double p;

  if (x.size() != size())
    SiconosVectorException::selfThrow("operator * (const MySiconosVector&): inconsistent sizes");

  if (x.getNum() == num)
  {
    if (x.getNum() == 1)
    {
      p = inner_prod(x.getDense(), *vect.Dense);
    }
    else if (x.getNum() == 2)
    {
      p = inner_prod(x.getSparse(), *vect.Sparse);
    }
  }
  else
  {
    SiconosVectorException::selfThrow("operator (const MySiconosVector&): use function inner_prod or outer_prod in order to multiply vectors of different type");
  }
  return p;
}

MySimpleVector add(const MySiconosVector &x, const MySiconosVector& m)
{
  DenseVect p;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("vector function add : inconsistent sizes");

  if (m.getNum() == 1)
  {
    DenseVect q;
    q = m.getDense();
    if (x.getNum() == 1)
    {
      p = x.getDense() + q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getSparse() + q;
    }
    else
      SiconosVectorException::selfThrow("vector function add : invalid type of vector");
  }
  else if (m.getNum() == 2)
  {
    SparseVect q;
    q = m.getSparse();
    if (x.getNum() == 1)
    {
      p = x.getDense() + q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getSparse() + q;
    }
    else
      SiconosVectorException::selfThrow("vector function add : invalid type of vector");
  }
  else
  {
    SiconosVectorException::selfThrow("vector function add : invalid type of vector");
  }
  return p;
}

MySimpleVector sub(const MySiconosVector &x, const MySiconosVector& m)
{
  DenseVect p;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("vector function sub : inconsistent sizes");

  if (m.getNum() == 1)
  {
    DenseVect q;
    q = m.getDense();
    if (x.getNum() == 1)
    {
      p = x.getDense() - q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getSparse() - q;
    }
    else
      SiconosVectorException::selfThrow("vector function sub : invalid type of vector");

  }
  else if (m.getNum() == 2)
  {
    SparseVect q;
    q = m.getSparse();
    if (x.getNum() == 1)
    {
      p = x.getDense() - q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getSparse() - q;
    }
    else
      SiconosVectorException::selfThrow("vector function sub : invalid type of vector");

  }
  else
  {
    SiconosVectorException::selfThrow("vector function sub : invalid type of vector");
  }
  return p;
}

double inner_prod(const MySiconosVector &x, const MySiconosVector& m)
{
  double p;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("vector function inner_prod : inconsistent sizes");

  if (m.getNum() == 1)
  {
    DenseVect q;
    q = m.getDense();
    if (x.getNum() == 1)
    {
      p = inner_prod(x.getDense(), q);
    }
    else if (x.getNum() == 2)
    {
      p = inner_prod(x.getSparse(), q);
    }
    else
      SiconosVectorException::selfThrow("vector function inner_prod : invalid type of vector");
  }
  else if (m.getNum() == 2)
  {
    SparseVect q;
    q = m.getSparse();
    if (x.getNum() == 1)
    {
      p = inner_prod(x.getDense(), q);
    }
    else if (x.getNum() == 2)
    {
      p = inner_prod(x.getSparse(), q);
    }
    else
      SiconosVectorException::selfThrow("vector function inner_prod : invalid type of vector");

  }
  else
  {
    SiconosVectorException::selfThrow("vector function inner_prod : invalid type of vector");
  }
  return p;

}


MySimpleMatrix outer_prod(const MySiconosVector &x, const MySiconosVector& m)
{
  DenseMat p;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("vector function outer_prod : inconsistent sizes");

  if (m.getNum() == 1)
  {
    DenseVect q;
    q = m.getDense();
    if (x.getNum() == 1)
    {
      p = outer_prod(x.getDense(), q);
    }
    else if (x.getNum() == 2)
    {
      p = outer_prod(x.getSparse(), q);
    }
    else
      SiconosVectorException::selfThrow("vector function outer_prod : invalid type of vector");
  }
  else if (m.getNum() == 2)
  {
    SparseVect q;
    q = m.getSparse();
    if (x.getNum() == 1)
    {
      p = outer_prod(x.getDense(), q);
    }
    else if (x.getNum() == 2)
    {
      p = outer_prod(x.getSparse(), q);
    }
    else
      SiconosVectorException::selfThrow("vector function outer_prod : invalid type of vector");

  }
  else
  {
    SiconosVectorException::selfThrow("vector function outer_prod : invalid type of vector");
  }
  return p;

}

MySimpleVector operator * (const MySimpleVector &m, double d)
{
  DenseVect p;
  SparseVect s;
  if (m.getNum() != 1 || m.getNum() != 2)
    SiconosVectorException::selfThrow("opertor * (const MySimpleVector&, double) : invalid type of vector");

  if (m.getNum() == 1)
  {
    p = m.getDense() * d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    s = m.getSparse() * d;
    return s;
  }
}

MySimpleVector operator * (const MySimpleVector &m, int d)
{
  DenseVect p;
  SparseVect s;
  if (m.getNum() != 1 || m.getNum() != 2)
    SiconosVectorException::selfThrow("opertor * (const MySimpleVector&, int) : invalid type of vector");

  if (m.getNum() == 1)
  {
    p = m.getDense() * d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    s = m.getSparse() * d;
    return s;
  }
}

MySimpleVector operator * (double d, const MySimpleVector &m)
{
  DenseVect p;
  SparseVect s;
  if (m.getNum() != 1 || m.getNum() != 2)
    SiconosVectorException::selfThrow("opertor * (double, const MySimpleVector&) : invalid type of vector");

  if (m.getNum() == 1)
  {
    p = d * m.getDense();
    return p;
  }
  else if (m.getNum() == 2)
  {
    s = d * m.getSparse();
    return s;
  }
}

MySimpleVector operator * (int d, const MySimpleVector &m)
{
  DenseVect p;
  SparseVect s;
  if (m.getNum() != 1 || m.getNum() != 2)
    SiconosVectorException::selfThrow("opertor * (int, const MySimpleVector&) : invalid type of vector");

  if (m.getNum() == 1)
  {
    p = d * m.getDense();
    return p;
  }
  else if (m.getNum() == 2)
  {
    s = d * m.getSparse();
    return s;
  }
}

MySimpleVector operator / (const MySimpleVector &m, double d)
{
  DenseVect p;
  SparseVect s;
  if (m.getNum() != 1 || m.getNum() != 2)
    SiconosVectorException::selfThrow("opertor / (const MySimpleVector&, double) : invalid type of vector");

  if (m.getNum() == 1)
  {
    p = m.getDense() / d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    s = m.getSparse() / d;
    return s;
  }
}

MySimpleVector operator / (const MySimpleVector &m, int d)
{
  DenseVect p;
  SparseVect s;
  if (m.getNum() != 1 || m.getNum() != 2)
    SiconosVectorException::selfThrow("opertor / (const MySimpleVector&, int) : invalid type of vector");

  if (m.getNum() == 1)
  {
    p = m.getDense() / d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    s = m.getSparse() / d;
    return s;
  }
}



