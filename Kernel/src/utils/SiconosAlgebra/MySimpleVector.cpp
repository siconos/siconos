
#include "MySimpleVector.h"
#include "MySimpleMatrix.h"
#include "ioVector.h"

using namespace boost::numeric::ublas;

//Default private
MySimpleVector::MySimpleVector(void): num(0)
{
  SetIsBlock(false);
}

/***************************** CONSTRUCTORS ****************************/
MySimpleVector::MySimpleVector(TYP typ)
{
  SetIsBlock(false);
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
  SetIsBlock(false);
  if (svect.GetNum() == 1)
  {
    vect.Dense = new DenseVect(svect.GetDense());
    num = 1;
  }
  else if (svect.GetNum() == 2)
  {
    vect.Sparse = new SparseVect(svect.GetSparse());
    num = 2;
  }
  else
    SiconosVectorException::selfThrow("constructor(const MySimpleVector) : invalid type given");
}

MySimpleVector::MySimpleVector(const MySiconosVector &svect)
{
  SetIsBlock(false);
  assert(svect.isBlock() == false);
  if (svect.GetNum() == 1)
  {
    vect.Dense = new DenseVect(svect.GetDense());
    num = 1;
  }
  else if (svect.GetNum() == 2)
  {
    vect.Sparse = new SparseVect(svect.GetSparse());
    num = 2;
  }
  else
    SiconosVectorException::selfThrow("constructor(const MySiconosVector) : invalid type given");
}

MySimpleVector::MySimpleVector(const DenseVect& m)
{
  SetIsBlock(false);
  vect.Dense = new DenseVect(m);
  num = 1;
}

MySimpleVector::MySimpleVector(const SparseVect& m)
{
  SetIsBlock(false);
  vect.Sparse = new SparseVect(m);
  num = 2;
}

MySimpleVector::MySimpleVector(TYP typ, int row)
{
  SetIsBlock(false);
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
    SiconosVectorException::selfThrow("constructor(TYP, int) : invalid type given");
  }
}

MySimpleVector::MySimpleVector(TYP typ, const std::vector<double> v, int row)
{
  SetIsBlock(false);
  if (typ == DENSE)
  {
    vect.Dense = new DenseVect(row, v);
    num = 1;
  }
  else
  {
    SiconosVectorException::selfThrow("constructor(TYP, std::vector<double>, int) : invalid type given");
  }
}

MySimpleVector::MySimpleVector(const std::string &file, bool ascii)
{

  SetIsBlock(false);
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
int  MySimpleVector::GetNum(void)const
{
  return num;
}
void MySimpleVector::SetNum(int n)
{
  num = n;
}

void MySimpleVector::zero(void)
{
  int size = (*this).size();
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

int MySimpleVector::size(void)const
{
  int n;
  if (num == 1)
    n = (vect.Dense)->size();
  if (num == 2)
    n = (vect.Sparse)->size();

  return n;
}

void MySimpleVector::resize(int n, bool preserve)
{

  if (num == 1)
    (vect.Dense)->resize(n, preserve);
  if (num == 2)
    (vect.Sparse)->resize(n, preserve);
}

const DenseVect MySimpleVector::GetDense(void)const
{

  if (num != 1)
    SiconosVectorException::selfThrow("DenseVect GetDense(int row, int col) : the current vector is not a Dense vector");

  return *vect.Dense;
}

const SparseVect MySimpleVector::GetSparse(void)const
{

  if (num != 2)
    SiconosVectorException::selfThrow("SparseVect GetSparse(int row, int col) : the current vector is not a Sparse vector");

  return *vect.Sparse;
}

const DenseVect* MySimpleVector::GetDensePtr(void)const
{

  if (num != 1)
    SiconosVectorException::selfThrow("DenseVect* GetDensePtr(int row, int col) : the current vector is not a Dense vector");

  return vect.Dense;
}

const SparseVect* MySimpleVector::GetSparsePtr(void)const
{

  if (num != 2)
    SiconosVectorException::selfThrow("SparseVect* GetSparsePtr(int row, int col) : the current vector is not a Sparse vector");

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

double& MySimpleVector::operator()(int row)
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
    SiconosVectorException::selfThrow("operator() (int) : invalid type given");

    break;
  }
  return d;
}

double MySimpleVector::operator()(int row)const
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
    SiconosVectorException::selfThrow("operator() (int) : invalid type given");

    break;
  }
  return d;
}

const MySimpleVector& MySimpleVector::operator = (const MySiconosVector& m)
{
  switch (num)
  {
  case 1:
    switch (m.GetNum())
    {
    case 1:
      *vect.Dense = m.GetDense();
      break;
    case 2:
      *vect.Dense = m.GetSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator = : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *vect.Sparse = m.GetDense();
      break;
    case 2:
      *vect.Sparse = m.GetSparse();
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
    switch (m.GetNum())
    {
    case 1:
      *vect.Dense = m.GetDense();
      break;
    case 2:
      *vect.Dense = m.GetSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator = : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *vect.Sparse = m.GetDense();
      break;
    case 2:
      *vect.Sparse = m.GetSparse();
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
    switch (m.GetNum())
    {
    case 1:
      *vect.Dense += m.GetDense();
      break;
    case 2:
      *vect.Dense += m.GetSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator += : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *vect.Sparse += m.GetDense();
      break;
    case 2:
      *vect.Sparse += m.GetSparse();
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
    switch (m.GetNum())
    {
    case 1:
      *vect.Dense -= m.GetDense();
      break;
    case 2:
      *vect.Dense -= m.GetSparse();
      break;
    default:
      SiconosVectorException::selfThrow("operator -= : invalid type given");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *vect.Sparse -= m.GetDense();
      break;
    case 2:
      *vect.Sparse -= m.GetSparse();
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

  if (x.GetNum() == m.GetNum())
  {
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + m.GetDense();
      return p;
    }
    else if (x.GetNum() == 2)
    {
      s = x.GetSparse() + m.GetSparse();
      return s;
    }
  }
  else
    SiconosVectorException::selfThrow("Vector addition: use function add in order to add vectors of different type");

}

MySimpleVector operator - (const MySimpleVector &x, const MySimpleVector &m)
{
  DenseVect p;
  SparseVect s;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("Vector subtraction: inconsistent sizes");

  if (x.GetNum() == m.GetNum())
  {
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - m.GetDense();
      return p;
    }
    else if (x.GetNum() == 2)
    {
      s = x.GetSparse() - m.GetSparse();
      return s;
    }
  }
  else
    SiconosVectorException::selfThrow("Vector subtraction: use function sub in order to subtract vectors of different type");
}

double MySimpleVector::operator * (const MySiconosVector &x)
{
  double p;

  if (x.size() != size())
    SiconosVectorException::selfThrow("operator * (const MySiconosVector&): inconsistent sizes");

  if (x.GetNum() == num)
  {
    if (x.GetNum() == 1)
    {
      p = inner_prod(x.GetDense(), *vect.Dense);
    }
    else if (x.GetNum() == 2)
    {
      p = inner_prod(x.GetSparse(), *vect.Sparse);
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

  if (m.GetNum() == 1)
  {
    DenseVect q;
    q = m.GetDense();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetSparse() + q;
    }
    else
      SiconosVectorException::selfThrow("vector function add : invalid type of vector");
  }
  else if (m.GetNum() == 2)
  {
    SparseVect q;
    q = m.GetSparse();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetSparse() + q;
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

  if (m.GetNum() == 1)
  {
    DenseVect q;
    q = m.GetDense();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetSparse() - q;
    }
    else
      SiconosVectorException::selfThrow("vector function sub : invalid type of vector");

  }
  else if (m.GetNum() == 2)
  {
    SparseVect q;
    q = m.GetSparse();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetSparse() - q;
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

  if (m.GetNum() == 1)
  {
    DenseVect q;
    q = m.GetDense();
    if (x.GetNum() == 1)
    {
      p = inner_prod(x.GetDense(), q);
    }
    else if (x.GetNum() == 2)
    {
      p = inner_prod(x.GetSparse(), q);
    }
    else
      SiconosVectorException::selfThrow("vector function inner_prod : invalid type of vector");
  }
  else if (m.GetNum() == 2)
  {
    SparseVect q;
    q = m.GetSparse();
    if (x.GetNum() == 1)
    {
      p = inner_prod(x.GetDense(), q);
    }
    else if (x.GetNum() == 2)
    {
      p = inner_prod(x.GetSparse(), q);
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

  if (m.GetNum() == 1)
  {
    DenseVect q;
    q = m.GetDense();
    if (x.GetNum() == 1)
    {
      p = outer_prod(x.GetDense(), q);
    }
    else if (x.GetNum() == 2)
    {
      p = outer_prod(x.GetSparse(), q);
    }
    else
      SiconosVectorException::selfThrow("vector function outer_prod : invalid type of vector");
  }
  else if (m.GetNum() == 2)
  {
    SparseVect q;
    q = m.GetSparse();
    if (x.GetNum() == 1)
    {
      p = outer_prod(x.GetDense(), q);
    }
    else if (x.GetNum() == 2)
    {
      p = outer_prod(x.GetSparse(), q);
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
  if (m.GetNum() == 1)
  {
    p = m.GetDense() * d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    s = m.GetSparse() * d;
    return s;
  }
  else
    SiconosVectorException::selfThrow("opertor * (const MySimpleVector&, double) : invalid type of vector");
}

MySimpleVector operator * (const MySimpleVector &m, int d)
{
  DenseVect p;
  SparseVect s;
  if (m.GetNum() == 1)
  {
    p = m.GetDense() * d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    s = m.GetSparse() * d;
    return s;
  }
  else
    SiconosVectorException::selfThrow("opertor * (const MySimpleVector&, int) : invalid type of vector");
}

MySimpleVector operator * (double d, const MySimpleVector &m)
{
  DenseVect p;
  SparseVect s;
  if (m.GetNum() == 1)
  {
    p = d * m.GetDense();
    return p;
  }
  else if (m.GetNum() == 2)
  {
    s = d * m.GetSparse();
    return s;
  }
  else
    SiconosVectorException::selfThrow("opertor * (double, const MySimpleVector&) : invalid type of vector");
}

MySimpleVector operator * (int d, const MySimpleVector &m)
{
  DenseVect p;
  SparseVect s;
  if (m.GetNum() == 1)
  {
    p = d * m.GetDense();
    return p;
  }
  else if (m.GetNum() == 2)
  {
    s = d * m.GetSparse();
    return s;
  }
  else
    SiconosVectorException::selfThrow("opertor * (int, const MySimpleVector&) : invalid type of vector");
}

MySimpleVector operator / (const MySimpleVector &m, double d)
{
  DenseVect p;
  SparseVect s;
  if (m.GetNum() == 1)
  {
    p = m.GetDense() / d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    s = m.GetSparse() / d;
    return s;
  }
  else
    SiconosVectorException::selfThrow("opertor / (const MySimpleVector&, double) : invalid type of vector");
}

MySimpleVector operator / (const MySimpleVector &m, int d)
{
  DenseVect p;
  SparseVect s;
  if (m.GetNum() == 1)
  {
    p = m.GetDense() / d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    s = m.GetSparse() / d;
    return s;
  }
  else
    SiconosVectorException::selfThrow("opertor / (const MySimpleVector&, int) : invalid type of vector");
}



