/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "ioVector.h"
#include "SiconosVectorException.h"
#include <boost/numeric/ublas/io.hpp>            // for >> 
#include <boost/numeric/ublas/vector_proxy.hpp>  // for project

/***************************** CONSTRUCTORS ****************************/
// Default (private)
SimpleVector::SimpleVector(TYP typ): SiconosVector(false), num(1)
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
    SiconosVectorException::selfThrow("SimpleVector:constructor(TYP) : invalid type given");
  zero();
}

SimpleVector::SimpleVector(unsigned int row, TYP typ): SiconosVector(false), num(1)
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
    SiconosVectorException::selfThrow("SimpleVector::constructor(TYP, unsigned int) : invalid type given");
  }
  zero();
}

SimpleVector::SimpleVector(unsigned int row, double val, TYP typ): SiconosVector(false), num(1)
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
    SiconosVectorException::selfThrow("SimpleVector::constructor(TYP, unsigned int) : invalid type given");
  }
  fill(val);
}

SimpleVector::SimpleVector(const std::vector<double> v, TYP typ): SiconosVector(false), num(1)
{
  if (typ == DENSE)
  {
    vect.Dense = new DenseVect(v.size(), v);
    num = 1;
  }
  else
  {
    SiconosVectorException::selfThrow("SimpleVector::constructor(TYP, std::vector<double>, unsigned int) : invalid type given");
  }
}

SimpleVector::SimpleVector(const SimpleVector &svect): SiconosVector(false), num(svect.getNum())
{
  if (num == 1)
    vect.Dense = new DenseVect(svect.getDense());

  else if (num == 4)
    vect.Sparse = new SparseVect(svect.getSparse());

  else
    SiconosVectorException::selfThrow("SimpleVector:constructor(const SimpleVector) : invalid type given");
}

SimpleVector::SimpleVector(const SiconosVector &v): SiconosVector(false), num(1)
{
  if (v.isBlock())
  {
    vect.Dense = new DenseVect(v.size());
    for (unsigned int i = 0; i < size(); ++i)
      (*vect.Dense)(i) = v(i);
  }
  else
  {
    num = v.getNum();
    if (num == 1)
      vect.Dense = new DenseVect(v.getDense());

    else if (num == 4)
      vect.Sparse = new SparseVect(v.getSparse());

    else
      SiconosVectorException::selfThrow("SimpleVector:constructor(const SiconosVector) : invalid type given");
  }
}

SimpleVector::SimpleVector(const DenseVect& m): SiconosVector(false), num(1)
{
  vect.Dense = new DenseVect(m);
}

SimpleVector::SimpleVector(const SparseVect& m): SiconosVector(false), num(4)
{
  vect.Sparse = new SparseVect(m);
}

SimpleVector::SimpleVector(const std::string &file, bool ascii): SiconosVector(false), num(1)
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
SimpleVector::~SimpleVector(void)
{
  if (num == 1)
    delete(vect.Dense);
  else if (num == 4)
    delete(vect.Sparse);
}

/******************************** METHODS ******************************/
unsigned int SimpleVector::getNum() const
{
  return num;
}

void SimpleVector::zero(void)
{
  unsigned int size = (*this).size();
  if (num == 1)
  {
    ublas::zero_vector<double> p(size);
    *vect.Dense = p;
  }
  else //if(num==4)
  {
    ublas::zero_vector<double> p(size);
    *vect.Sparse = p;
  }
}

double* SimpleVector::getArray(unsigned int) const
{
  if (num == 4)
    SiconosVectorException::selfThrow("SimpleVector::getArray() : not yet implemented for sparse vector.");

  return &(((*vect.Dense).data())[0]);
}

void SimpleVector::getBlock(unsigned int pos, SiconosVector& vOut) const
{
  unsigned int end = vOut.size();
  if ((pos + end) > size())
    SiconosVectorException::selfThrow("SimpleVector::getBlock(pos,vOut): vOut size+pos is out of range.");

  if (num == 1)
    *(vOut.getDensePtr()) = project((*vect.Dense), ublas::range(pos, pos + end));
  else
    *(vOut.getSparsePtr()) = project((*vect.Sparse), ublas::range(pos, pos + end));
}

const double SimpleVector::normInf(void)const
{
  double d = 0;
  if (num == 1)
    d = norm_inf(*vect.Dense);
  else if (num == 4)
    d = norm_inf(*vect.Sparse);
  return d;
}


const double SimpleVector::norm(void)const
{
  double d = 0;
  if (num == 1)
    d = norm_2(*vect.Dense);
  else if (num == 4)
    d = norm_2(*vect.Sparse);
  return d;
}

unsigned int SimpleVector::size(void)const
{
  unsigned int n = 0;
  if (num == 1)
    n = (vect.Dense)->size();
  if (num == 4)
    n = (vect.Sparse)->size();

  return n;
}

void SimpleVector::resize(unsigned int n, bool preserve)
{

  if (num == 1)
    (vect.Dense)->resize(n, preserve);
  if (num == 4)
    (vect.Sparse)->resize(n, preserve);
}

const DenseVect SimpleVector::getDense(unsigned int)const
{

  if (num != 1)
    SiconosVectorException::selfThrow("SimpleVector::getDense(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return *vect.Dense;
}

const SparseVect SimpleVector::getSparse(unsigned int)const
{

  if (num != 4)
    SiconosVectorException::selfThrow("SimpleVector::getSparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return *vect.Sparse;
}

DenseVect* SimpleVector::getDensePtr(unsigned int) const
{
  if (num != 1)
    SiconosVectorException::selfThrow("SimpleVector::getDensePtr(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return vect.Dense;
}

SparseVect* SimpleVector::getSparsePtr(unsigned int)const
{

  if (num != 4)
    SiconosVectorException::selfThrow("SimpleVector::getSparsePtr(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return vect.Sparse;
}

void SimpleVector::display()const
{
  std::cout << "vect: " ;
  if (num == 1)
    std::cout << *vect.Dense << std::endl;
  else if (num == 4)
    std::cout << *vect.Sparse << std::endl;
}

void SimpleVector::fill(double value)
{
  if (num == 4)
  {
    for (unsigned int i = 0; i < (vect.Sparse)->size(); ++i)
      (vect.Sparse)->push_back(i, value);
  }
  else
  {
    for (unsigned int i = 0; i < (vect.Dense)->size(); ++i)
      (vect.Dense)->insert_element(i, value);
  }
}

std::string SimpleVector::toString() const
{
  std::stringstream sstr;
  std::string s;
  if (num == 1)
    sstr << *vect.Dense;
  else
    sstr << *vect.Sparse;
  sstr >> s;
  s = s.substr(4, s.size() - 5); // Remove "[size](" at the beginning of the string
  std::string::size_type pos;
  while ((pos = s.find(",")) != std::string::npos) // Replace "," by " " in the string
    s[pos] = ' ';
  return s;
}

/***************************** OPERATORS ******************************/

double SimpleVector::getValue(unsigned int row)
{
  if (num == 1)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row);
}

void SimpleVector::setValue(unsigned int row, double value)
{
  if (num == 1)
    (*vect.Dense)(row) = value ;
  else
    (*vect.Sparse)(row) = value;
}

double& SimpleVector::operator()(unsigned int row)
{
  if (num == 1)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

double SimpleVector::operator()(unsigned int row)const
{
  double d = 0;
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

SimpleVector& SimpleVector::operator = (const SiconosVector& m)
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
      SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse = m.getSparse();
    else
      SiconosVectorException::selfThrow("SimpleVector::operator = : can not set sparse = dense.");
    break;
  default:
    SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator = (const SimpleVector& m)
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
      SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse = m.getSparse();
    else
      SiconosVectorException::selfThrow("SimpleVector::operator = : can not set sparse = dense.");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator = (const DenseVect& d)
{
  if (num != 1)
    SiconosVectorException::selfThrow("SimpleVector::operator = DenseVect : current vector is not dense.");
  if (d.size() != size())
    SiconosVectorException::selfThrow("SimpleVector::operator = DenseVect : inconsistent size.");

  *vect.Dense = d;

  return *this;
}

SimpleVector& SimpleVector::operator = (const SparseVect& sp)
{
  if (num != 4)
    SiconosVectorException::selfThrow("SimpleVector::operator = SparseVect : current vector is not sparse.");
  if (sp.size() != size())
    SiconosVectorException::selfThrow("SimpleVector::operator = SparseVect : inconsistent size.");

  *vect.Sparse = sp;

  return *this;
}

bool operator == (const SiconosVector &m, const SiconosVector &x)
{
  assert(m.isBlock() == false && x.isBlock() == false);
  return ((m - x).norm() < tolerance);
}

SimpleVector& SimpleVector::operator += (const SiconosVector& m)
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
      SiconosVectorException::selfThrow("SimpleVector::operator += : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse += m.getSparse();
    else SiconosVectorException::selfThrow("SimpleVector::operator += : can not add a dense to a sparse.");
    break;
  default:
    SiconosVectorException::selfThrow("SimpleVector::operator += : invalid type given");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator -= (const SiconosVector& m)
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
      SiconosVectorException::selfThrow("SimpleVector::operator -= : invalid type given");
      break;
    }
    break;
  case 4:
    if (m.getNum() == 4)
      *vect.Sparse -= m.getSparse();
    else SiconosVectorException::selfThrow("SimpleVector::operator -= : can not sub a dense to a sparse.");
    break;
  default:
    SiconosVectorException::selfThrow("SimpleVector::operator -= : invalid type given");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator *= (double m)
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
    SiconosVectorException::selfThrow("SimpleVector::operator *=  : invalid type given");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator *= (int m)
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
    SiconosVectorException::selfThrow("SimpleVector::operator *=  : invalid type given");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator /= (double m)
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
    SiconosVectorException::selfThrow("SimpleVector::operator /=  : invalid type given");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator /= (int m)
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
    SiconosVectorException::selfThrow("SimpleVector::operator /=  : invalid type given");
    break;
  }
  return *this;
}

const double inner_prod(const SiconosVector &x, const SiconosVector &m)
{
  double p = 0;

  if (x.size() != m.size())
    SiconosVectorException::selfThrow("SimpleVector::operator * (const SiconosVector&): inconsistent sizes");

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

SimpleVector operator + (const SimpleVector &x, const SimpleVector &m)
{
  if (x.size() != m.size())
    SiconosVectorException::selfThrow("SimpleVector::Vector addition: inconsistent sizes");

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

SimpleVector operator - (const SimpleVector &x, const SimpleVector &m)
{
  if (x.size() != m.size())
    SiconosVectorException::selfThrow("SimpleVector::Vector subtraction: inconsistent sizes");

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
SimpleMatrix outer_prod(const SiconosVector &x, const SiconosVector& m)
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

SimpleVector operator * (const SimpleVector &m, double d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (const SimpleVector&, double) : invalid type of vector");

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

SimpleVector multScal(const SiconosVector&m , double d)
{
  unsigned int num = m.getNum();
  if (num == 1)
  {
    DenseVect p = d * m.getDense();
    return p;
  }
  else // if(num==4)
  {
    SparseVect s = d * m.getSparse();
    return s;
  }
}


SimpleVector operator * (const SimpleVector &m, int d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (const SimpleVector&, int) : invalid type of vector");

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

SimpleVector operator * (double d, const SimpleVector &m)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (double, const SimpleVector&) : invalid type of vector");

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

SimpleVector operator * (int d, const SimpleVector &m)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor * (int, const SimpleVector&) : invalid type of vector");

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

SimpleVector operator / (const SimpleVector &m, double d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor / (const SimpleVector&, double) : invalid type of vector");

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

SimpleVector operator / (const SimpleVector &m, int d)
{
  if (m.getNum() != 1 && m.getNum() != 4)
    SiconosVectorException::selfThrow("opertor / (const SimpleVector&, int) : invalid type of vector");

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



