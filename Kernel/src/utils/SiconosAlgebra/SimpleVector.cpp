/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include <boost/numeric/ublas/io.hpp>            // for >> 
#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <boost/numeric/bindings/atlas/cblas1.hpp>


#include "SimpleVector.hpp"
#include "SimpleMatrix.hpp"
#include "ioVector.hpp"


// =================================================
//                CONSTRUCTORS
// =================================================

// Default
SimpleVector::SimpleVector(): SiconosVector()
{
  _dense = true;
  vect.Dense = new DenseVect(ublas::zero_vector<double>());
}

// parameters: dimension and type.
SimpleVector::SimpleVector(unsigned int row, Siconos::UBLAS_TYPE typ): SiconosVector()
{
  if (typ == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(ublas::zero_vector<double>(row));
  }
  else if (typ == Siconos::DENSE)
  {
    _dense = true;
    vect.Dense = new DenseVect(ublas::zero_vector<double>(row));
  }
  else
  {
    SiconosVectorException::selfThrow("SimpleVector::constructor(Siconos::UBLAS_TYPE, unsigned int) failed, invalid type given");
  }
}

// parameters: dimension, default value for all components and type.
SimpleVector::SimpleVector(unsigned int row, double val, Siconos::UBLAS_TYPE typ): SiconosVector()
{
  if (typ == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(row);
    fill(val);
  }
  else if (typ == Siconos::DENSE)
  {
    _dense = true;
    vect.Dense = new DenseVect(ublas::scalar_vector<double>(row, val));
  }
  else
  {
    SiconosVectorException::selfThrow("SimpleVector::constructor(Siconos::UBLAS_TYPE, unsigned int) : invalid type given");
  }
}

// parameters: a vector (stl) of double and the type.
SimpleVector::SimpleVector(const std::vector<double>& v, Siconos::UBLAS_TYPE typ): SiconosVector()
{
  if (typ != Siconos::DENSE)
    SiconosVectorException::selfThrow("SimpleVector::constructor(Siconos::UBLAS_TYPE, std::vector<double>, unsigned int) : invalid type given");

  _dense = true;
  vect.Dense = new DenseVect(v.size());
  std::copy(v.begin(), v.end(), (vect.Dense)->begin());
}

// Copy
SimpleVector::SimpleVector(const SimpleVector &svect): SiconosVector()
{
  if (ask<IsDense>(svect)) // dense
  {
    _dense = true;
    vect.Dense = new DenseVect(svect.size());
    noalias(*vect.Dense) = (*svect.dense());
    // std::copy((vect.Dense)->begin(), (vect.Dense)->end(), (svect.dense())->begin());
  }
  else //sparse
  {
    _dense = false;
    vect.Sparse = new SparseVect(svect.size());
    noalias(*vect.Sparse) = (*svect.sparse());
    //std::copy((vect.Sparse)->begin(), (vect.Sparse)->end(), (svect.sparse())->begin());
  }

  // Note FP: using constructor + noalias = (or std::copy) is more
  // efficient than a call to ublas::vector copy constructor, this for
  // large or small vectors.
}

// Copy
SimpleVector::SimpleVector(const SiconosVector &v): SiconosVector() // Dense by default
{
  unsigned int numV = v.getNum();
  if (numV == 0) // ie if v is a BlockVector, "this" is set to dense.
  {
    // num = 1; default value
    _dense = true;
    vect.Dense = new DenseVect(v.size());

    VectorOfVectors::const_iterator it;
    unsigned int pos = 0;
    for (it = v.begin(); it != v.end(); ++it)
    {
      setBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else if (numV == 1) // dense
  {
    _dense = true;
    vect.Dense = new DenseVect(v.size());
    noalias(*vect.Dense) = (*v.dense());
    //std::copy((v.dense())->begin(), (v.dense())->end(), (vect.Dense)->begin());
  }
  else // sparse
  {
    _dense = false;
    vect.Sparse = new SparseVect(v.size());
    std::copy((v.sparse())->begin(), (v.sparse())->end(), (vect.Sparse)->begin());
  }
}

SimpleVector::SimpleVector(const DenseVect& m): SiconosVector()
{
  _dense = true;
  vect.Dense = new DenseVect(m.size());
  noalias(*vect.Dense) = m;

}

SimpleVector::SimpleVector(const SparseVect& m): SiconosVector()
{
  _dense = false;
  vect.Sparse = new SparseVect(m.size());
  noalias(*vect.Sparse) = m;
}

SimpleVector::SimpleVector(const std::string &file, bool ascii): SiconosVector()
{
  _dense = true;
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

SimpleVector::~SimpleVector()
{
  if (_dense)
    delete(vect.Dense);
  else delete(vect.Sparse);
}


// =================================================
//        get Ublas component (dense or sparse)
// =================================================

const DenseVect SimpleVector::getDense(unsigned int) const
{
  if (!_dense)
    SiconosVectorException::selfThrow("SimpleVector::getDense(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return *vect.Dense;
}

const SparseVect SimpleVector::getSparse(unsigned int)const
{

  if (_dense)
    SiconosVectorException::selfThrow("SimpleVector::getSparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return *vect.Sparse;
}

SparseVect* SimpleVector::sparse(unsigned int)const
{

  if (_dense)
    SiconosVectorException::selfThrow("SimpleVector::sparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return vect.Sparse;
}

double* SimpleVector::getArray(unsigned int) const
{
  if (!_dense)
    SiconosVectorException::selfThrow("SimpleVector::getArray() : not yet implemented for sparse vector.");

  assert(vect.Dense);

  return &(((*vect.Dense).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void SimpleVector::zero()
{
  if (_dense)
    atlas::set(0.0, *vect.Dense);

  else
  {
    assert(vect.Sparse);
    *vect.Sparse *= 0.0;
  }

}

void SimpleVector::setVector(unsigned int, const SiconosVector& newV)
{
  if (newV.size() != size())
    SiconosVectorException::selfThrow("SimpleVector::setVector(num,v), unconsistent sizes.");

  *this = newV ;
}

void SimpleVector::fill(double value)
{
  if (!_dense)
  {
    for (unsigned int i = 0; i < (vect.Sparse)->size(); ++i)
      (vect.Sparse)->push_back(i, value);
  }
  else
    atlas::set(value, *vect.Dense);
}

//=======================
// set vector dimension
//=======================

void SimpleVector::resize(unsigned int n, bool preserve)
{
  if (_dense)
    (vect.Dense)->resize(n, preserve);
  else
    (vect.Sparse)->resize(n, preserve);
}

//=======================
//       get norm
//=======================

double SimpleVector::normInf() const
{
  if (_dense)
    return norm_inf(*vect.Dense);
  else //if(num==4)
    return norm_inf(*vect.Sparse);
}

double SimpleVector::norm2() const
{
  if (_dense)
    return ublas::norm_2(*vect.Dense);
  else //if(num==4)
    return ublas::norm_2(*vect.Sparse);
}
//======================================
// get sum of all elements of the vector
//=====================================
double SimpleVector::sum() const
{
  if (_dense)
    return ublas::sum(*vect.Dense);
  else
    return ublas::sum(*vect.Sparse);
}

//=====================
// screen display
//=====================

void SimpleVector::display()const
{
  if (_dense)
    std::cout << *vect.Dense << std::endl;
  else if (vect.Sparse)
    std::cout << *vect.Sparse << std::endl;
}

//============================
// Convert vector to a string
//============================

const std::string SimpleVector::toString() const
{
  std::stringstream sstr;
  std::string s;
  if (_dense)
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

//=============================
// Elements access (get or set)
//=============================

double SimpleVector::getValue(unsigned int row) const
{
  if (row >= size())
    SiconosVectorException::selfThrow("SimpleVector::getValue(index) : Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row);
}

void SimpleVector::setValue(unsigned int row, double value)
{
  if (row >= size())
    SiconosVectorException::selfThrow("SimpleVector::setValue(index, value) : Index out of range");
  if (_dense)
    (*vect.Dense)(row) = value ;
  else
    (*vect.Sparse)(row) = value;
}

double& SimpleVector::operator()(unsigned int row)
{
  if (row >= size())
    SiconosVectorException::selfThrow("SimpleVector::operator ( index ): Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

double SimpleVector::operator()(unsigned int row) const
{
  if (row >= size())
    SiconosVectorException::selfThrow("SimpleVector::operator ( index ): Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return ((*vect.Sparse)(row)).ref();
}

//============================================
// Access (get or set) to blocks of elements
//============================================

void SimpleVector::setBlock(unsigned int index, const SiconosVector& vIn)
{
  // Set current vector elements, starting from position "index", to the values of vector vIn
  // vIn may be a BlockVector.

  // Exceptions ...
  if (&vIn == this)
    SiconosVectorException::selfThrow("SimpleVector::this->setBlock(pos,vIn): vIn = this.");

  if (index > size())
    SiconosVectorException::selfThrow("SimpleVector::setBlock : invalid ranges");

  unsigned int end = vIn.size() + index;
  if (end > size())
    SiconosVectorException::selfThrow("SimpleVector::setBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();
  if (numVin == 0) // if vIn is a BlockVector
  {
    VectorOfVectors::const_iterator it;
    unsigned int pos = index;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      setBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else // vIn a SimpleVector ...
  {
    if (numVin != getNum()) SiconosVectorException::selfThrow("SimpleVector::setBlock: inconsistent types.");

    if (_dense)
      noalias(ublas::subrange(*vect.Dense, index, end)) = *vIn.dense();
    else
      noalias(ublas::subrange(*vect.Sparse, index, end)) = *vIn.sparse();
  }
}

void setBlock(const SiconosVector& vIn, SP::SiconosVector vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut)
{
  // To copy a subBlock of vIn (from position startIn to startIn+sizeB) into vOut (from pos. startOut to startOut+sizeB).
  if (&vIn == vOut.get()) // useless op. => nothing to be done
  {}// SiconosVectorException::selfThrow("");
  else
  {
    // Check dim ...
    unsigned int sizeIn = vIn.size();
    unsigned int sizeOut = vOut->size();

    if (startIn >= sizeIn)
      SiconosVectorException::selfThrow("vector setBlock(v1,v2,...): start position in input vector is out of range.");

    if (startOut >= sizeOut)
      SiconosVectorException::selfThrow("vector setBlock(v1,v2,...): start position in output vector is out of range.");

    unsigned int endIn = startIn + sizeB;
    unsigned int endOut = startOut + sizeB;

    if (endIn > sizeIn)
      SiconosVectorException::selfThrow("vector setBlock(v1,v2,...): end position in input vector is out of range.");
    if (endOut > sizeOut)
      SiconosVectorException::selfThrow("vector setBlock(v1,v2,...): end position in output vector is out of range.");

    unsigned int numIn = vIn.getNum();
    unsigned int numOut = vOut->getNum();

    if (numIn == 0 && numOut == 0)
      SiconosVectorException::selfThrow("vector setBlock(v1,v2,...): not yet implemented for v1 and v2 both BlockVectors. Try to use setBlock on the sub-vectors?");
    else if (numOut == 0) // vOut is block ...
    {
      // We look for the block of vOut that include index startOut
      unsigned int blockOutStart = 0;
      const SP::Index tabOut = vOut->tabIndex();
      while (startOut >= (*tabOut)[blockOutStart] && blockOutStart < tabOut->size())
        blockOutStart ++;
      // Relative position in the block blockOutStart.
      unsigned int posOut = startOut;
      if (blockOutStart != 0)
        posOut -= (*tabOut)[blockOutStart - 1];

      // We look for the block of vOut that include index endOut
      unsigned int blockOutEnd = blockOutStart;
      while (endOut > (*tabOut)[blockOutEnd] && blockOutEnd < tabOut->size())
        blockOutEnd ++;

      // => the block to be set runs from block number blockOutStart to block number blockOutEnd.

      if (blockOutEnd == blockOutStart) //
      {
        setBlock(vIn, vOut->vector(blockOutStart), sizeB, startIn, posOut);
      }
      else // More that one block of vOut are concerned
      {

        // The current considered block ...
        SP::SiconosVector currentBlock = vOut->vector(blockOutStart);

        // Size of the subBlock of vOut to be set.
        unsigned int subSizeB = currentBlock->size() - posOut;
        unsigned int posIn = startIn;

        // Set first sub-block (currentBlock) values, between index posOut and posOut+subSizeB,
        // with vIn values from posIn to posIn+subSizeB.
        setBlock(vIn, currentBlock, subSizeB, posIn, posOut);

        // Other blocks, except number blockOutEnd.
        unsigned int currentBlockNum = blockOutStart + 1;
        while (currentBlockNum != blockOutEnd)
        {
          posIn += subSizeB;
          currentBlock =  vOut->vector(currentBlockNum);
          subSizeB = currentBlock->size();
          setBlock(vIn, currentBlock, subSizeB, posIn, 0);
          currentBlockNum++;
        }
        // set last subBlock ...
        currentBlock =  vOut->vector(blockOutEnd);

        posIn += subSizeB;

        // Size of the considered sub-block
        subSizeB = endOut - (*tabOut)[blockOutEnd - 1];

        setBlock(vIn, currentBlock, subSizeB, posIn, 0);
      }
    }
    else if (numIn == 0) // vIn is block
    {

      // We look for the block of vIn that include index startIn
      unsigned int blockInStart = 0;
      const SP::Index tabIn = vIn.tabIndex();
      while (startIn >= (*tabIn)[blockInStart] && blockInStart < tabIn->size())
        blockInStart ++;
      // Relative position in the block blockInStart.
      unsigned int posIn = startIn;
      if (blockInStart != 0)
        posIn -= (*tabIn)[blockInStart - 1];

      // We look for the block of vIn that include index endIn
      unsigned int blockInEnd = blockInStart;
      while (endIn > (*tabIn)[blockInEnd] && blockInEnd < tabIn->size())
        blockInEnd ++;

      // => the block to be set runs from block number blockInStart to block number blockInEnd.

      if (blockInEnd == blockInStart) //
      {
        setBlock(*vIn.vector(blockInStart), vOut, sizeB, posIn, startOut);
      }
      else // More that one block of vIn are concerned
      {

        // The current considered block ...
        SPC::SiconosVector currentBlock = vIn.vector(blockInStart);

        // Size of the subBlock of vIn to be set.
        unsigned int subSizeB = currentBlock->size() - posIn;
        unsigned int posOut = startOut;

        // Set vOut values, between index posOut and posOut+subSizeB,
        // with first sub-block (currentBlock) of vIn values from posIn to posIn+subSizeB.
        setBlock(*currentBlock, vOut, subSizeB, posIn, posOut);

        // Other blocks, except number blockInEnd.
        unsigned int currentBlockNum = blockInStart + 1;
        while (currentBlockNum != blockInEnd)
        {
          posOut += subSizeB;
          currentBlock =  vIn.vector(currentBlockNum);
          subSizeB = currentBlock->size();
          setBlock(*currentBlock, vOut, subSizeB, 0, posOut);
          currentBlockNum++;
        }
        // set last subBlock ...
        currentBlock =  vIn.vector(blockInEnd);

        posOut += subSizeB;

        // Relative position of index endIn in vIn[blockInEnd]
        subSizeB = endIn - (*tabIn)[blockInEnd - 1];

        setBlock(*currentBlock, vOut, subSizeB, 0, posOut);
      }
    }
    else // neither vIn nor vOut is a BlockVector
    {

      if (numIn == numOut)
      {
        if (numIn == 1) // vIn / vOut are Dense
          noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) = ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
        else // if(numIn == 4)// vIn / vOut are Sparse
          noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) = ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
      }
      else // vIn and vout of different types ...
      {
        if (numIn == 1) // vIn Dense
          noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) = ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
        else // if(numIn == 4)// vIn Sparse
          noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) = ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
      }
    }
  }
}

void SimpleVector::addBlock(unsigned int index, const SiconosVector& vIn)
{
  // Add vIn to the current vector, starting from position "index".
  // vIn may be a BlockVector.

  //if ( num != 1 ) SiconosVectorException::selfThrow("SimpleVector::addBlock : vector should be dense");

  if (&vIn == this)
    SiconosVectorException::selfThrow("SimpleVector::this->addBlock(pos,vIn): vIn = this.");

  unsigned int end = vIn.size();
  if ((index + end) > size()) SiconosVectorException::selfThrow("SimpleVector::addBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();
  if (numVin == 0) // if vIn is a BlockVector
  {
    VectorOfVectors::const_iterator it;
    unsigned int pos = index;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      addBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else // vIn is a SimpleVector ...
  {
    if (numVin != getNum()) SiconosVectorException::selfThrow("SimpleVector::addBlock : inconsistent types.");

    if (_dense)
      noalias(ublas::subrange(*vect.Dense, index, index + end)) += *vIn.dense();
    else
      noalias(ublas::subrange(*vect.Sparse, index, index + end)) += *vIn.sparse();
  }
}

void SimpleVector::subBlock(unsigned int index, const SiconosVector& vIn)
{
  // Add vIn from the current vector, starting from position "index".
  // vIn may be a BlockVector.

  //  if ( num != 1 ) SiconosVectorException::selfThrow("SimpleVector::subBlock : vector should be dense");

  unsigned int end = vIn.size();
  if ((index + end) > size()) SiconosVectorException::selfThrow("SimpleVector::subBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();
  if (numVin == 0) // if vIn is a BlockVector
  {
    VectorOfVectors::const_iterator it;
    unsigned int pos = index;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      subBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else //  vIn a SimpleVector ...
  {
    if (numVin != getNum()) SiconosVectorException::selfThrow("SimpleVector::subBlock : inconsistent types.");

    if (_dense)
      noalias(ublas::subrange(*vect.Dense, index, index + end)) -= *vIn.dense();
    else
      noalias(ublas::subrange(*vect.Sparse, index, index + end)) -= *vIn.sparse();
  }
}

//===============
//  Assignment
//===============

SimpleVector& SimpleVector::operator = (const SiconosVector& vIn)
{
  if (&vIn == this) return *this; // auto-assignment.

  if (size() != vIn.size())
    SiconosVectorException::selfThrow("SimpleVector::operator = failed: inconsistent sizes.");

  unsigned int vInNum = vIn.getNum();

  if (vInNum == 0) // if vIn is a BlockVector
  {
    VectorOfVectors::const_iterator it;
    unsigned int pos = 0;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      setBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else // vIn is a SimpleVector
  {
    switch (getNum())
    {
    case 1:
      switch (vInNum)
      {
      case 1:
        //atlas::copy(*vIn.dense(),*vect.Dense);
        noalias(*vect.Dense) = *vIn.dense();
        break;
      case 4:
        noalias(*vect.Dense) = *vIn.sparse();
        break;
      default:
        SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
        break;
      }
      break;
    case 4:
      if (vInNum == 4)
        noalias(*vect.Sparse) = *vIn.sparse();
      else
        SiconosVectorException::selfThrow("SimpleVector::operator = : can not set sparse = dense.");
      break;
    default:
      SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
      break;
    }
  }
  return *this;
}

SimpleVector& SimpleVector::operator = (const SimpleVector& vIn)
{
  if (&vIn == this) return *this; // auto-assignment.

  if (size() != vIn.size())
    SiconosVectorException::selfThrow("SimpleVector::operator = failed: inconsistent sizes.");

  unsigned int vInNum = vIn.getNum();
  switch (getNum())
  {
  case 1:
    switch (vInNum)
    {
    case 1:
      //atlas::copy(*vIn.dense(),*vect.Dense);
      noalias(*vect.Dense) = *vIn.dense();

      break;
    case 4:
      noalias(*vect.Dense) = *vIn.sparse();
      break;
    default:
      SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
      break;
    }
    break;
  case 4:
    if (vInNum == 4)
      noalias(*vect.Sparse) = *vIn.sparse();
    else
      SiconosVectorException::selfThrow("SimpleVector::operator = : can not set sparse = dense.");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator = (const DenseVect& d)
{
  if (!_dense)
    SiconosVectorException::selfThrow("SimpleVector::operator = DenseVect : forbidden: the current vector is not dense.");
  if (d.size() != size())
    SiconosVectorException::selfThrow("SimpleVector::operator = DenseVect : inconsistent size.");

  atlas::copy(d, *vect.Dense);
  return *this;
}

SimpleVector& SimpleVector::operator = (const SparseVect& sp)
{
  if (_dense)
    SiconosVectorException::selfThrow("SimpleVector::operator = SparseVect : current vector is not sparse.");
  if (sp.size() != size())
    SiconosVectorException::selfThrow("SimpleVector::operator = SparseVect : inconsistent size.");

  noalias(*vect.Sparse) = sp;

  return *this;
}

//=================================
// Op. and assignment (+=, -= ... )
//=================================

SimpleVector& SimpleVector::operator += (const SiconosVector& vIn)
{
  if (&vIn == this) // alias
  {
    // Note: using this *= 2.0 is much more time-consuming.
    switch (getNum())
    {
    case 1:
      *vect.Dense += *vect.Dense;
      break;
    case 4:
      *vect.Sparse += *vect.Sparse;
      break;
    default:
      SiconosVectorException::selfThrow("SimpleVector::operator += : invalid type given");
      break;
    }
    return *this;
  }

  unsigned int vInNum = vIn.getNum();
  if (vInNum == 0) // vIn block
  {
    if (size() != vIn.size())
      SiconosVectorException::selfThrow("SimpleVector::operator += failed: inconsistent sizes.");

    VectorOfVectors::const_iterator it;
    unsigned int pos = 0;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      addBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else  // vIn Simple
  {
    switch (getNum())
    {
    case 1:
      switch (vInNum)
      {
      case 1:
        noalias(*vect.Dense) += *vIn.dense();
        break;
      case 4:
        noalias(*vect.Dense) += *vIn.sparse();
        break;
      default:
        SiconosVectorException::selfThrow("SimpleVector::operator += : invalid type given");
        break;
      }
      break;
    case 4:
      if (vInNum == 4)
        noalias(*vect.Sparse) += *vIn.sparse();
      else SiconosVectorException::selfThrow("SimpleVector::operator += : can not add a dense to a sparse.");
      break;
    default:
      SiconosVectorException::selfThrow("SimpleVector::operator += : invalid type given");
      break;
    }
  }
  return *this;
}

SimpleVector& SimpleVector::operator -= (const SiconosVector& vIn)
{
  if (&vIn == this)
  {
    this->zero();
    return *this;
  }

  unsigned int vInNum = vIn.getNum();
  if (vInNum == 0) // vIn block
  {
    if (size() != vIn.size())
      SiconosVectorException::selfThrow("SimpleVector::operator -= failed: inconsistent sizes.");

    VectorOfVectors::const_iterator it;
    unsigned int pos = 0;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      subBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else // if vIn is not a BlockVector
  {
    unsigned int vInNum = vIn.getNum();
    switch (getNum())
    {
    case 1:
      switch (vInNum)
      {
      case 1:
        noalias(*vect.Dense) -= *vIn.dense();
        break;
      case 4:
        noalias(*vect.Dense) -= *vIn.sparse();
        break;
      default:
        SiconosVectorException::selfThrow("SimpleVector::operator -= : invalid type given");
        break;
      }
      break;
    case 4:
      if (vInNum == 4)
        noalias(*vect.Sparse) -= *vIn.sparse();
      else SiconosVectorException::selfThrow("SimpleVector::operator -= : can not sub a dense to a sparse.");
      break;
    default:
      SiconosVectorException::selfThrow("SimpleVector::operator -= : invalid type given");
      break;
    }
  }
  return *this;
}

//===============
// Comparison
//===============

bool operator == (const SiconosVector &m, const SiconosVector &x)
{
  assert(m.isBlock() == false && x.isBlock() == false);
  return ((m - x).norm2() < tolerance);
}

//==================
// y = scalar * x
//==================

SimpleVector operator * (const  SiconosVector&m, double d)
{
  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block
  {
    SimpleVector tmp(m); // copy ...
    tmp *= d;
    return tmp;
  }
  else if (numM == 1)
  {
    // Copy m into p and call atlas::scal(d,p), p = d*p.
    DenseVect p = *m.dense();
    atlas::scal(d, p);
    return p;
  }
  else// if(numM==4)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

SimpleVector operator * (double d, const  SiconosVector&m)
{
  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block
  {
    SimpleVector tmp(m); // copy ...
    tmp *= d;
    return tmp;
  }
  else if (numM == 1)
  {
    // Copy m into p and call atlas::scal(d,p), p = d*p.
    DenseVect p = *m.dense();
    atlas::scal(d, p);
    return p;
  }
  else// if(numM==4)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

SimpleVector operator / (const SimpleVector &m, double d)
{
  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block
  {
    SimpleVector tmp(m); // copy ...
    tmp /= d;
    return tmp;
  }

  else if (numM == 1)
  {
    DenseVect p = *m.dense();
    atlas::scal((1.0 / d), p);
    return p;
  }

  else// if(numM==4){
    return (SparseVect)(*m.sparse() / d);
}

//====================
//  Vectors addition
//====================

SimpleVector operator + (const  SiconosVector& x, const  SiconosVector& y)
{
  if (x.size() != y.size())
    SiconosVectorException::selfThrow("SiconosVector, x + y: inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY && numX != 0) // x, y SimpleVector of the same type
  {
    if (numX == 1)
    {
      //    atlas::xpy(*x.dense(),p);
      //    return p;
      return (DenseVect)(*x.dense() + *y.dense());
    }
    else
      return (SparseVect)(*x.sparse() + *y.sparse());
  }

  else if (numX != 0 && numY != 0  && numX != numY) // x, y SimpleVector with y and x of different types
  {
    if (numX == 1)
      return (DenseVect)(*x.dense() + *y.sparse());
    else
      return (DenseVect)(*x.sparse() + *y.dense());
  }

  else if (numX == 0) // if x is block and y simple or block
  {
    SimpleVector tmp(y);
    tmp += x;
    return tmp;
  }

  else // if y is block and x simple
  {
    SimpleVector tmp(x);
    tmp += y;
    return tmp;
  }
}

void add(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  // Computes z = x + y in an "optimized" way (in comparison with operator +)

  if (x.size() != y.size() || x.size() != z.size())
    SiconosVectorException::selfThrow("add(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  unsigned int numZ = z.getNum();

  if (&z == &x) // x, and z are the same object.
  {
    z += y;
  }
  else if (&z == &y) // y and z are the same object, different from x
  {
    z += x;
  }
  else // No common memory between x,y and z
  {

    if (numZ != 0) // z is a SimpleVector
    {
      if (numX == numY && numX != 0) // x, y SimpleVector of the same type
      {
        if (numX == 1)
        {
          if (numZ != 1)
            SiconosVectorException::selfThrow("SiconosVector addition, add(x,y,z) failed - Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() + *y.dense() ;
        }
        else
        {
          if (numZ == 1)
            noalias(*z.dense()) = *x.sparse() + *y.sparse() ;
          else
            noalias(*z.sparse()) = *x.sparse() + *y.sparse() ;
        }
      }
      else if (numX != 0 && numY != 0) // x and y of different types => z must be dense.
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector addition, add(x,y,z) failed - z can not be sparse.");
        if (numX == 1)
          noalias(*z.dense()) = *x.dense() + *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() + *y.dense() ;
      }
      else if (numX == 0) // y simple or block
      {
        z = y;
        z += x;
      }
      else // x simple, y block
      {
        z = x;
        z += y;
      }
    }
    else // z is a BlockVector
    {
      if (numX != 0)
      {
        z = x;
        z += y;
      }
      else
      {
        z = y;
        z += x;
      }
    }
  }
}

//======================
//  Vectors subtraction
//======================

SimpleVector operator - (const  SiconosVector& x, const  SiconosVector& y)
{
  if (x.size() != y.size())
    SiconosVectorException::selfThrow("SiconosVector, x - y: inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY && numX != 0) // x, y SimpleVector of the same type
  {
    if (numX == 1)
    {
      //    atlas::xpy(*x.dense(),p);
      //    return p;
      return (DenseVect)(*x.dense() - *y.dense());
    }
    else
      return (SparseVect)(*x.sparse() - *y.sparse());
  }

  else if (numX != 0 && numY != 0  && numX != numY) // x, y SimpleVector with y and x of different types
  {
    if (numX == 1)
      return (DenseVect)(*x.dense() - *y.sparse());
    else
      return (DenseVect)(*x.sparse() - *y.dense());
  }

  else // if y OR x is block
  {
    SimpleVector tmp(x);
    tmp -= y;
    return tmp;
  }
}

void sub(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  // Computes z = x - y in an "optimized" way (in comparison with operator +)

  if (x.size() != y.size() || x.size() != z.size())
    SiconosVectorException::selfThrow("sub(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  unsigned int numZ = z.getNum();

  if (&z == &x) // x and z are the same object.
  {
    z -= y;
  }
  else if (&z == &y) // y and z are the same object
  {
    if (numZ == 0 ||  numX == 0) // z or x is block
    {
      z *= -1.0; // ...
      z += x ;
    }
    else
    {
      if (numX == 1)
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector subtraction, sub(x,y,z) failed - Subtraction of two dense vectors into a sparse.");
        *z.dense() = *x.dense() - *y.dense() ;
      }
      else
      {
        if (numZ == 1)
          *z.dense() = *x.sparse() - *y.dense() ;
        else
          *z.sparse() = *x.sparse() - *y.sparse() ;
      }
    }
  }
  else // No common memory between x or y and z
  {

    if (numZ != 0) // z is a SimpleVector
    {
      if (numX == numY && numX != 0) // x, y SimpleVector of the same type
      {
        if (numX == 1)
        {
          if (numZ != 1)
            SiconosVectorException::selfThrow("SiconosVector addition, sub(x,y,z) failed - Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() - *y.dense() ;
        }
        else
        {
          if (numZ == 1)
            noalias(*z.dense()) = *x.sparse() - *y.sparse() ;
          else
            noalias(*z.sparse()) = *x.sparse() - *y.sparse() ;
        }
      }
      else if (numX != 0 && numY != 0) // x and y of different types => z must be dense.
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector addition, sub(x,y,z) failed - z can not be sparse.");
        if (numX == 1)
          noalias(*z.dense()) = *x.dense() - *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() - *y.dense() ;
      }
      else // x simple, y block
      {
        z = x;
        z -= y;
      }
    }
    else // z is a BlockVector
    {
      z = x;
      z -= y;
    }
  }
}

void axpby(double a, const SiconosVector& x, double b, SiconosVector& y)
{
  // Computes y = ax + by

  if (x.size() != y.size())
    SiconosVectorException::selfThrow("axpby(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY) // x and y of the same type
  {
    if (numX == 1) // all dense
      atlas::axpby(a, *x.dense(), b, *y.dense());

    else if (numX == 0) // ie if x and y are block
    {
      if (isComparableTo(x, y))
      {
        VectorOfVectors::iterator itY;
        VectorOfVectors::const_iterator itX = x.begin();
        for (itY = y.begin(); itY != y.end(); ++itY)
          axpby(a, **itX++, b, **itY);
      }
      else
      {
        // bad case ...
        y *= b;
        SimpleVector tmp(x);
        tmp *= a;
        y += tmp;
      }
    }

    else // all sparse
    {
      *y.sparse() *= b;
      if (&y != &x)
        noalias(*y.sparse()) += a**x.sparse();
      else
        *y.sparse() += a**x.sparse();
    }
  }

  else // x and y of different types
  {
    y *= b;
    if (numY == 0) // y  block, x simple
    {
      SimpleVector tmp(x);
      tmp *= a;
      y += tmp;
    }
    else if (numX == 0) // x block, y simple
    {
      SimpleVector tmp(x);
      if (numY == 1)
        atlas::axpy(a, *tmp.dense(), *y.dense());
      else
        SiconosVectorException::selfThrow("axpby failed: try to add block to sparse vector.");
    }
    else // x and y simple but of different types
    {
      if (numX == 1)
        *y.sparse() += a**x.dense();
      else
        *y.dense() +=  a**x.sparse();
    }
  }
}

void axpy(double a, const SiconosVector& x, SiconosVector& y)
{
  // Computes y = ax + y

  if (x.size() != y.size())
    SiconosVectorException::selfThrow("axpy(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY) // x and y of the same type
  {
    if (numX == 1) // all dense
      atlas::axpy(a, *x.dense(), *y.dense());

    else if (numX == 0) // ie if x and y are block
    {
      if (isComparableTo(x, y))
      {
        VectorOfVectors::iterator itY;
        VectorOfVectors::const_iterator itX = x.begin();
        for (itY = y.begin(); itY != y.end(); ++itY)
          axpy(a, **itX++, **itY);
      }
      else
      {
        // bad case ...
        SimpleVector tmp(x);
        tmp *= a;
        y += tmp;
      }
    }

    else // all sparse
    {
      if (&y != &x)
        noalias(*y.sparse()) += a**x.sparse();
      else
        *y.sparse() += a**x.sparse();
    }
  }

  else // x and y of different types
  {
    if (numY == 0) // y  block, x simple
    {
      SimpleVector tmp(x);
      tmp *= a;
      y += tmp;
    }
    else if (numX == 0) // x block, y simple
    {
      SimpleVector tmp(x);
      if (numY == 1)
        atlas::axpy(a, *tmp.dense(), *y.dense());
      else
        SiconosVectorException::selfThrow("axpby failed: try to add block to sparse vector.");
    }
    else // x and y simple but of different types
    {
      if (numX == 1)
        *y.sparse() += a**x.dense();
      else
        *y.dense() +=  a**x.sparse();
    }
  }
}

double inner_prod(const SiconosVector &x, const SiconosVector &m)
{
  if (x.size() != m.size())
    SiconosVectorException::selfThrow("inner_prod: inconsistent sizes");

  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numX == 0 || numM == 0)
    SiconosVectorException::selfThrow("inner_prod: not implemented for BlockVectors.");

  if (numX == numM)
  {
    if (numM == 1)
      return atlas::dot(*x.dense(), *m.dense());
    else
      return inner_prod(*x.sparse(), *m.sparse());
  }
  else if (numM == 1)
    return inner_prod(*x.sparse(), *m.dense());
  else
    return inner_prod(*x.dense(), *m.sparse());
}

// outer_prod(v,w) = trans(v)*w
SimpleMatrix outer_prod(const SiconosVector &x, const SiconosVector& m)
{
  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();
  if (numX == 0 || numM == 0)
    SiconosVectorException::selfThrow("outer_prod: not implemented for BlockVectors.");
  if (numM == 1)
  {
    if (numX == 1)
      return (DenseMat)(outer_prod(*x.dense(), *m.dense()));

    else// if(numX == 4)
      return (DenseMat)(outer_prod(*x.sparse(), *m.dense()));
  }
  else // if(numM == 4)
  {
    if (numX == 1)
      return (DenseMat)(outer_prod(*x.dense(), *m.sparse()));

    else //if(numX == 4)
      return (DenseMat)(outer_prod(*x.sparse(), *m.sparse()));
  }
}

void scal(double a, const SiconosVector & x, SiconosVector & y, bool init)
{
  // To compute y = a *x (init = true) or y += a*x (init = false)

  if (&x == &y)
  {
    if (init)
      y *= a;
    else
    {
      y *= (1.0 + a);
    }
  }
  else
  {
    unsigned int sizeX = x.size();
    unsigned int sizeY = y.size();

    if (sizeX != sizeY)
      SiconosVectorException::selfThrow("scal(a,SiconosVector,SiconosVector) failed, sizes are not consistent.");

    unsigned int numY = y.getNum();
    unsigned int numX = x.getNum();
    if (numX == numY)
    {
      if (numX == 0) // ie if both are block vectors
      {
        if (isComparableTo(x, y)) // if x and y are "block-consistent"
        {
          VectorOfVectors::const_iterator itX = x.begin();
          VectorOfVectors::iterator itY ;
          for (itY = y.begin(); itY != y.end(); ++itY)
            scal(a, **itX++, **itY++, init);
        }
        else
        {
          if (init)
          {
            for (unsigned int i = 0; i < x.size(); ++i)
              y(i) = a * x(i);
          }
          else
          {
            for (unsigned int i = 0; i < x.size(); ++i)
              y(i) += a * x(i);
          }
        }
      }
      else if (numX == 1) // ie if both are Dense
      {
        if (init)
          //atlas::axpby(a,*x.dense(),0.0,*y.dense());
          noalias(*y.dense()) = a * *x.dense();
        else
          noalias(*y.dense()) += a * *x.dense();
      }
      else  // if both are sparse
      {
        if (init)
          noalias(*y.sparse()) = a**x.sparse();
        else
          noalias(*y.sparse()) += a**x.sparse();
      }
    }
    else
    {
      if (numY == 0 || numX == 0) // if y or x is block
      {
        if (init)
        {
          y = x;
          y *= a;
        }
        else
        {
          SimpleVector tmp(x);
          tmp *= a;
          y += tmp;
        }
      }
      else
      {
        if (numY == 1) // if y is dense
        {
          if (init)
            noalias(*y.dense()) = a**x.sparse();
          else
            noalias(*y.dense()) += a**x.sparse();

        }
        else
          SiconosVectorException::selfThrow("SiconosVector::scal(a,dense,sparse) not allowed.");
      }
    }
  }
}

void subscal(double a, const SiconosVector & x, SiconosVector & y, const Index& coord, bool init)
{
  // To compute sub_y = a *sub_x (init = true) or sub_y += a*sub_x (init = false)
  // Coord  = [r0x r1x r0y r1y];
  // subX is the sub-vector of x, for row numbers between r0x and r1x-1.
  // The same for y with riy.


  // Check dimensions
  unsigned int dimX = coord[1] - coord[0];
  unsigned int dimY = coord[3] - coord[2];
  if (dimY != dimX)
    SiconosVectorException::selfThrow("subscal(a,x,y,...) error: inconsistent sizes between (sub)x and (sub)y.");
  if (dimY > y.size() || dimX > x.size())
    SiconosVectorException::selfThrow("subscal(a,x,y,...) error: input index too large.");

  unsigned int numY = y.getNum();
  unsigned int numX = x.getNum();

  if (&x == &y) // if x and y are the same object
  {
    if (numX == 0) // if x and y are block vectors
    {
      if (coord[0] != coord[2] || coord[1] != coord[3])
        SiconosVectorException::selfThrow("subscal(a,x,y,...) error: x=y=blockVector and try to affect different positions in x and y; not yet implemented!");

      VectorOfVectors::const_iterator it;
      // Number of the subvector of x that handles element at position coord[0]
      std::size_t firstBlockNum = x.getNumVectorAtPos(coord[0]);
      // Number of the subvector of x that handles element at position coord[1]
      std::size_t lastBlockNum = x.getNumVectorAtPos(coord[1]);
      Index subCoord = coord;
      SP::SiconosVector tmp = y[firstBlockNum];
      std::size_t subSize =  x[firstBlockNum]->size(); // Size of the sub-vector
      const SP::Index xTab = x.tabIndex();
      if (firstBlockNum != 0)
        subCoord[0] -= (*xTab)[firstBlockNum - 1];
      subCoord[1] =  std::min(coord[1] - (*xTab)[firstBlockNum - 1], subSize);
      subCoord[2] = subCoord[0];
      subCoord[3] = subCoord[1];
      if (firstBlockNum == lastBlockNum)
      {
        subscal(a, *tmp, *tmp, subCoord, init);
      }
      else
      {
        unsigned int xPos = 0 ; // Position in x of the current sub-vector of x
        bool firstLoop = true;
        for (it = x.begin(); it != x.end(); ++it)
        {
          if ((*it)->getNum() == 0)
            SiconosVectorException::selfThrow("subscal(a,x,y) error: not yet implemented for x block of blocks ...");
          if (xPos >= firstBlockNum && xPos <= lastBlockNum)
          {
            tmp = y[xPos];
            if (firstLoop)
            {
              subscal(a, *tmp, *tmp, subCoord, init);
              firstLoop = false;
            }
            else
            {
              subSize = tmp->size();
              subCoord[0] = 0;
              subCoord[1] = std::min(coord[1] - (*xTab)[xPos - 1], subSize);
              subCoord[2] = subCoord[0];
              subCoord[3] = subCoord[1];
              subscal(a, *tmp, *tmp, subCoord, init);
            }
          }
          xPos++;
        }
      }
    }

    else if (numX == 1) // Dense
    {
      ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2])
      {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
    else //if (numX == 4) // Sparse
    {
      ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2])
      {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
  }
  else
  {
    if (numX == numY)
    {
      if (numX == 0) // ie if both are block vectors
      {
        SiconosVectorException::selfThrow("subscal(a,x,y) error: not yet implemented for x and y block vectors");
      }
      else if (numX == 1) // ie if both are Dense
      {
        ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
      else  // if both are sparse
      {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
    }
    else // x and y of different types ...
    {
      if (numX == 0) // x a block vector
      {
        VectorOfVectors::const_iterator it;
        // Number of the subvector of x that handles element at position coord[0]
        std::size_t firstBlockNum = x.getNumVectorAtPos(coord[0]);
        // Number of the subvector of x that handles element at position coord[1]
        std::size_t lastBlockNum = x.getNumVectorAtPos(coord[1]);
        Index subCoord = coord;

        SPC::SiconosVector tmp = x[firstBlockNum];

        std::size_t subSize =  x[firstBlockNum]->size(); // Size of the sub-vector
        const SP::Index xTab = x.tabIndex();
        if (firstBlockNum != 0)
          subCoord[0] -= (*xTab)[firstBlockNum - 1];
        subCoord[1] =  std::min(coord[1] - (*xTab)[firstBlockNum - 1], subSize);
        subCoord[3] = subCoord[2] + subCoord[1] - subCoord[0];
        if (firstBlockNum == lastBlockNum)
        {
          subscal(a, *tmp, y, subCoord, init);
        }
        else
        {
          std::size_t xPos = 0 ; // Position in x of the current sub-vector of x
          bool firstLoop = true;
          for (it = x.begin(); it != x.end(); ++it)
          {
            if ((*it)->getNum() == 0)
              SiconosVectorException::selfThrow("subscal(a,x,y) error: not yet implemented for x block of blocks ...");
            if (xPos >= firstBlockNum && xPos <= lastBlockNum)
            {
              tmp = x[xPos];
              if (firstLoop)
              {
                subscal(a, *tmp, y, subCoord, init);
                firstLoop = false;
              }
              else
              {
                subSize = tmp->size();
                subCoord[0] = 0;
                subCoord[1] = std::min(coord[1] - (*xTab)[xPos - 1], subSize);
                subCoord[2] = subCoord[3];
                subCoord[3] = subCoord[2] + subCoord[1] - subCoord[0];
                subscal(a, *tmp, y, subCoord, init);
              }
            }
            xPos++;
          }
        }
      }

      else if (numY == 0) // y a block vector
      {
        VectorOfVectors::const_iterator it;
        // Number of the subvector of y that handles element at position coord[2]
        std::size_t firstBlockNum = y.getNumVectorAtPos(coord[2]);
        // Number of the subvector of x that handles element at position coord[3]
        std::size_t lastBlockNum = y.getNumVectorAtPos(coord[3]);
        Index subCoord = coord;
        SP::SiconosVector tmp = y[firstBlockNum];
        std::size_t subSize =  y[firstBlockNum]->size(); // Size of the sub-vector
        const SP::Index yTab = y.tabIndex();
        if (firstBlockNum != 0)
          subCoord[2] -= (*yTab)[firstBlockNum - 1];
        subCoord[3] =  std::min(coord[3] - (*yTab)[firstBlockNum - 1], subSize);
        subCoord[1] = subCoord[0] + subCoord[3] - subCoord[2];
        if (firstBlockNum == lastBlockNum)
        {
          subscal(a, x, *tmp, subCoord, init);
        }
        else
        {
          std::size_t yPos = 0 ; // Position in x of the current sub-vector of x
          bool firstLoop = true;
          for (it = y.begin(); it != y.end(); ++it)
          {
            if ((*it)->getNum() == 0)
              SiconosVectorException::selfThrow("subscal(a,x,y) error: not yet implemented for y block of blocks ...");
            if (yPos >= firstBlockNum && yPos <= lastBlockNum)
            {
              tmp = y[yPos];
              if (firstLoop)
              {
                subscal(a, x, *tmp, subCoord, init);
                firstLoop = false;
              }
              else
              {
                subSize = tmp->size();
                subCoord[2] = 0;
                subCoord[3] = std::min(coord[3] - (*yTab)[yPos - 1], subSize);
                subCoord[0] = subCoord[1];
                subCoord[1] = subCoord[0] + subCoord[3] - subCoord[2];
                subscal(a, x, *tmp, subCoord, init);
              }
            }
            yPos++;
          }
        }
      }
      else if (numY == 1) // y dense, x sparse
      {
        ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
      else // y sparse, x dense => fails
        SiconosVectorException::selfThrow("SiconosVector::subscal(a,dense,sparse) not allowed.");
    }
  }
}
void cross_product(const SiconosVector& V1, const SiconosVector& V2, SiconosVector& VOUT)
{
  if (V1.size() != 3 || V2.size() != 3 || VOUT.size() != 3)
    SiconosVectorException::selfThrow("SiconosVector::cross_product allowed only with dim 3.");

  double aux = V1.getValue(1) * V2.getValue(2) - V1.getValue(2) * V2.getValue(1);
  VOUT.setValue(0, aux);

  aux = V1.getValue(2) * V2.getValue(0) - V1.getValue(0) * V2.getValue(2);
  VOUT.setValue(1, aux);

  aux = V1.getValue(0) * V2.getValue(1) - V1.getValue(1) * V2.getValue(0);
  VOUT.setValue(2, aux);

}

//

void abs_wise(const SiconosVector& V, SiconosVector& Vabs)
{
  for (unsigned int it = 0; it < V.size(); ++it)
  {
    Vabs.setValue(it, std::abs(V.getValue(it)));
  };
}

//

void getMax(const SiconosVector& V, double& maxvalue, unsigned int& idmax)
{
  maxvalue = V.getValue(0);
  idmax = 0;
  for (unsigned int it = 1; it < V.size(); ++it)
  {
    if (V.getValue(it) > maxvalue)
    {
      maxvalue = V.getValue(it);
      idmax = it;
    };
  };
}

//

void getMin(const SiconosVector& V, double& minvalue, unsigned int& idmin)
{
  minvalue = V.getValue(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it)
  {
    if (V.getValue(it) < minvalue)
    {
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}

//
/*
SimpleVector abs_wise(const SimpleVector& V){
  SimpleVector Vabs(V.size());
  for (int it = 0; it < V.size(); ++it){
    Vabs.setValue(it,std::abs(V.getValue(it)));
  };
  return Vabs;
}
//
void getMax(const SimpleVector& V, double& maxvalue, unsigned int& idmax){
  maxvalue = V.getValue(0);
  idmax = 0;
  for (unsigned int it = 1; it < V.size(); ++it){
    if (V.getValue(it) > maxvalue){
    maxvalue = V.getValue(it);
    idmax = it;
    };
  };
}
//
void getMin(const SimpleVector& V, double& minvalue, unsigned int& idmin){
  minvalue = V.getValue(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it){
    if (V.getValue(it) < minvalue){
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}
*/
