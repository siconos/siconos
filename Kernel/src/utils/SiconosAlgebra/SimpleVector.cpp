/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include <boost/numeric/ublas/io.hpp>            // for >> 
#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <boost/numeric/bindings/atlas/cblas1.hpp>

// =================================================
//                CONSTRUCTORS
// =================================================

// parameters: dimension and type.
SimpleVector::SimpleVector(unsigned int row, UBLAS_TYPE typ): SiconosVector(1, row)
{
  if (typ == SPARSE)
  {
    vect.Sparse = new SparseVect(ublas::zero_vector<double>(row));
    num = 4;
  }
  else if (typ == DENSE)
  {
    vect.Dense = new DenseVect(ublas::zero_vector<double>(row));
    // num = 1; // set by default
  }
  else
  {
    SiconosVectorException::selfThrow("SimpleVector::constructor(UBLAS_TYPE, unsigned int) failed, invalid type given");
  }
}

// parameters: dimension, default value for all components and type.
SimpleVector::SimpleVector(unsigned int row, double val, UBLAS_TYPE typ): SiconosVector(1, row)
{
  if (typ == SPARSE)
  {
    vect.Sparse = new SparseVect(row);
    num = 4;
    fill(val);
  }
  else if (typ == DENSE)
  {
    vect.Dense = new DenseVect(ublas::scalar_vector<double>(row, val));
    // num = 1; // set by default
  }
  else
  {
    SiconosVectorException::selfThrow("SimpleVector::constructor(UBLAS_TYPE, unsigned int) : invalid type given");
  }
}

// parameters: a vector (stl) of double and the type.
SimpleVector::SimpleVector(const std::vector<double>& v, UBLAS_TYPE typ): SiconosVector(1, v.size())
{
  if (typ != DENSE)
    SiconosVectorException::selfThrow("SimpleVector::constructor(UBLAS_TYPE, std::vector<double>, unsigned int) : invalid type given");

  vect.Dense = new DenseVect(v.size());
  std::copy(v.begin(), v.end(), (vect.Dense)->begin());
}

// Copy
SimpleVector::SimpleVector(const SimpleVector &svect): SiconosVector(svect.getNum(), svect.size())
{
  if (num == 1) // dense
  {
    vect.Dense = new DenseVect(svect.size());
    noalias(*vect.Dense) = (*svect.getDensePtr());
    // std::copy((vect.Dense)->begin(), (vect.Dense)->end(), (svect.getDensePtr())->begin());
  }
  else //sparse
  {
    vect.Sparse = new SparseVect(svect.size());
    noalias(*vect.Sparse) = (*svect.getSparsePtr());
    //std::copy((vect.Sparse)->begin(), (vect.Sparse)->end(), (svect.getSparsePtr())->begin());
  }

  // Note FP: using constructor + noalias = (or std::copy) is more efficient than a call to ublas::vector copy constructor,
  // this for large or small vectors.
}

// Copy
SimpleVector::SimpleVector(const SiconosVector &v): SiconosVector(1, v.size()) // Dense by default
{
  unsigned int numV = v.getNum();
  if (numV == 0) // ie if v is a BlockVector, "this" is set to dense.
  {
    // num = 1; default value
    vect.Dense = new DenseVect(sizeV);

    ConstBlockVectIterator it;
    unsigned int pos = 0;
    for (it = v.begin(); it != v.end(); ++it)
    {
      setBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else if (numV == 1) // dense
  {
    // num = 1; default value
    vect.Dense = new DenseVect(v.size());
    noalias(*vect.Dense) = (*v.getDensePtr());
    //std::copy((v.getDensePtr())->begin(), (v.getDensePtr())->end(), (vect.Dense)->begin());
  }
  else // sparse
  {
    num = 4;
    vect.Sparse = new SparseVect(v.size());
    std::copy((v.getSparsePtr())->begin(), (v.getSparsePtr())->end(), (vect.Sparse)->begin());
  }
}

SimpleVector::SimpleVector(const DenseVect& m): SiconosVector(1, m.size())
{
  vect.Dense = new DenseVect(m.size());
  noalias(*vect.Dense) = m;

}

SimpleVector::SimpleVector(const SparseVect& m): SiconosVector(4, m.size())
{
  vect.Sparse = new SparseVect(m.size());
  noalias(*vect.Sparse) = m;
}

SimpleVector::SimpleVector(const std::string &file, bool ascii): SiconosVector(1)
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
  sizeV = (vect.Dense)->size();
}

SimpleVector::~SimpleVector()
{
  if (num == 1)
  {
    delete(vect.Dense);
    vect.Dense = NULL;
  }

  else if (num == 4)
    delete(vect.Sparse);
}

// =================================================
//        get Ublas component (dense or sparse)
// =================================================

const DenseVect SimpleVector::getDense(unsigned int) const
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

SparseVect* SimpleVector::getSparsePtr(unsigned int)const
{

  if (num != 4)
    SiconosVectorException::selfThrow("SimpleVector::getSparsePtr(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return vect.Sparse;
}

double* SimpleVector::getArray(unsigned int) const
{
  if (num == 4)
    SiconosVectorException::selfThrow("SimpleVector::getArray() : not yet implemented for sparse vector.");

  return &(((*vect.Dense).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void SimpleVector::zero()
{
  if (num == 1)
    atlas::set(0.0, *vect.Dense);

  else //if(num==4)
    *vect.Sparse *= 0.0;
}

void SimpleVector::setVector(unsigned int, const SiconosVector& newV)
{
  if (newV.size() != sizeV)
    SiconosVectorException::selfThrow("SimpleVector::setVector(num,v), unconsistent sizes.");

  *this = newV ;
}

void SimpleVector::fill(double value)
{
  if (num == 4)
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
  sizeV = n;
  if (num == 1)
    (vect.Dense)->resize(n, preserve);
  if (num == 4)
    (vect.Sparse)->resize(n, preserve);
}

//=======================
//       get norm
//=======================

const double SimpleVector::normInf() const
{
  if (num == 1)
    return norm_inf(*vect.Dense);
  else //if(num==4)
    return norm_inf(*vect.Sparse);
}

const double SimpleVector::norm2() const
{
  if (num == 1)
    return ublas::norm_2(*vect.Dense);
  else //if(num==4)
    return ublas::norm_2(*vect.Sparse);
}

//=====================
// screen display
//=====================

void SimpleVector::display()const
{
  if (num == 1)
    std::cout << *vect.Dense << std::endl;
  else if (num == 4)
    std::cout << *vect.Sparse << std::endl;
}

//============================
// Convert vector to a string
//============================

const std::string SimpleVector::toString() const
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

//=============================
// Elements access (get or set)
//=============================

const double SimpleVector::getValue(unsigned int row) const
{
  if (row >= sizeV)
    SiconosVectorException::selfThrow("SimpleVector::getValue(index) : Index out of range");

  if (num == 1)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row);
}

void SimpleVector::setValue(unsigned int row, double value)
{
  if (row >= sizeV)
    SiconosVectorException::selfThrow("SimpleVector::setValue(index, value) : Index out of range");
  if (num == 1)
    (*vect.Dense)(row) = value ;
  else
    (*vect.Sparse)(row) = value;
}

double& SimpleVector::operator()(unsigned int row)
{
  if (row >= sizeV)
    SiconosVectorException::selfThrow("SimpleVector::operator ( index ): Index out of range");

  if (num == 1)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

const double SimpleVector::operator()(unsigned int row) const
{
  if (row >= sizeV)
    SiconosVectorException::selfThrow("SimpleVector::operator ( index ): Index out of range");

  if (num == 1)
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

  if (index > sizeV)
    SiconosVectorException::selfThrow("SimpleVector::setBlock : invalid ranges");

  unsigned int end = vIn.size() + index;
  if (end > sizeV)
    SiconosVectorException::selfThrow("SimpleVector::setBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();
  if (numVin == 0) // if vIn is a BlockVector
  {
    ConstBlockVectIterator it;
    unsigned int pos = index;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      setBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else // vIn a SimpleVector ...
  {
    if (numVin != num) SiconosVectorException::selfThrow("SimpleVector::setBlock: inconsistent types.");

    if (num == 1)
      noalias(ublas::subrange(*vect.Dense, index, end)) = *vIn.getDensePtr();
    else
      noalias(ublas::subrange(*vect.Sparse, index, end)) = *vIn.getSparsePtr();
  }
}

void setBlock(const SiconosVector * vIn, SiconosVector * vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut)
{
  // To copy a subBlock of vIn (from position startIn to startIn+sizeB) into vOut (from pos. startOut to startOut+sizeB).
  if (vIn == vOut) // useless op. => nothing to be done
  {}// SiconosVectorException::selfThrow("");
  else
  {
    // Check dim ...
    unsigned int sizeIn = vIn->size();
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

    unsigned int numIn = vIn->getNum();
    unsigned int numOut = vOut->getNum();

    if (numIn == 0 && numOut == 0)
      SiconosVectorException::selfThrow("vector setBlock(v1,v2,...): not yet implemented for v1 and v2 both BlockVectors. Try to use setBlock on the sub-vectors?");
    else if (numOut == 0) // vOut is block ...
    {
      // We look for the block of vOut that include index startOut
      unsigned int blockOutStart = 0;
      const Index * tabOut = vOut->getTabIndexPtr();
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
        setBlock(vIn, vOut->getVectorPtr(blockOutStart), sizeB, startIn, posOut);
      }
      else // More that one block of vOut are concerned
      {

        // The current considered block ...
        SiconosVector * currentBlock = vOut->getVectorPtr(blockOutStart);

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
          currentBlock =  vOut->getVectorPtr(currentBlockNum);
          subSizeB = currentBlock->size();
          setBlock(vIn, currentBlock, subSizeB, posIn, 0);
          currentBlockNum++;
        }
        // set last subBlock ...
        currentBlock =  vOut->getVectorPtr(blockOutEnd);

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
      const Index * tabIn = vIn->getTabIndexPtr();
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
        setBlock(vIn->getVectorPtr(blockInStart), vOut, sizeB, posIn, startOut);
      }
      else // More that one block of vIn are concerned
      {

        // The current considered block ...
        const SiconosVector * currentBlock = vIn->getVectorPtr(blockInStart);

        // Size of the subBlock of vIn to be set.
        unsigned int subSizeB = currentBlock->size() - posIn;
        unsigned int posOut = startOut;

        // Set vOut values, between index posOut and posOut+subSizeB,
        // with first sub-block (currentBlock) of vIn values from posIn to posIn+subSizeB.
        setBlock(currentBlock, vOut, subSizeB, posIn, posOut);

        // Other blocks, except number blockInEnd.
        unsigned int currentBlockNum = blockInStart + 1;
        while (currentBlockNum != blockInEnd)
        {
          posOut += subSizeB;
          currentBlock =  vIn->getVectorPtr(currentBlockNum);
          subSizeB = currentBlock->size();
          setBlock(currentBlock, vOut, subSizeB, 0, posOut);
          currentBlockNum++;
        }
        // set last subBlock ...
        currentBlock =  vIn->getVectorPtr(blockInEnd);

        posOut += subSizeB;

        // Relative position of index endIn in vIn[blockInEnd]
        subSizeB = endIn - (*tabIn)[blockInEnd - 1];

        setBlock(currentBlock, vOut, subSizeB, 0, posOut);
      }
    }
    else // neither vIn nor vOut is a BlockVector
    {

      if (numIn == numOut)
      {
        if (numIn == 1) // vIn / vOut are Dense
          noalias(ublas::subrange(*vOut->getDensePtr(), startOut, endOut)) = ublas::subrange(*vIn->getDensePtr(), startIn, startIn + sizeB);
        else // if(numIn == 4)// vIn / vOut are Sparse
          noalias(ublas::subrange(*vOut->getSparsePtr(), startOut, endOut)) = ublas::subrange(*vIn->getSparsePtr(), startIn, startIn + sizeB);
      }
      else // vIn and vout of different types ...
      {
        if (numIn == 1) // vIn Dense
          noalias(ublas::subrange(*vOut->getSparsePtr(), startOut, endOut)) = ublas::subrange(*vIn->getDensePtr(), startIn, startIn + sizeB);
        else // if(numIn == 4)// vIn Sparse
          noalias(ublas::subrange(*vOut->getDensePtr(), startOut, endOut)) = ublas::subrange(*vIn->getSparsePtr(), startIn, startIn + sizeB);
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
  if ((index + end) > sizeV) SiconosVectorException::selfThrow("SimpleVector::addBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();
  if (numVin == 0) // if vIn is a BlockVector
  {
    ConstBlockVectIterator it;
    unsigned int pos = index;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      addBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else // vIn is a SimpleVector ...
  {
    if (numVin != num) SiconosVectorException::selfThrow("SimpleVector::addBlock : inconsistent types.");

    if (num == 1)
      noalias(ublas::subrange(*vect.Dense, index, index + end)) += *vIn.getDensePtr();
    else
      noalias(ublas::subrange(*vect.Sparse, index, index + end)) += *vIn.getSparsePtr();
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
    ConstBlockVectIterator it;
    unsigned int pos = index;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      subBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else //  vIn a SimpleVector ...
  {
    if (numVin != num) SiconosVectorException::selfThrow("SimpleVector::subBlock : inconsistent types.");

    if (num == 1)
      noalias(ublas::subrange(*vect.Dense, index, index + end)) -= *vIn.getDensePtr();
    else
      noalias(ublas::subrange(*vect.Sparse, index, index + end)) -= *vIn.getSparsePtr();
  }
}

//===============
//  Assignment
//===============

SimpleVector& SimpleVector::operator = (const SiconosVector& vIn)
{
  if (&vIn == this) return *this; // auto-assignment.

  if (sizeV != vIn.size())
    SiconosVectorException::selfThrow("SimpleVector::operator = failed: inconsistent sizes.");

  unsigned int vInNum = vIn.getNum();

  if (vInNum == 0) // if vIn is a BlockVector
  {
    ConstBlockVectIterator it;
    unsigned int pos = 0;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      setBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else // vIn is a SimpleVector
  {
    switch (num)
    {
    case 1:
      switch (vInNum)
      {
      case 1:
        //atlas::copy(*vIn.getDensePtr(),*vect.Dense);
        noalias(*vect.Dense) = *vIn.getDensePtr();
        break;
      case 4:
        noalias(*vect.Dense) = *vIn.getSparsePtr();
        break;
      default:
        SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
        break;
      }
      break;
    case 4:
      if (vInNum == 4)
        noalias(*vect.Sparse) = *vIn.getSparsePtr();
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

  if (sizeV != vIn.size())
    SiconosVectorException::selfThrow("SimpleVector::operator = failed: inconsistent sizes.");

  unsigned int vInNum = vIn.getNum();
  switch (num)
  {
  case 1:
    switch (vInNum)
    {
    case 1:
      //atlas::copy(*vIn.getDensePtr(),*vect.Dense);
      noalias(*vect.Dense) = *vIn.getDensePtr();

      break;
    case 4:
      noalias(*vect.Dense) = *vIn.getSparsePtr();
      break;
    default:
      SiconosVectorException::selfThrow("SimpleVector::operator = : invalid type given");
      break;
    }
    break;
  case 4:
    if (vInNum == 4)
      noalias(*vect.Sparse) = *vIn.getSparsePtr();
    else
      SiconosVectorException::selfThrow("SimpleVector::operator = : can not set sparse = dense.");
    break;
  }
  return *this;
}

SimpleVector& SimpleVector::operator = (const DenseVect& d)
{
  if (num != 1)
    SiconosVectorException::selfThrow("SimpleVector::operator = DenseVect : forbidden: the current vector is not dense.");
  if (d.size() != sizeV)
    SiconosVectorException::selfThrow("SimpleVector::operator = DenseVect : inconsistent size.");

  atlas::copy(d, *vect.Dense);
  return *this;
}

SimpleVector& SimpleVector::operator = (const SparseVect& sp)
{
  if (num != 4)
    SiconosVectorException::selfThrow("SimpleVector::operator = SparseVect : current vector is not sparse.");
  if (sp.size() != sizeV)
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
    switch (num)
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
    if (sizeV != vIn.size())
      SiconosVectorException::selfThrow("SimpleVector::operator += failed: inconsistent sizes.");

    ConstBlockVectIterator it;
    unsigned int pos = 0;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      addBlock(pos, **it);
      pos += (*it)->size();
    }
  }
  else  // vIn Simple
  {
    switch (num)
    {
    case 1:
      switch (vInNum)
      {
      case 1:
        noalias(*vect.Dense) += *vIn.getDensePtr();
        break;
      case 4:
        noalias(*vect.Dense) += *vIn.getSparsePtr();
        break;
      default:
        SiconosVectorException::selfThrow("SimpleVector::operator += : invalid type given");
        break;
      }
      break;
    case 4:
      if (vInNum == 4)
        noalias(*vect.Sparse) += *vIn.getSparsePtr();
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
    if (sizeV != vIn.size())
      SiconosVectorException::selfThrow("SimpleVector::operator -= failed: inconsistent sizes.");

    ConstBlockVectIterator it;
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
    switch (num)
    {
    case 1:
      switch (vInNum)
      {
      case 1:
        noalias(*vect.Dense) -= *vIn.getDensePtr();
        break;
      case 4:
        noalias(*vect.Dense) -= *vIn.getSparsePtr();
        break;
      default:
        SiconosVectorException::selfThrow("SimpleVector::operator -= : invalid type given");
        break;
      }
      break;
    case 4:
      if (vInNum == 4)
        noalias(*vect.Sparse) -= *vIn.getSparsePtr();
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
    DenseVect p = *m.getDensePtr();
    atlas::scal(d, p);
    return p;
  }
  else// if(numM==4)
  {
    return (SparseVect)(*m.getSparsePtr() * d);
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
    DenseVect p = *m.getDensePtr();
    atlas::scal(d, p);
    return p;
  }
  else// if(numM==4)
  {
    return (SparseVect)(*m.getSparsePtr() * d);
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
    DenseVect p = *m.getDensePtr();
    atlas::scal((1.0 / d), p);
    return p;
  }

  else// if(numM==4){
    return (SparseVect)(*m.getSparsePtr() / d);
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
      //    atlas::xpy(*x.getDensePtr(),p);
      //    return p;
      return (DenseVect)(*x.getDensePtr() + *y.getDensePtr());
    }
    else
      return (SparseVect)(*x.getSparsePtr() + *y.getSparsePtr());
  }

  else if (numX != 0 && numY != 0  && numX != numY) // x, y SimpleVector with y and x of different types
  {
    if (numX == 1)
      return (DenseVect)(*x.getDensePtr() + *y.getSparsePtr());
    else
      return (DenseVect)(*x.getSparsePtr() + *y.getDensePtr());
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
          noalias(*z.getDensePtr()) = *x.getDensePtr() + *y.getDensePtr() ;
        }
        else
        {
          if (numZ == 1)
            noalias(*z.getDensePtr()) = *x.getSparsePtr() + *y.getSparsePtr() ;
          else
            noalias(*z.getSparsePtr()) = *x.getSparsePtr() + *y.getSparsePtr() ;
        }
      }
      else if (numX != 0 && numY != 0) // x and y of different types => z must be dense.
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector addition, add(x,y,z) failed - z can not be sparse.");
        if (numX == 1)
          noalias(*z.getDensePtr()) = *x.getDensePtr() + *y.getSparsePtr();
        else
          noalias(*z.getDensePtr()) = *x.getSparsePtr() + *y.getDensePtr() ;
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
      //    atlas::xpy(*x.getDensePtr(),p);
      //    return p;
      return (DenseVect)(*x.getDensePtr() - *y.getDensePtr());
    }
    else
      return (SparseVect)(*x.getSparsePtr() - *y.getSparsePtr());
  }

  else if (numX != 0 && numY != 0  && numX != numY) // x, y SimpleVector with y and x of different types
  {
    if (numX == 1)
      return (DenseVect)(*x.getDensePtr() - *y.getSparsePtr());
    else
      return (DenseVect)(*x.getSparsePtr() - *y.getDensePtr());
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
        *z.getDensePtr() = *x.getDensePtr() - *y.getDensePtr() ;
      }
      else
      {
        if (numZ == 1)
          *z.getDensePtr() = *x.getSparsePtr() - *y.getDensePtr() ;
        else
          *z.getSparsePtr() = *x.getSparsePtr() - *y.getSparsePtr() ;
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
          noalias(*z.getDensePtr()) = *x.getDensePtr() - *y.getDensePtr() ;
        }
        else
        {
          if (numZ == 1)
            noalias(*z.getDensePtr()) = *x.getSparsePtr() - *y.getSparsePtr() ;
          else
            noalias(*z.getSparsePtr()) = *x.getSparsePtr() - *y.getSparsePtr() ;
        }
      }
      else if (numX != 0 && numY != 0) // x and y of different types => z must be dense.
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector addition, sub(x,y,z) failed - z can not be sparse.");
        if (numX == 1)
          noalias(*z.getDensePtr()) = *x.getDensePtr() - *y.getSparsePtr();
        else
          noalias(*z.getDensePtr()) = *x.getSparsePtr() - *y.getDensePtr() ;
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
      atlas::axpby(a, *x.getDensePtr(), b, *y.getDensePtr());

    else if (numX == 0) // ie if x and y are block
    {
      if (isComparableTo(&x, &y))
      {
        BlockVectIterator itY;
        ConstBlockVectIterator itX = x.begin();
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
      *y.getSparsePtr() *= b;
      if (&y != &x)
        noalias(*y.getSparsePtr()) += a**x.getSparsePtr();
      else
        *y.getSparsePtr() += a**x.getSparsePtr();
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
        atlas::axpy(a, *tmp.getDensePtr(), *y.getDensePtr());
      else
        SiconosVectorException::selfThrow("axpby failed: try to add block to sparse vector.");
    }
    else // x and y simple but of different types
    {
      if (numX == 1)
        *y.getSparsePtr() += a**x.getDensePtr();
      else
        *y.getDensePtr() +=  a**x.getSparsePtr();
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
      atlas::axpy(a, *x.getDensePtr(), *y.getDensePtr());

    else if (numX == 0) // ie if x and y are block
    {
      if (isComparableTo(&x, &y))
      {
        BlockVectIterator itY;
        ConstBlockVectIterator itX = x.begin();
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
        noalias(*y.getSparsePtr()) += a**x.getSparsePtr();
      else
        *y.getSparsePtr() += a**x.getSparsePtr();
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
        atlas::axpy(a, *tmp.getDensePtr(), *y.getDensePtr());
      else
        SiconosVectorException::selfThrow("axpby failed: try to add block to sparse vector.");
    }
    else // x and y simple but of different types
    {
      if (numX == 1)
        *y.getSparsePtr() += a**x.getDensePtr();
      else
        *y.getDensePtr() +=  a**x.getSparsePtr();
    }
  }
}

const double inner_prod(const SiconosVector &x, const SiconosVector &m)
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
      return atlas::dot(*x.getDensePtr(), *m.getDensePtr());
    else
      return inner_prod(*x.getSparsePtr(), *m.getSparsePtr());
  }
  else if (numM == 1)
    return inner_prod(*x.getSparsePtr(), *m.getDensePtr());
  else
    return inner_prod(*x.getDensePtr(), *m.getSparsePtr());
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
      return (DenseMat)(outer_prod(*x.getDensePtr(), *m.getDensePtr()));

    else// if(numX == 4)
      return (DenseMat)(outer_prod(*x.getSparsePtr(), *m.getDensePtr()));
  }
  else // if(numM == 4)
  {
    if (numX == 1)
      return (DenseMat)(outer_prod(*x.getDensePtr(), *m.getSparsePtr()));

    else //if(numX == 4)
      return (DenseMat)(outer_prod(*x.getSparsePtr(), *m.getSparsePtr()));
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
        if (isComparableTo(&x, &y)) // if x and y are "block-consistent"
        {
          ConstBlockVectIterator itX = x.begin();
          BlockVectIterator itY ;
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
          //atlas::axpby(a,*x.getDensePtr(),0.0,*y.getDensePtr());
          noalias(*y.getDensePtr()) = a * *x.getDensePtr();
        else
          noalias(*y.getDensePtr()) += a * *x.getDensePtr();
      }
      else  // if both are sparse
      {
        if (init)
          noalias(*y.getSparsePtr()) = a**x.getSparsePtr();
        else
          noalias(*y.getSparsePtr()) += a**x.getSparsePtr();
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
            noalias(*y.getDensePtr()) = a**x.getSparsePtr();
          else
            noalias(*y.getDensePtr()) += a**x.getSparsePtr();

        }
        else
          SiconosVectorException::selfThrow("SiconosVector::scal(a,dense,sparse) not allowed.");
      }
    }
  }
}

void subscal(double a, const SiconosVector & x, SiconosVector & y, const std::vector<unsigned int>& coord, bool init)
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

      ConstBlockVectIterator it;
      // Number of the subvector of x that handles element at position coord[0]
      unsigned int firstBlockNum = x.getNumVectorAtPos(coord[0]);
      // Number of the subvector of x that handles element at position coord[1]
      unsigned int lastBlockNum = x.getNumVectorAtPos(coord[1]);
      std::vector<unsigned int> subCoord = coord;
      SiconosVector * tmp = y[firstBlockNum];
      unsigned int subSize =  x[firstBlockNum]->size(); // Size of the sub-vector
      const Index * xTab = x.getTabIndexPtr();
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
            SiconosMatrixException::selfThrow("subscal(a,x,y) error: not yet implemented for x block of blocks ...");
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
      ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2])
      {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<DenseVect> subX(*x.getDensePtr(), ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
    else //if (numX == 4) // Sparse
    {
      ublas::vector_range<SparseVect> subY(*y.getSparsePtr(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2])
      {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<SparseVect> subX(*x.getSparsePtr(), ublas::range(coord[0], coord[1]));
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
        SiconosMatrixException::selfThrow("subscal(a,x,y) error: not yet implemented for x and y block vectors");
      }
      else if (numX == 1) // ie if both are Dense
      {
        ublas::vector_range<DenseVect> subX(*x.getDensePtr(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[2], coord[3]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
      else  // if both are sparse
      {
        ublas::vector_range<SparseVect> subX(*x.getSparsePtr(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<SparseVect> subY(*y.getSparsePtr(), ublas::range(coord[2], coord[3]));

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
        ConstBlockVectIterator it;
        // Number of the subvector of x that handles element at position coord[0]
        unsigned int firstBlockNum = x.getNumVectorAtPos(coord[0]);
        // Number of the subvector of x that handles element at position coord[1]
        unsigned int lastBlockNum = x.getNumVectorAtPos(coord[1]);
        std::vector<unsigned int> subCoord = coord;
        const SiconosVector * tmp = x[firstBlockNum];
        unsigned int subSize =  x[firstBlockNum]->size(); // Size of the sub-vector
        const Index * xTab = x.getTabIndexPtr();
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
          unsigned int xPos = 0 ; // Position in x of the current sub-vector of x
          bool firstLoop = true;
          for (it = x.begin(); it != x.end(); ++it)
          {
            if ((*it)->getNum() == 0)
              SiconosMatrixException::selfThrow("subscal(a,x,y) error: not yet implemented for x block of blocks ...");
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
        ConstBlockVectIterator it;
        // Number of the subvector of y that handles element at position coord[2]
        unsigned int firstBlockNum = y.getNumVectorAtPos(coord[2]);
        // Number of the subvector of x that handles element at position coord[3]
        unsigned int lastBlockNum = y.getNumVectorAtPos(coord[3]);
        std::vector<unsigned int> subCoord = coord;
        SiconosVector * tmp = y[firstBlockNum];
        unsigned int subSize =  y[firstBlockNum]->size(); // Size of the sub-vector
        const Index * yTab = y.getTabIndexPtr();
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
          unsigned int yPos = 0 ; // Position in x of the current sub-vector of x
          bool firstLoop = true;
          for (it = y.begin(); it != y.end(); ++it)
          {
            if ((*it)->getNum() == 0)
              SiconosMatrixException::selfThrow("subscal(a,x,y) error: not yet implemented for y block of blocks ...");
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
        ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[2], coord[3]));
        ublas::vector_range<SparseVect> subX(*x.getSparsePtr(), ublas::range(coord[0], coord[1]));

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
