/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "SiconosConfig.h"

#include "BlockCSRMatrix.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "NonSmoothLaw.hpp"

#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "SimulationGraphs.hpp"
#include "SparseBlockMatrix.h" // From numerics, for SparseBlockStructuredMatrix
#include "Tools.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include "debug.h"

// Default constructor: empty matrix
BlockCSRMatrix::BlockCSRMatrix():
  _nr(0),
  _blockCSR(new CompressedRowMat()),
  _sparseBlockStructuredMatrix(new SparseBlockStructuredMatrix()),
  _diagsize0(new IndexInt()),
  _diagsize1(new IndexInt()),
  rowPos(new IndexInt()),
  colPos(new IndexInt())
{}

// Constructor with dimensions
BlockCSRMatrix::BlockCSRMatrix(unsigned int nRow):
  _nr(nRow),
  // Only square-blocks matrices for the moment (ie nRow = nr = nrol)

  // Allocate memory and fill in the matrix rowPos, rowCol ... are
  // initialized with nr to reserve at first step the maximum possible
  // (according to given nr) space in memory.  Thus a future resize
  // will not require memory allocation or copy.
  _blockCSR(new CompressedRowMat(_nr, _nr)),
  _sparseBlockStructuredMatrix(new SparseBlockStructuredMatrix()),
  _diagsize0(new IndexInt(_nr)),
  _diagsize1(new IndexInt(_nr)),
  rowPos(new IndexInt(_nr)),
  colPos(new IndexInt(_nr))
{}

// Basic constructor
BlockCSRMatrix::BlockCSRMatrix(InteractionsGraph& indexSet):
  _nr(indexSet.size()),
  _blockCSR(new CompressedRowMat(_nr, _nr)),
  _sparseBlockStructuredMatrix(new SparseBlockStructuredMatrix()),
  _diagsize0(new IndexInt(_nr)),
  _diagsize1(new IndexInt(_nr)),
  rowPos(new IndexInt(_nr)),
  colPos(new IndexInt(_nr))
{
  DEBUG_BEGIN("BlockCSRMatrix::BlockCSRMatrix(SP::InteractionsGraph indexSet)\n");
  fill(indexSet);
  DEBUG_END("BlockCSRMatrix::BlockCSRMatrix(SP::InteractionsGraph indexSet)\n");
}

BlockCSRMatrix::~BlockCSRMatrix()
{}

// Fill the SparseMat
void BlockCSRMatrix::fill(InteractionsGraph& indexSet)
{
  // ======> Aim: find inter1 and inter2 both in indexSets[level] and which
  // have common DynamicalSystems.  Then get the corresponding matrix
  // from map blocks.

  // Number of blocks in a row = number of active constraints.
  _nr = indexSet.size();

  // (re)allocate memory for ublas matrix
  _blockCSR->resize(_nr, _nr, false);

  _diagsize0->resize(_nr);
  _diagsize1->resize(_nr);

  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===


  int sizeV = 0;

  InteractionsGraph::VIterator vi, viend;
  for (std11::tie(vi, viend) = indexSet.vertices();
       vi != viend; ++vi)
  {
    SP::Interaction inter = indexSet.bundle(*vi);

    assert(inter->nonSmoothLaw()->size() > 0);

    sizeV  += inter->nonSmoothLaw()->size();
    (*_diagsize0)[indexSet.index(*vi)] = sizeV;
    (*_diagsize1)[indexSet.index(*vi)] = sizeV;
    assert((*_diagsize0)[indexSet.index(*vi)] > 0);
    assert((*_diagsize1)[indexSet.index(*vi)] > 0);

    (*_blockCSR)(indexSet.index(*vi), indexSet.index(*vi)) =
      indexSet.properties(*vi).block->getArray();
  }

  InteractionsGraph::EIterator ei, eiend;
  for (std11::tie(ei, eiend) = indexSet.edges();
       ei != eiend; ++ei)
  {
    InteractionsGraph::VDescriptor vd1 = indexSet.source(*ei);
    InteractionsGraph::VDescriptor vd2 = indexSet.target(*ei);
    SP::Interaction inter1 = indexSet.bundle(vd1);
    SP::Interaction inter2 = indexSet.bundle(vd2);

    assert(indexSet.index(vd1) < _nr);
    assert(indexSet.index(vd2) < _nr);

    assert(indexSet.is_vertex(inter2));

    assert(vd2 == indexSet.descriptor(inter2));
    assert(indexSet.index(vd2) == indexSet.index(indexSet.descriptor(inter2)));


    unsigned int pos = indexSet.index(vd1);
    unsigned int col = indexSet.index(vd2);

    assert(pos != col);

    (*_blockCSR)(std::min(pos, col), std::max(pos, col)) =
      indexSet.properties(*ei).upper_block->getArray();

    (*_blockCSR)(std::max(pos, col), std::min(pos, col)) =
      indexSet.properties(*ei).lower_block->getArray();
  }
  DEBUG_EXPR(display(););
}

void BlockCSRMatrix::fillM(InteractionsGraph& indexSet)
{
  /* on adjoint graph a dynamical system may be on several edges */
  std::map<SP::DynamicalSystem, bool> involvedDS;
  InteractionsGraph::EIterator ei, eiend;
  for(std11::tie(ei, eiend) = indexSet.edges();
      ei != eiend; ++ei)
  {
    if (Type::value(*indexSet.bundle(*ei)) != Type::NewtonEulerDS)
    {
      RuntimeException::selfThrow("BlockCSRMatrix::fillM only for Newton EulerDS");
    }

    _nr = 0;
    
    if (involvedDS.find(indexSet.bundle(*ei)) == involvedDS.end())
    {
      _nr++;
      involvedDS[indexSet.bundle(*ei)] = true;
      _blockCSR->resize(_nr, _nr, false);

      (*_blockCSR)(_nr-1, _nr-1) = std11::static_pointer_cast<NewtonEulerDS>
        (indexSet.bundle(*ei))->mass()->getArray();
    }
  }
  
  _diagsize0->resize(involvedDS.size());
  _diagsize1->resize(involvedDS.size());

  /* here we suppose NewtonEuler with 6 dofs */
  /* it cannot be another case at this point */
  unsigned int index, ac;
  for (index = 0, ac = 6; 
       index < involvedDS.size();
       ++index, ac+=6)
  {
    (*_diagsize0)[index] = ac;
    (*_diagsize1)[index] = ac;
  }
  
}

void BlockCSRMatrix::fillH(InteractionsGraph& indexSet)
{
  /* on adjoint graph a dynamical system may be on several edges */
  std::map<SP::DynamicalSystem, unsigned int> involvedDS;
  InteractionsGraph::EIterator ei, eiend;
  {
    unsigned int index;
    for(std11::tie(ei, eiend) = indexSet.edges(), index=0;
        ei != eiend; ++ei, ++index)
    {
      if (involvedDS.find(indexSet.bundle(*ei)) == involvedDS.end())
      {
        if (Type::value(*indexSet.bundle(*ei)) != Type::NewtonEulerDS)
        {
          RuntimeException::selfThrow("BlockCSRMatrix::fillH only for Newton EulerDS");
        }
        involvedDS[indexSet.bundle(*ei)] = index;
      }
    }
  }

  _nr = involvedDS.size();

  _blockCSR->resize(_nr, _nr, false);

  InteractionsGraph::VIterator vi, viend;
  for(std11::tie(vi, viend) = indexSet.vertices();
      vi != viend; ++vi)
  {

    SP::DynamicalSystem first = SP::DynamicalSystem();
    unsigned int pos=0, col=0;
    InteractionsGraph::EDescriptor ed1, ed2;
    InteractionsGraph::OEIterator oei, oeiend;
    for(std11::tie(oei, oeiend) = indexSet.out_edges(*vi);
        oei != oeiend; ++oei)
    {
      if (!first)
      {
        first = indexSet.bundle(*oei);
        col = involvedDS[first];
        pos = involvedDS[first];
      }
      else
      {
        if (indexSet.bundle(*oei) != first)
        {
          pos = involvedDS[indexSet.bundle(*oei)];
        }
      }
    }

    (*_blockCSR)(std::min(pos, col), std::max(pos, col)) = 
      std11::static_pointer_cast<NewtonEulerR>(indexSet.bundle(*vi)->relation())->jachqT()->getArray();
    
    (*_blockCSR)(std::max(pos, col), std::min(pos, col)) = 
      std11::static_pointer_cast<NewtonEulerR>(indexSet.bundle(*vi)->relation())->jachqT()->getArray();
    
  }
  
  _diagsize0->resize(involvedDS.size());
  _diagsize1->resize(involvedDS.size());
  
  /* only NewtonEulerFrom3DLocalFrameR */
  unsigned int index, ac0, ac1;
  for (index= 0, ac0 = 6, ac1 = 3; 
       index < involvedDS.size();
       ++index, ac0 +=6, ac1 +=3)
  {
    (*_diagsize0)[index] = ac0;
    (*_diagsize1)[index] = ac1;
  }
  
}


// convert _blockCSR to numerics structure
void BlockCSRMatrix::convert()
{
  _sparseBlockStructuredMatrix->blocknumber0 = _nr;
  _sparseBlockStructuredMatrix->blocknumber1 = _nr;  // nc not always set
  _sparseBlockStructuredMatrix->nbblocks = (*_blockCSR).nnz();
  // Next copies: pointer links!!
  _sparseBlockStructuredMatrix->blocksize0 =  _diagsize0->data();
  _sparseBlockStructuredMatrix->blocksize1 =  _diagsize1->data(); // nr = nc

  // boost
  _sparseBlockStructuredMatrix->filled1 = (*_blockCSR).filled1();
  _sparseBlockStructuredMatrix->filled2 = (*_blockCSR).filled2();
  _sparseBlockStructuredMatrix->index1_data = _blockCSR->index1_data().begin();
  if (_nr > 0)
  {
    _sparseBlockStructuredMatrix->index2_data = _blockCSR->index2_data().begin();
    _sparseBlockStructuredMatrix->block =  _blockCSR->value_data().begin();
  };

  //   // Loop through the non-null blocks
  //   for (SpMatIt1 i1 = _blockCSR->begin1(); i1 != _blockCSR->end1(); ++i1)
  //     {
  //       for (SpMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2)
  //  {
  //    block[i] = *i2;
  //  }
  //     }
}

// Display data
void BlockCSRMatrix::display() const
{
  std::cout << "----- Sparse Block Matrix with "
            << _nr << " blocks in a row/col and "
            << _blockCSR->nnz()
            << " non-null blocks" <<std::endl;
  std::cout << "filled1 (index of the last non empty line + 1):" << _blockCSR->filled1() <<std::endl;
  std::cout << "filled2 (number of non null blocks):" << _blockCSR->filled2() <<std::endl;
  std::cout << "_blockCSR->index1_data().size()" << _blockCSR->index1_data().size() << std::endl;
  print(_blockCSR->index1_data().begin(), _blockCSR->index1_data().end(),"index1_data", "\t");

  assert(_blockCSR->index2_data().size() >= _blockCSR->filled2() );

  std::cout << "_blockCSR->index2_data().size()" << _blockCSR->index2_data().size() << std::endl;
  print(_blockCSR->index2_data().begin(), _blockCSR->index2_data().end(), "index2_data (column number for each block)", "\t");

  std::cout << "last column number  "<<   _blockCSR->index2_data()[_blockCSR->filled2()-1] <<  " for block   " << _blockCSR->filled2() << std::endl;
  print(_diagsize0->begin(), _diagsize0->end(),"_diagsize0 , sum of row sizes of the diagonal blocks", "\t" );
  print(_diagsize1->begin(), _diagsize1->end(),"_diagsize1 , sum of col sizes of the diagonal blocks", "\t"  );
}

unsigned int BlockCSRMatrix::getNbNonNullBlocks() const
{
  return _blockCSR->nnz();
};
