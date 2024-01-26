/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "BlockCSRMatrix.hpp"

#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "Interaction.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NonSmoothLaw.hpp"
#include "NumericsToolsNamespace.h"  // for SparseBlockStructuredmatrix
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "Tools.hpp"     // For print
#include "TypeName.hpp"  // for DS type visitor

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

// Only square-blocks matrices for the moment (ie nRow = nr = nrol)

// Allocate memory and fill in the matrix rowPos, rowCol ... are
// initialized with nr to reserve at first step the maximum possible
// (according to given nr) space in memory.  Thus a future resize
// will not require memory allocation or copy.
siconos::simulation::BlockCSRMatrix::BlockCSRMatrix(unsigned int nRow) : _nr(nRow), _nc{_nr} {
  _blockCSR = std::make_shared<CompressedRowMat>(_nr, _nr);
  _sparseBlockStructuredMatrix = std::make_shared<SparseBlockStructuredMatrix>();
  _diagsize0 = std::make_shared<std::vector<unsigned int>>(_nr);
  _diagsize1 = std::make_shared<std::vector<unsigned int>>(_nr);
  rowPos = std::make_shared<std::vector<unsigned int>>(_nr);
  colPos = std::make_shared<std::vector<unsigned int>>(_nr);
}

// Basic constructor
siconos::simulation::BlockCSRMatrix::BlockCSRMatrix(
    siconos::graphs::InteractionsGraph& indexSet)
    : BlockCSRMatrix(indexSet.size()) {
  DEBUG_BEGIN("siconos::simulation::BlockCSRMatrix::BlockCSRMatrix(auto indexSet)\n");
  fill(indexSet);
  DEBUG_END(
      "siconos::simulation::BlockCSRMatrix::BlockCSRMatrix(std::shared_ptr<siconos::graphs::"
      "siconos::graphs::"
      "InteractionsGraph> indexSet)\n");
}

// Fill the SparseMat
void siconos::simulation::BlockCSRMatrix::fill(siconos::graphs::InteractionsGraph& indexSet) {
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

  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  for (std::tie(vi, viend) = indexSet.vertices(); vi != viend; ++vi) {
    std::shared_ptr<siconos::modeling::Interaction> inter = indexSet.bundle(*vi);

    assert(inter->nonSmoothLaw()->size() > 0);

    sizeV += inter->nonSmoothLaw()->size();
    (*_diagsize0)[indexSet.index(*vi)] = sizeV;
    (*_diagsize1)[indexSet.index(*vi)] = sizeV;
    assert((*_diagsize0)[indexSet.index(*vi)] > 0);
    assert((*_diagsize1)[indexSet.index(*vi)] > 0);

    (*_blockCSR)(indexSet.index(*vi), indexSet.index(*vi)) =
        indexSet.properties(*vi).block->getArray();
  }

  siconos::graphs::InteractionsGraph::EIterator ei, eiend;
  for (std::tie(ei, eiend) = indexSet.edges(); ei != eiend; ++ei) {
    auto vd1 = indexSet.source(*ei);
    auto vd2 = indexSet.target(*ei);
    auto inter1 = indexSet.bundle(vd1);
    auto inter2 = indexSet.bundle(vd2);

    assert(indexSet.index(vd1) < _nr);
    assert(indexSet.index(vd2) < _nr);

    assert(indexSet.is_vertex(inter2));

    assert(vd2 == indexSet.descriptor(inter2));
    assert(indexSet.index(vd2) == indexSet.index(indexSet.descriptor(inter2)));

    auto pos = indexSet.index(vd1);
    auto col = indexSet.index(vd2);

    assert(pos != col);

    (*_blockCSR)(std::min(pos, col), std::max(pos, col)) =
        indexSet.properties(*ei).upper_block->getArray();

    (*_blockCSR)(std::max(pos, col), std::min(pos, col)) =
        indexSet.properties(*ei).lower_block->getArray();
  }
  DEBUG_EXPR(display(););
}

void siconos::simulation::BlockCSRMatrix::fillW(siconos::graphs::InteractionsGraph& indexSet) {
  /* on adjoint graph a dynamical system may be on several edges */
  std::map<std::shared_ptr<siconos::modeling::DynamicalSystem>, bool> involvedDS;
  siconos::graphs::InteractionsGraph::EIterator ei, eiend;
  for (std::tie(ei, eiend) = indexSet.edges(); ei != eiend; ++ei) {
    if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(
            indexSet.bundle(*ei))) {
      _nr = 0;

      if (involvedDS.find(indexSet.bundle(*ei)) == involvedDS.end()) {
        _nr++;
        involvedDS[indexSet.bundle(*ei)] = true;
        _blockCSR->resize(_nr, _nr, false);

        (*_blockCSR)(_nr - 1, _nr - 1) = neds->mass()->getArray();
      }
    } else {
      THROW_EXCEPTION("siconos::simulation::BlockCSRMatrix::fillW only for Newton EulerDS");
    }
  }

  _diagsize0->resize(involvedDS.size());
  _diagsize1->resize(involvedDS.size());

  /* here we suppose NewtonEuler with 6 dofs */
  /* it cannot be another case at this point */
  unsigned int index, ac;
  for (index = 0, ac = 6; index < involvedDS.size(); ++index, ac += 6) {
    (*_diagsize0)[index] = ac;
    (*_diagsize1)[index] = ac;
  }
}

void siconos::simulation::BlockCSRMatrix::fillH(siconos::graphs::InteractionsGraph& indexSet) {
  /* on adjoint graph a dynamical system may be on several edges */
  std::map<std::shared_ptr<siconos::modeling::DynamicalSystem>, unsigned int> involvedDS;
  siconos::graphs::InteractionsGraph::EIterator ei, eiend;
  {
    unsigned int index;
    for (std::tie(ei, eiend) = indexSet.edges(), index = 0; ei != eiend; ++ei, ++index) {
      if (involvedDS.find(indexSet.bundle(*ei)) == involvedDS.end()) {
        if (siconos::types::type_value(*indexSet.bundle(*ei)) !=
            siconos::modeling::Type::NewtonEulerDS) {
          THROW_EXCEPTION(
              "siconos::simulation::BlockCSRMatrix::fillH only for Newton EulerDS");
        }
        involvedDS[indexSet.bundle(*ei)] = index;
      }
    }
  }

  _nr = involvedDS.size();

  _blockCSR->resize(_nr, _nr, false);

  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  for (std::tie(vi, viend) = indexSet.vertices(); vi != viend; ++vi) {
    std::shared_ptr<siconos::modeling::DynamicalSystem> first{nullptr};
    unsigned int pos = 0, col = 0;
    siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
    siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = indexSet.out_edges(*vi); oei != oeiend; ++oei) {
      if (!first) {
        first = indexSet.bundle(*oei);
        col = involvedDS[first];
        pos = involvedDS[first];
      } else {
        if (indexSet.bundle(*oei) != first) {
          pos = involvedDS[indexSet.bundle(*oei)];
        }
      }
    }

    (*_blockCSR)(std::min(pos, col), std::max(pos, col)) =
        std::static_pointer_cast<siconos::modeling::NewtonEulerR>(
            indexSet.bundle(*vi)->relation())
            ->jachqT()
            ->getArray();

    (*_blockCSR)(std::max(pos, col), std::min(pos, col)) =
        std::static_pointer_cast<siconos::modeling::NewtonEulerR>(
            indexSet.bundle(*vi)->relation())
            ->jachqT()
            ->getArray();
  }

  _diagsize0->resize(involvedDS.size());
  _diagsize1->resize(involvedDS.size());

  /* only NewtonEuler3DR */
  unsigned int index, ac0, ac1;
  for (index = 0, ac0 = 6, ac1 = 3; index < involvedDS.size(); ++index, ac0 += 6, ac1 += 3) {
    (*_diagsize0)[index] = ac0;
    (*_diagsize1)[index] = ac1;
  }
}

// convert _blockCSR to numerics structure
void siconos::simulation::BlockCSRMatrix::convert() {
  DEBUG_BEGIN("void siconos::simulation::BlockCSRMatrix::convert()\n");
  _sparseBlockStructuredMatrix->blocknumber0 = _nr;
  _sparseBlockStructuredMatrix->blocknumber1 = _nr;  // nc not always set
  _sparseBlockStructuredMatrix->nbblocks = (*_blockCSR).nnz();
  // Next copies: pointer links!!
  _sparseBlockStructuredMatrix->blocksize0 = _diagsize0->data();
  _sparseBlockStructuredMatrix->blocksize1 = _diagsize1->data();  // nr = nc

  // boost
  _sparseBlockStructuredMatrix->filled1 = (*_blockCSR).filled1();
  _sparseBlockStructuredMatrix->filled2 = (*_blockCSR).filled2();
  _sparseBlockStructuredMatrix->index1_data = _blockCSR->index1_data().begin();
  if (_nr > 0) {
    _sparseBlockStructuredMatrix->index2_data = _blockCSR->index2_data().begin();
    _sparseBlockStructuredMatrix->block = _blockCSR->value_data().begin();
  };
  if (_sparseBlockStructuredMatrix->diagonal_blocks) {
    free(_sparseBlockStructuredMatrix->diagonal_blocks);
    _sparseBlockStructuredMatrix->diagonal_blocks = nullptr;
  }
  //   // Loop through the non-null blocks
  //   for (SparseMat::iterator1 i1 = _blockCSR->begin1(); i1 != _blockCSR->end1(); ++i1)
  //     {
  //       for (SparseMat::iterator2 i2 = i1.begin(); i2 != i1.end(); ++i2)
  //  {
  //    block[i] = *i2;
  //  }
  //     }
  DEBUG_END("void siconos::simulation::BlockCSRMatrix::convert()\n");
}

// Display data
void siconos::simulation::BlockCSRMatrix::display() const {
  std::cout << "----- Sparse Block Matrix with " << _nr << " blocks in a row/col and "
            << _blockCSR->nnz() << " non-null blocks\n";
  ;
  std::cout << "filled1 (index of the last non empty line + 1):" << _blockCSR->filled1()
            << "\n";
  std::cout << "filled2 (number of non null blocks):" << _blockCSR->filled2() << "\n";
  std::cout << "_blockCSR->index1_data().size()" << _blockCSR->index1_data().size() << "\n";
  // siconos::tools::print(_blockCSR->index1_data().begin(), _blockCSR->index1_data().end(),
  //                       "index1_data", "\t");
  siconos::tools::print("index1_data\t", _blockCSR->index1_data());

  assert(_blockCSR->index2_data().size() >= _blockCSR->filled2());

  std::cout << "_blockCSR->index2_data().size()" << _blockCSR->index2_data().size() << "\n";
  // siconos::tools::print(_blockCSR->index2_data().begin(), _blockCSR->index2_data().end(),
  //                       "index2_data (column number for each block)", "\t");
  siconos::tools::print("index2_data (column number for each block)\t",
                        _blockCSR->index2_data());

  std::cout << "last column number  " << _blockCSR->index2_data()[_blockCSR->filled2() - 1]
            << " for block   " << _blockCSR->filled2() << "\n";
  // siconos::tools::print(_diagsize0->begin(), _diagsize0->end(),
  //                       "_diagsize0 , sum of row sizes of the diagonal blocks", "\t");
  // siconos::tools::print(_diagsize1->begin(), _diagsize1->end(),
  //                       "_diagsize1 , sum of col sizes of the diagonal blocks", "\t");
  siconos::tools::print("_diagsize0 , sum of row sizes of the diagonal blocks\t", *_diagsize0);
  siconos::tools::print("_diagsize1 , sum of col sizes of the diagonal blocks\t", *_diagsize1);
}

unsigned int siconos::simulation::BlockCSRMatrix::getNbNonNullBlocks() const {
  return _blockCSR->nnz();
};
