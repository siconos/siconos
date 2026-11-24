// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SiconosMemory.hpp"

#include <assert.h>

#include <iostream>

#include "SiconosVector.hpp"

// From the size of the container (number of saved vectors) and
// the size of the vectors.
siconos::algebra::SiconosMemory::SiconosMemory(
    const siconos::algebra::blocks::size_type size,
    const siconos::algebra::blocks::size_type vectorSize)
    : siconos::algebra::blocks::Vector(), _indx(size - 1) {
  for (siconos::algebra::blocks::size_type i = 0; i < size; i++) {
    push_back(siconos::algebra::SiconosVector(vectorSize));
  }
}

// Copy constructor
siconos::algebra::SiconosMemory::SiconosMemory(const SiconosMemory& Mem)
    : siconos::algebra::blocks::Vector(),
      _nbVectorsInMemory(Mem.nbVectorsInMemory()),
      _indx(Mem.size() - 1) {
  for (siconos::algebra::blocks::size_type i = 0; i < Mem.size(); i++) {
    push_back(Mem[i]);
  }
}

siconos::algebra::SiconosMemory& siconos::algebra::SiconosMemory::operator=(
    const SiconosMemory& V) {
  if (size() != V.size()) {
    this->resize(V.size());  // => copy construction of old content
  }

  for (siconos::algebra::blocks::size_type i = 0; i < V.size(); i++) {
    (*this)[i].resize(V[i].size(), true);
    (*this)[i] = V[i];  // copy
  }
  _indx = V._indx;
  _nbVectorsInMemory = V._nbVectorsInMemory;
  return *this;
}

// (Re)set the size of an existing SiconosMemory
void siconos::algebra::SiconosMemory::setMemorySize(
    const siconos::algebra::blocks::size_type steps,
    const siconos::algebra::Index vectorSize) {
  _nbVectorsInMemory = 0;
  _indx = steps - 1;
  auto sizev = size();
  for (decltype(sizev) i = 0; i < sizev; i++) {
    this->at(i).resize(vectorSize, true);
  }
  for (decltype(sizev) i = size(); i < steps; i++) {
    this->push_back(siconos::algebra::SiconosVector(vectorSize));
  }
}

// --- GETTERS/SETTERS ---

const siconos::algebra::SiconosVector& siconos::algebra::SiconosMemory::getSiconosVector(
    const siconos::algebra::blocks::size_type index) const {
  assert(index < _nbVectorsInMemory && "getSiconosVector(index) : inconsistent index value");
  return this->at((_indx + 1 + index) % this->size());
}

siconos::algebra::SiconosVector& siconos::algebra::SiconosMemory::getSiconosVectorMutable(
    const siconos::algebra::blocks::size_type index) {
  assert(index < _nbVectorsInMemory && "getSiconosVector(index) : inconsistent index value");
  return *(siconos::algebra::SiconosVector*)(&this->at((_indx + 1 + index) % this->size()));
}

void siconos::algebra::SiconosMemory::swap(const siconos::algebra::SiconosVector& v) {
  // Be robust to empty memory
  if (size() == 0) return;

  // If _nbVectorsInMemory is this->size(), we remove the last element.
  (*this)[_indx] = v;
  _nbVectorsInMemory = std::min(_nbVectorsInMemory + 1, this->size());
  if (_indx > 0)
    _indx--;
  else
    _indx = this->size() - 1;
}

void siconos::algebra::SiconosMemory::swap(
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  // Be robust to empty memory
  // Be robust to null pointer
  if (size() == 0 || !v) return;

  // If _nbVectorsInMemory is this->size(), we remove the last element.
  (*this)[_indx] = *v;
  _nbVectorsInMemory = std::min(_nbVectorsInMemory + 1, this->size());
  if (_indx > 0)
    _indx--;
  else
    _indx = this->size() - 1;
}

void siconos::algebra::print(const SiconosMemory& vec) {
  std::cout << " ====== Memory vector display ======= \n";
  std::cout << "| size : " << vec.size() << "\n";
  std::cout << "| _nbVectorsInMemory : " << vec.nbVectorsInMemory() << "\n";
  for (siconos::algebra::blocks::size_type i = 0; i < vec.nbVectorsInMemory(); i++) {
    std::cout << "vector number " << i << ": address = " << vec.at(i) << " | \n ";
  }
  std::cout << " ===================================== \n";
}
