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
#include "SiconosMemory.hpp"
#include "BlockVector.hpp"
#include "SiconosVector.hpp"

#include <iostream>


// --- CONSTRUCTORS ---


// from data: _size
SiconosMemory::SiconosMemory(const unsigned int size, const unsigned int vectorSize):
  _size(size), _nbVectorsInMemory(0)
{
  _vectorMemory.reset(new MemoryContainer);
  _vectorMemory->resize(_size);
  _indx = _size-1;
  for (unsigned int i = 0; i < _size; i++)
  {
    (*_vectorMemory)[i].reset(new SiconosVector(vectorSize));
  }
}

// copy of a std::vector of siconos vectors
SiconosMemory::SiconosMemory(const MemoryContainer& V):
  _size(V.size()), _nbVectorsInMemory(V.size())
{
  _indx = _size-1;
  _vectorMemory.reset(new MemoryContainer);
  _vectorMemory->resize(_size);
  for (unsigned int i = 0; i < _size; i++)
  {
    (*_vectorMemory)[i].reset(new SiconosVector(*V[i]));
  }
}

// copy of a std::vector of siconos vectors  + _size
SiconosMemory::SiconosMemory(const unsigned int newMemorySize, const MemoryContainer& V):
  _size(newMemorySize), _nbVectorsInMemory(V.size())
{
  _indx = _size-1;
  _vectorMemory.reset(new MemoryContainer);
  _vectorMemory->resize(_size);
  if (newMemorySize < V.size())
    SiconosMemoryException::selfThrow("SiconosMemory(int _size, vector<SP::SiconosVector> V) : V.size > _size");
  else
  {
    for (unsigned int i = 0; i < V.size(); i++)
    {
      (*_vectorMemory)[i].reset(new SiconosVector(*V[i]));
    }
  }
}

//Copy constructor
SiconosMemory::SiconosMemory(const SiconosMemory& Mem)
{
  _size = Mem.getMemorySize();
  _indx = _size-1;
  _nbVectorsInMemory = Mem.nbVectorsInMemory();
  _vectorMemory.reset(new MemoryContainer);
  _vectorMemory->resize(_size);
  // XXX is this correct ?
  const MemoryContainer VtoCopy = *(Mem.vectorMemory());
  for (unsigned int i = 0; i < VtoCopy.size(); i++)
  {
    (*_vectorMemory)[i].reset(new SiconosVector(*VtoCopy[i]));
  }

}

// Destructor
SiconosMemory::~SiconosMemory()
{
  _vectorMemory->clear();
}

// --- GETTERS/SETTERS ---

void SiconosMemory::setVectorMemory(const MemoryContainer& V)
{
  _size = V.size();
  _indx = _size-1;
  _vectorMemory->resize(_size);
  _nbVectorsInMemory = _size;
  _vectorMemory->clear();
  for (unsigned int i = 0; i < V.size(); i++)
  {
    (*_vectorMemory)[i].reset(new SiconosVector(*V[i]));
  }
}

SP::SiconosVector SiconosMemory::getSiconosVector(const unsigned int index) const
{
  assert(index < _nbVectorsInMemory && "getSiconosVector(index) : inconsistent index value");
  return _vectorMemory->at( (_indx + 1 + index) % _size);
}

void SiconosMemory::swap(const SiconosVector& v)
{
  // In case of no memory (e.g. ZeroOrderHoldOSI), do nothing.
  if (_size == 0)
    return;
  // If vectorMemory size is _size, we remove its last element.
  *(*_vectorMemory)[_indx] = v;
  _nbVectorsInMemory = std::min(_nbVectorsInMemory+1, _size);
  if (_indx > 0)
    _indx--;
  else
    _indx = _size-1;
}

void SiconosMemory::display() const
{
  std::cout << " ====== Memory vector display ======= " <<std::endl;
  std::cout << "| _size : " << _size <<std::endl;
  std::cout << "| _nbVectorsInMemory : " << _nbVectorsInMemory <<std::endl;
  std::cout << "| vectorMemory size : " << _vectorMemory->size() <<std::endl;
  for (unsigned int i = 0; i < _nbVectorsInMemory; i++)
  {
    std::cout << "vector number " << i << ": adress = " << (*_vectorMemory)[i] << " | " <<std::endl; ;
    (*_vectorMemory)[i]->display();
  }
  std::cout << " ===================================== " <<std::endl;
}
