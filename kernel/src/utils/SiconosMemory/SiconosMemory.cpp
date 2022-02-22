/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

// From the size of the container (number of saved vectors) and
// the size of the vectors.
SiconosMemory::SiconosMemory(const unsigned int size, const unsigned int vectorSize): MemoryContainer(), _indx(size-1)
{
  for(unsigned int i = 0; i < size; i++)
  {
    push_back(SiconosVector(vectorSize));
  }
}

//Copy constructor
SiconosMemory::SiconosMemory(const SiconosMemory& Mem)
  : MemoryContainer(),
    _nbVectorsInMemory(Mem.nbVectorsInMemory()),
    _indx(Mem.size()-1)
{
  for(unsigned int i = 0; i < Mem.size(); i++)
  {
    push_back(Mem[i]);
  }
}

void SiconosMemory::operator=(const SiconosMemory& V)
{
  if(size() != V.size())
  {
    this->resize(V.size()); // => copy construction of old content
  }

  for(unsigned int i = 0; i < V.size(); i++)
  {
    (*this)[i].resize(V[i].size(), true);
    (*this)[i] = V[i]; // copy
  }
  _indx = V._indx;
  _nbVectorsInMemory = V._nbVectorsInMemory;
}


// (Re)set the size of an existing SiconosMemory
void SiconosMemory::setMemorySize(const unsigned int steps,
                                  const unsigned int vectorSize)
{
  _nbVectorsInMemory = 0;
  _indx = steps-1;
  for(unsigned int i = 0; i < size(); i++)
  {
    this->at(i).resize(vectorSize, true);
  }
  for(unsigned int i = size(); i < steps; i++)
  {
    this->push_back(SiconosVector(vectorSize));
  }
}

// --- GETTERS/SETTERS ---

const SiconosVector& SiconosMemory::getSiconosVector(const unsigned int index) const
{
  assert(index < _nbVectorsInMemory && "getSiconosVector(index) : inconsistent index value");
  return this->at((_indx + 1 + index) % this->size());
}

SiconosVector& SiconosMemory::getSiconosVectorMutable(const unsigned int index)
{
  assert(index < _nbVectorsInMemory && "getSiconosVector(index) : inconsistent index value");
  return *(SiconosVector*)(&this->at((_indx + 1 + index) % this->size()));
}

void SiconosMemory::swap(const SiconosVector& v)
{
  // Be robust to empty memory
  if(size()==0)
    return;

  // If _nbVectorsInMemory is this->size(), we remove the last element.
  (*this)[_indx] = v;
  _nbVectorsInMemory = std::min(_nbVectorsInMemory+1, this->size());
  if(_indx > 0)
    _indx--;
  else
    _indx = this->size()-1;
}

void SiconosMemory::swap(SP::SiconosVector v)
{
  // Be robust to empty memory
  // Be robust to null pointer
  if(size()==0 || !v)
    return;

  // If _nbVectorsInMemory is this->size(), we remove the last element.
  (*this)[_indx] = *v;
  _nbVectorsInMemory = std::min(_nbVectorsInMemory+1, this->size());
  if(_indx > 0)
    _indx--;
  else
    _indx = this->size()-1;
}

void SiconosMemory::display() const
{
  std::cout << " ====== Memory vector display ======= " <<std::endl;
  std::cout << "| size : " << this->size() <<std::endl;
  std::cout << "| _nbVectorsInMemory : " << _nbVectorsInMemory <<std::endl;
  for(unsigned int i = 0; i < _nbVectorsInMemory; i++)
  {
    std::cout << "vector number " << i << ": address = "
              << &this->at(i) << " | " << std::endl;
    this->at(i).display();
  }
  std::cout << " ===================================== " <<std::endl;
}
