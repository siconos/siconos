/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "SiconosMemory.hpp"
#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosMemoryXML.hpp"



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

// from xml file + optional value of _size
SiconosMemory::SiconosMemory(SP::SiconosMemoryXML memXML):
  _size(0), _nbVectorsInMemory(0), _memoryXML(memXML)
{

  if (!_memoryXML)
    SiconosMemoryException::selfThrow("SiconosMemory, xml constructor: xml file==NULL");

  _size = _memoryXML->getSiconosMemorySize();
  _vectorMemory.reset(new MemoryContainer);
  _vectorMemory->resize(_size);

  if (!_memoryXML->hasMemory())
    SiconosMemoryException::selfThrow("SiconosMemory, xml constructor: no memory node found.");

  // get memory from xml file
  _vectorMemory =  _memoryXML->getSiconosMemoryVector();
  _nbVectorsInMemory = _vectorMemory->size();
  _indx = _size-1;
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
  // this was not always initialised
  if (Mem.getSiconosMemoryXML())
    _memoryXML.reset(new SiconosMemoryXML(*(Mem.getSiconosMemoryXML())));

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

void SiconosMemory::saveMemorySizeToXML()
{
  if (_memoryXML) _memoryXML->setSiconosMemorySize(_size);
  else SiconosMemoryException::selfThrow("SiconosMemory::saveMemorySizeToXML() - _memoryXML object == NULL");
}
