/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "SiconosMemory.h"
#include "BlockVector.h"
#include "SimpleVector.h"

using namespace std;

// --- CONSTRUCTORS ---


// from data: memorySize
SiconosMemory::SiconosMemory(const unsigned int newValue):
  memorySize(newValue), nbVectorsInMemory(0), memoryXML(NULL)
{}

// from xml file + optional value of memorySize
SiconosMemory::SiconosMemory(SiconosMemoryXML *memXML, const unsigned int newMemorySize):
  memorySize(newMemorySize), nbVectorsInMemory(0), memoryXML(memXML)
{
  if (memoryXML != NULL)
  {
    // Convention: memorySize==1 (default value) means read its value in xml file (only for this constructor)
    // if hasMemory() is true.
    if (memorySize == 1 && memoryXML->hasMemory())
      memorySize = memoryXML->getSiconosMemorySize();

    // get memory from xml file
    deque<SiconosVector*> V;
    bool vAllocated = false; // true if memory is allocated for vector reading in SiconosMemoryXML::getSiconosMemoryVector()
    if (memoryXML->hasMemory())
    {
      V =  memoryXML->getSiconosMemoryVector();
      vAllocated = true;
    }
    // if size of V overpass memorysize:
    if (memorySize < V.size())
    {
      cout << "Warning: xml SiconosMemory constructor, size of V(xml) greater than memory size, excess values will be lost" << endl;
      // unused nodes of the DOM tree are deleted
      memoryXML->deleteUnusedMemoryNodes(memorySize);
      nbVectorsInMemory = memorySize;
    }
    else
      nbVectorsInMemory = V.size();

    unsigned int i;
    for (i = 0; i < nbVectorsInMemory ; i++)
    {
      if (V[i]->isBlock())
        vectorMemory.push_back(new BlockVector(*V[i]));

      else
        vectorMemory.push_back(new SimpleVector(*V[i]));
      isVectorMemoryAllocated.push_back(true);
    }

    if (vAllocated)
    {
      deque<SiconosVector*>::iterator it;
      for (it = V.begin(); it != V.end(); it++)
        delete(*it);
    }
  }
  else
    SiconosMemoryException::selfThrow("SiconosMemory, xml constructor: xml file==NULL");
}

// copy of a std::vector of siconos vectors
SiconosMemory::SiconosMemory(const deque<SiconosVector*>& V):
  memorySize(V.size()), nbVectorsInMemory(V.size()), memoryXML(NULL)
{
  unsigned int sizeV = V.size();

  memorySize = sizeV;
  nbVectorsInMemory = sizeV;
  for (unsigned int i = 0; i < sizeV; i++)
  {
    if (V[i]->isBlock())
      vectorMemory.push_back(new BlockVector(*V[i]));
    else
      vectorMemory.push_back(new SimpleVector(*V[i]));
    isVectorMemoryAllocated.push_back(true);
  }
}

// copy of a std::vector of siconos vectors  + memorySize
SiconosMemory::SiconosMemory(const unsigned int newMemorySize, const  deque<SiconosVector*>& V):
  memorySize(newMemorySize), nbVectorsInMemory(V.size()), memoryXML(NULL)
{
  if (newMemorySize < V.size())
    SiconosMemoryException::selfThrow("SiconosMemory(int memorySize, vector<SiconosVector*> V) : V.size > memorySize");
  else
  {
    for (unsigned int i = 0; i < V.size(); i++)
    {
      if (V[i]->isBlock())
        vectorMemory.push_back(new BlockVector(*V[i]));
      else
        vectorMemory.push_back(new SimpleVector(*V[i]));
      isVectorMemoryAllocated.push_back(true);
    }
  }
}

// copy
SiconosMemory::SiconosMemory(const SiconosMemory&  source):
  memorySize(source.memorySize), nbVectorsInMemory(source.nbVectorsInMemory), memoryXML(NULL)
{
  for (unsigned int i = 0; i < nbVectorsInMemory; i++)
  {
    if (source.vectorMemory[i]->isBlock())
      vectorMemory.push_back(new BlockVector(*(source.vectorMemory[i])));
    else
      vectorMemory.push_back(new SimpleVector(*(source.vectorMemory[i])));
    isVectorMemoryAllocated.push_back(true);
  }
}

// Destructor
SiconosMemory::~SiconosMemory()
{
  for (unsigned int i = 0; i < nbVectorsInMemory; i++)
  {
    if (isVectorMemoryAllocated[i])
    {
      delete vectorMemory[i];
      vectorMemory[i] = 0;
    }
  }
  vectorMemory.clear();
  isVectorMemoryAllocated.clear();
}

// --- GETTERS/SETTERS ---

void SiconosMemory::setVectorMemory(const deque<SiconosVector*>& V)
{
  for (unsigned int i = 0; i < vectorMemory.size(); i++)
    if (isVectorMemoryAllocated[i])
    {
      delete vectorMemory[i];
      vectorMemory[i] = 0;
    }

  unsigned int sizeV = V.size();

  memorySize = sizeV;
  nbVectorsInMemory = sizeV;
  vectorMemory.clear();
  isVectorMemoryAllocated.clear();
  for (unsigned int i = 0; i < V.size(); i++)
  {
    if (V[i]->isBlock())
      vectorMemory.push_back(new BlockVector(*V[i]));
    else
      vectorMemory.push_back(new SimpleVector(*V[i]));
    isVectorMemoryAllocated.push_back(true);
  }
}

SiconosVector* SiconosMemory::getSiconosVector(const unsigned int index) const
{
  if (index >= nbVectorsInMemory)
    SiconosMemoryException::selfThrow("getSiconosVector(index) : inconsistent index value");
  return vectorMemory[index];
}


void SiconosMemory::swap(SiconosVector* v)
{
  unsigned int i;
  SiconosVector* tmp;
  double tmp2;

  // if it remains space in the vector
  if (nbVectorsInMemory < memorySize)
  {
    // allocate memory for the new vector
    if (v->isBlock())
      vectorMemory.push_front(new BlockVector(*v));
    else
      vectorMemory.push_front(new SimpleVector(*v));
    isVectorMemoryAllocated.push_front(true);
    nbVectorsInMemory ++;
  }
  else
  {
    *vectorMemory[nbVectorsInMemory - 1] = *v;
    // Permutations to reorganise the vector
    tmp = vectorMemory[nbVectorsInMemory - 1];
    tmp2 = isVectorMemoryAllocated[nbVectorsInMemory - 1];
    for (i = nbVectorsInMemory - 1; i > 0; i--)
    {
      vectorMemory[i] = vectorMemory[i - 1];
      isVectorMemoryAllocated[i] = isVectorMemoryAllocated[i - 1];
    }
    // copy of v into the first position
    vectorMemory[0] = tmp;
    isVectorMemoryAllocated[0] = tmp2;
  }
}



void SiconosMemory::display() const
{
  cout << " ====== Memory vector display ======= " << endl;
  cout << "| memorySize : " << memorySize << endl;
  cout << "| nbVectorsInMemory : " << nbVectorsInMemory << endl;
  cout << "| vectorMemory size : " << vectorMemory.size() << endl;
  for (unsigned int i = 0; i < nbVectorsInMemory; i++)
  {
    cout << "vector number " << i << ": adress = " << vectorMemory[i] << " | " << endl; ;
    vectorMemory[i]->display();
  }
  cout << " ===================================== " << endl;
}

// operator =
// Warning: memory vectors should have the same memorySize
// Warning: no memory allocation is required.
// If necessary, use copy constructor instead
SiconosMemory& SiconosMemory::operator = (const SiconosMemory& source)
{
  // error if vector have not the same size
  if (memorySize != source.memorySize)
    SiconosMemoryException::selfThrow("SiconosMemory, operator =, vectors have not the same size.");

  // clean memory before assignment
  for (unsigned int i = 0; i < nbVectorsInMemory ; i++)
    if (isVectorMemoryAllocated[i])
    {
      delete vectorMemory[i];
      vectorMemory[i] = 0;
      isVectorMemoryAllocated[i] = false;
    }

  isVectorMemoryAllocated.resize(source.isVectorMemoryAllocated.size());
  vectorMemory.resize(source.vectorMemory.size());
  nbVectorsInMemory = source.nbVectorsInMemory;

  // !! no memory allocation !!
  vectorMemory = source.vectorMemory;
  // Warning: do not copy isVectorMemoryAllocatedIn to avoid double memory deallocation
  memoryXML = source.memoryXML;

  return *this;
}

// default (private) constructor
SiconosMemory::SiconosMemory(): memorySize(0), nbVectorsInMemory(0),
  memoryXML(NULL)
{}
