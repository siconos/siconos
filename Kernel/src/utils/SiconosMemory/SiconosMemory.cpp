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
#include "SiconosMemory.h"
#include "BlockVector.h"
#include "SimpleVector.h"

using namespace std;

// --- CONSTRUCTORS ---


// from data: memorySize
SiconosMemory::SiconosMemory(const unsigned int newValue):
  memorySize(newValue), nbVectorsInMemory(0)
{}

// from xml file + optional value of memorySize
SiconosMemory::SiconosMemory(SP::SiconosMemoryXML memXML, const unsigned int newMemorySize):
  memorySize(newMemorySize), nbVectorsInMemory(0), memoryXML(memXML)
{
  if (memoryXML)
  {
    // Convention: memorySize==1 (default value) means read its value in xml file (only for this constructor)
    // if hasMemory() is true.
    if (memorySize == 1 && memoryXML->hasMemory())
      memorySize = memoryXML->getSiconosMemorySize();

    // get memory from xml file
    deque<SP::SiconosVector> V;
    if (memoryXML->hasMemory())
    {
      V =  memoryXML->getSiconosMemoryVector();
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
        vectorMemory.push_back(SP::BlockVector(new BlockVector(*V[i])));

      else
        vectorMemory.push_back(SP::SimpleVector(new SimpleVector(*V[i])));
    }

  }
  else
    SiconosMemoryException::selfThrow("SiconosMemory, xml constructor: xml file==NULL");
}

// copy of a std::vector of siconos vectors
SiconosMemory::SiconosMemory(const deque<SP::SiconosVector>& V):
  memorySize(V.size()), nbVectorsInMemory(V.size())
{
  unsigned int sizeV = V.size();

  memorySize = sizeV;
  nbVectorsInMemory = sizeV;
  for (unsigned int i = 0; i < sizeV; i++)
  {
    if (V[i]->isBlock())
      vectorMemory.push_back(SP::BlockVector(new BlockVector(*V[i])));
    else
      vectorMemory.push_back(SP::SimpleVector(new SimpleVector(*V[i])));
  }
}

// copy of a std::vector of siconos vectors  + memorySize
SiconosMemory::SiconosMemory(const unsigned int newMemorySize, const  deque<SP::SiconosVector>& V):
  memorySize(newMemorySize), nbVectorsInMemory(V.size())
{
  if (newMemorySize < V.size())
    SiconosMemoryException::selfThrow("SiconosMemory(int memorySize, vector<SP::SiconosVector> V) : V.size > memorySize");
  else
  {
    for (unsigned int i = 0; i < V.size(); i++)
    {
      if (V[i]->isBlock())
        vectorMemory.push_back(SP::BlockVector(new BlockVector(*V[i])));
      else
        vectorMemory.push_back(SP::SimpleVector(new SimpleVector(*V[i])));
    }
  }
}

// copy
SiconosMemory::SiconosMemory(const SiconosMemory&  source):
  memorySize(source.memorySize), nbVectorsInMemory(source.nbVectorsInMemory)
{
  for (unsigned int i = 0; i < nbVectorsInMemory; i++)
  {
    if (source.vectorMemory[i]->isBlock())
      vectorMemory.push_back(SP::BlockVector(new BlockVector(*(source.vectorMemory[i]))));
    else
      vectorMemory.push_back(SP::SimpleVector(new SimpleVector(*(source.vectorMemory[i]))));
  }
}

// Destructor
SiconosMemory::~SiconosMemory()
{
  vectorMemory.clear();
}

// --- GETTERS/SETTERS ---

void SiconosMemory::setVectorMemory(const deque<SP::SiconosVector>& V)
{
  unsigned int sizeV = V.size();

  memorySize = sizeV;
  nbVectorsInMemory = sizeV;
  vectorMemory.clear();
  for (unsigned int i = 0; i < V.size(); i++)
  {
    if (V[i]->isBlock())
      vectorMemory.push_back(SP::BlockVector(new BlockVector(*V[i])));
    else
      vectorMemory.push_back(SP::SimpleVector(new SimpleVector(*V[i])));
  }
}

SP::SiconosVector SiconosMemory::getSiconosVector(const unsigned int index) const
{
  assert(index < nbVectorsInMemory &&
         "getSiconosVector(index) : inconsistent index value");
  return vectorMemory[index];
}


void SiconosMemory::swap(SP::SiconosVector v)
{
  unsigned int i;
  SP::SiconosVector tmp;
  double tmp2;

  // if it remains space in the vector
  if (nbVectorsInMemory < memorySize)
  {
    // allocate memory for the new vector
    if (v->isBlock())
      vectorMemory.push_front(SP::BlockVector(new BlockVector(*v)));
    else
      vectorMemory.push_front(SP::SimpleVector(new SimpleVector(*v)));
    nbVectorsInMemory ++;
  }
  else
  {
    *vectorMemory[nbVectorsInMemory - 1] = *v;
    // Permutations to reorganise the vector
    tmp = vectorMemory[nbVectorsInMemory - 1];
    for (i = nbVectorsInMemory - 1; i > 0; i--)
    {
      vectorMemory[i] = vectorMemory[i - 1];
    }
    // copy of v into the first position
    vectorMemory[0] = tmp;
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

  vectorMemory.resize(source.vectorMemory.size());
  nbVectorsInMemory = source.nbVectorsInMemory;

  // !! no memory allocation !!
  vectorMemory = source.vectorMemory;
  // Warning: do not copy isVectorMemoryAllocatedIn to avoid double memory deallocation
  memoryXML = source.memoryXML;

  return *this;
}

// default (private) constructor
SiconosMemory::SiconosMemory(): memorySize(0), nbVectorsInMemory(0)
{}
