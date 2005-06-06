#include "SiconosMemory.h"
using namespace std;

SiconosMemory::SiconosMemory(): memorySize(0), nbVectorsInMemory(0),
  memoryXML(NULL)
{
  IN("SiconosMemory()\n");
  vectorMemory.resize(memorySize);
  OUT("SiconosMemory()\n");
}

SiconosMemory::SiconosMemory(const int& newValue):
  memorySize(newValue), nbVectorsInMemory(0), memoryXML(NULL)
{
  IN("SiconosMemory(int memorySize)\n");
  vectorMemory.clear();
  vectorMemory.resize(memorySize);
  for (int i = 0; i < memorySize; i++)
    (vectorMemory)[i] = new SimpleVector();
  OUT("SiconosMemory(int memorySize)\n");
}

SiconosMemory::SiconosMemory(const int& newMemorySize, SiconosMemoryXML *memXML):
  memorySize(newMemorySize), nbVectorsInMemory(0), memoryXML(memXML)
{
  IN("SiconosMemory(int memorySize, SiconosMemoryXML *memoryXML)\n");
  vectorMemory.clear();
  vectorMemory.resize(memorySize);
  for (int i = 0; i < memorySize; i++)
    (vectorMemory)[i] = new SimpleVector();

  if (memoryXML != NULL) setVectorMemory(memorySize, memoryXML->getSiconosMemoryVector());

  OUT("SiconosMemory(int memorySize, SiconosMemoryXML *memoryXML)\n");
}

SiconosMemory::SiconosMemory(SiconosMemoryXML *memXML):
  memorySize(0), nbVectorsInMemory(0), memoryXML(memXML)
{
  IN("SiconosMemory(SiconosMemoryXML *memoryXML)\n");
  vectorMemory.resize(memorySize);
  setVectorMemory(memoryXML->getSiconosMemorySize(), memoryXML->getSiconosMemoryVector());
  memoryXML = memoryXML;
  OUT("SiconosMemory(SiconosMemoryXML *memoryXML)\n");
}


SiconosMemory::SiconosMemory(const vector<SiconosVector*>& V):
  memorySize(V.size()), nbVectorsInMemory(memorySize), memoryXML(NULL)
{
  IN("SiconosMemory(vector<SiconosVector*> V)\n");
  vectorMemory.clear();
  vectorMemory.resize(memorySize);
  for (int i = 0; i < memorySize; i++)
  {
    (vectorMemory)[i] = new SimpleVector();
    *(vectorMemory[i]) = *(V[i]);
  }
  OUT("SiconosMemory(vector<SiconosVector*> V)\n");
}

SiconosMemory::SiconosMemory(const int& newMemorySize, const  vector<SiconosVector*>& V):
  memorySize(newMemorySize), nbVectorsInMemory(V.size()), memoryXML(NULL)
{
  IN("SiconosMemory(int memorySize, vector<SiconosVector*> V)\n");
  if (newMemorySize < V.size())
  {
    // exception : out of range
    SiconosMemoryException::selfThrow("SiconosMemory(int memorySize, vector<SiconosVector*> V) : V.size > memorySize");
  }
  else
  {
    vectorMemory.clear();
    vectorMemory.resize(memorySize);
    for (int i = 0; i < memorySize; i++)
    {
      (vectorMemory)[i] = new SimpleVector();
      if (i < nbVectorsInMemory)
        *(vectorMemory[i]) = *(V[i]);
    }
  }
  OUT("SiconosMemory(int memorySize, vector<SiconosVector*> V)\n");
}

SiconosMemory::SiconosMemory(const SiconosMemory&  source):
  memorySize(source.memorySize), nbVectorsInMemory(source.nbVectorsInMemory), memoryXML(NULL)
{
  IN("SiconosMemory(const SiconosMemory&  source) \n");
  int i;
  vectorMemory.clear();
  vectorMemory.resize(memorySize);
  for (i = 0; i < memorySize; i++)
    if (i < nbVectorsInMemory)
      (vectorMemory)[i] = new SimpleVector(*(source.vectorMemory[i]));
    else
      (vectorMemory)[i] = new SimpleVector();
  OUT("SiconosMemory(const SiconosMemory&  source) \n");
}


SiconosMemory::~SiconosMemory()
{
  IN("~SiconosMemory()\n");
  for (int i = 0; i < memorySize; i++)
    delete(vectorMemory)[i];
  vectorMemory.clear();
  OUT("~SiconosMemory()\n");
}


/* ********************************************** */

void SiconosMemory::setVectorMemory(const vector<SiconosVector*>& V)
{
  IN("SiconosMemory::setVectorMemory(vector<SiconosVector*> V)\n");

  for (int i = 0; i < memorySize; i++)
    delete(vectorMemory)[i];

  memorySize = V.size();
  nbVectorsInMemory = memorySize;
  vectorMemory.clear();
  vectorMemory.resize(memorySize);
  for (int i = 0; i < memorySize; i++)
  {
    (vectorMemory)[i] = new SimpleVector();
    *(vectorMemory[i]) = *(V[i]);
  }

  OUT("SiconosMemory::setVectorMemory(vector<SiconosVector*> V)\n");
}

void SiconosMemory::setVectorMemory(const int& mSize, const vector<SiconosVector*>& V)
{
  IN("SiconosMemory::setVectorMemory(int mSize, vector<SiconosVector*> V)\n");

  if (mSize < V.size())
  {
    // exception : out of range
    //SiconosMemoryException::selfThrow("setVectorMemory(int mSize, vector<SiconosVector*> V) : V.size > mSize");
    cout << "setVectorMemory(int mSize, vector<SiconosVector*> V) : V.size > mSize.\nThe saved SiconosVector which position in the vector of SiconosVector is greater than mSize, will be lost." << endl;

    /*
     * the node in the DOM tree which are unused will be deleted
     */
    if (memoryXML != NULL)
      memoryXML->deleteUnusedMemoryNodes(mSize);
    else SiconosMemoryException::selfThrow("setVectorMemory(int mSize, vector<SiconosVector*> V) : V.size > mSize");
  }
  else
  {
    int i;

    if (memorySize != mSize)
    {
      for (i = 0; i < memorySize; i++)
        delete(vectorMemory)[i];
      memorySize = mSize;

      vectorMemory.clear();
      vectorMemory.resize(memorySize);
      for (i = 0; i < memorySize; i++)
        (vectorMemory)[i] = new SimpleVector();
    }

    nbVectorsInMemory = V.size();
    for (i = 0; i < memorySize; i++)
    {
      if (i < nbVectorsInMemory)
        *(vectorMemory[i]) = *(V[i]);
    }
  }

  OUT("SiconosMemory::setVectorMemory(int mSize, vector<SiconosVector*> V)\n");
}


SiconosVector* SiconosMemory::getSiconosVector(const int& index) const
{
  IN("SiconosMemory::getSiconosVector(int index) const\n");

  if ((index < 0) || (index >= memorySize))
  {
    // exception : out of range
    SiconosMemoryException::selfThrow("getSiconosVector(int index) : index > memorySize");
  }
  else if (index >= nbVectorsInMemory)
  {
    // Warning : not significant value
    cout << "WARNING --- SiconosMemory::getSiconosVector(int index) : index > nbVectorsInMemory" << endl;
    OUT("SiconosMemory::getSiconosVector(int index) const\n");
    return NULL;
  }
  else
  {
    OUT("SiconosMemory::getSiconosVector(int index) const\n");
    return (vectorMemory)[index];
  }
}


void SiconosMemory::swap(SiconosVector* v)
{
  IN("SiconosMemory::swap(SiconosVector* v)\n");
  int i;
  SiconosVector* tmp;

  // we add always in head of vectorMemory
  // direct swapping
  if (nbVectorsInMemory < memorySize)
  {
    nbVectorsInMemory ++;
  }
  else if (memorySize <= 0)
  {
    // exception : out of range
    SiconosMemoryException::selfThrow("swap(SiconosVector* v) : memorySize <= 0");
  }

  tmp = (vectorMemory)[nbVectorsInMemory - 1];

  for (i = nbVectorsInMemory - 1; i >= 0; i--)
  {
    if (i == 0)
    {
      // where we put v (we copy it)
      (vectorMemory)[i] = tmp;
      *(vectorMemory)[i] = *v;
    }
    else
    {
      // we move the old vectors in memory
      (vectorMemory)[i] = (vectorMemory)[i - 1];
    }
  }

  OUT("SiconosMemory::swap(SiconosVector* v)\n");
}



void SiconosMemory::display() const
{
  cout << "| memorySize : " << memorySize << endl;
  cout << "| nbVectorsInMemory : " << nbVectorsInMemory << endl;
  cout << "| vectorMemory size : " << vectorMemory.size() << endl;
  for (int i = 0; i < nbVectorsInMemory; i++)
    //for (int i = 0; i < memorySize; i++)
  {
    cout << "Memory " << i << " >>> " << (vectorMemory)[i] << " | ";
    vectorMemory[i]->display();
  }
}


SiconosMemory& SiconosMemory::operator = (const SiconosMemory& source)
{
  IN("SiconosMemory::operator = (const SiconosMemory& source)\n");

  int i;

  for (i = 0; i < memorySize; i++)
    delete(vectorMemory)[i];

  memorySize = source.memorySize;

  vectorMemory.clear();
  vectorMemory.resize(memorySize);
  for (i = 0; i < memorySize; i++)
    (vectorMemory)[i] = new SimpleVector();

  nbVectorsInMemory = source.nbVectorsInMemory;
  for (i = 0; i < memorySize; i++)
  {
    if (i < nbVectorsInMemory)
      *(vectorMemory[i]) = *(source.vectorMemory[i]);
  }
  memoryXML = source.memoryXML;

  OUT("SiconosMemory::operator = (const SiconosMemory& source)\n");
  return *this;
}


void SiconosMemory::fillMemoryWithMemoryXML()
{
  IN("SiconosMemory::fillMemoryWithMemoryXML()\n");
  int size = memoryXML->getSiconosMemorySize();
  vector<SiconosVector*> v = memoryXML->getSiconosMemoryVector();
  setVectorMemory(size, v);
  OUT("SiconosMemory::fillMemoryWithMemoryXML()\n");
}


