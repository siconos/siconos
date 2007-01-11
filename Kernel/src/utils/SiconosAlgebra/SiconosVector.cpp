#include "SiconosVector.h"
#include "SiconosVectorException.h"

// Default (private) constructor
SiconosVector::SiconosVector(bool isblock): isBlockVector(isblock)
{}

SiconosVector::~SiconosVector() {}

Index SiconosVector::getTabIndex() const
{
  SiconosVectorException::selfThrow("SiconosVector::getTabIndex() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
  // fake to avoid error on warning.
  Index tmp;
  return tmp;
}

void SiconosVector::add(const  SiconosVector&)
{
  SiconosVectorException::selfThrow("SiconosVector::add() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
}

void SiconosVector::addPtr(SiconosVector*)
{
  SiconosVectorException::selfThrow("SiconosVector::addPtr() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
}
