#include "SiconosVector.h"
#include "SiconosVectorException.h"

// Default (private) constructor
SiconosVector::SiconosVector(bool isblock): isBlockVector(isblock), sizeV(0)
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

void SiconosVector::xpy(const SiconosVector &, const SiconosVector &)
{
  SiconosVectorException::selfThrow("SiconosVector::xpy() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
}

void SiconosVector::scal(double, const SiconosVector&)
{
  SiconosVectorException::selfThrow("SiconosVector::scal() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
}


void SiconosVector::axpy(double, const SiconosVector&, const SiconosVector&)
{
  SiconosVectorException::selfThrow("SiconosVector::axpy() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
}


void SiconosVector::axpby(double, const SiconosVector&, double, const SiconosVector&)
{
  SiconosVectorException::selfThrow("SiconosVector::axpby() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
}

