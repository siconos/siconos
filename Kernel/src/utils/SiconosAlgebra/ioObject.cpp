#include "ioObject.h"
#include "RuntimeException.h"

ioObject::ioObject(const std::string& m): FileName("NoName.dat"), Mode(m)
{}

ioObject::ioObject(const std::string& file, const std::string& m): FileName(file), Mode(m) {}

ioObject::~ioObject(void) {};

const bool ioObject::read(SiconosMatrix& m)const
{
  RuntimeException::selfThrow("ioObject::read(SiconosMatrix) - not implemented");
  return false;
}

const bool ioObject::write(const SiconosMatrix& m, const std::string&) const
{
  RuntimeException::selfThrow("ioObject::write(SiconosMatrix) - not implemented");
  return false;
}

const bool ioObject::read(SiconosVector& m) const
{
  RuntimeException::selfThrow("ioObject::read(SiconosVector&) - not implemented");
  return false;
}

const bool ioObject::write(const SiconosVector& m, const std::string&) const
{
  RuntimeException::selfThrow("ioObject::write(SiconosVector&) - not implemented");
  return false;
}

