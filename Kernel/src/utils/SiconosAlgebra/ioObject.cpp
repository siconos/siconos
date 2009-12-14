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
#include "ioObject.hpp"
#include "RuntimeException.hpp"

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

