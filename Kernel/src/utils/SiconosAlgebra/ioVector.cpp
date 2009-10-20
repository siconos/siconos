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
#include "ioVector.h"
#include "SimpleVector.h"
#include "SiconosVectorException.h"
#include <boost/numeric/ublas/io.hpp>
#include<fstream>

// Default private
ioVector::ioVector(): ioObject() {}

ioVector::ioVector(const std::string& file, const std::string& mode): ioObject(file, mode) {}

ioVector::~ioVector(void) {}

const bool ioVector::read(SiconosVector& m) const
{

  std::ifstream infile;
  if (Mode == "ascii")
    infile.open(FileName.c_str(), std::ifstream::in);
  else if (Mode == "binary")
    infile.open(FileName.c_str(), std::ifstream::binary);
  else
    SiconosVectorException::selfThrow(" ioVector::read : Fail to open file \"" + FileName + "\"");

  if (!infile.good())
    SiconosVectorException::selfThrow("ioVector::read error : Fail to open \"" + FileName + "\"");

  if (m.isBlock())
    SiconosVectorException::selfThrow(" ioVector::read : read BlockVector is not implemented");

  infile.precision(15);

  DenseVect *p = m.dense();

  // Read the dimension of the vector in the first line of the input file
  // Just use to check that sizes are consistents.
  unsigned int dim;
  infile >> dim;

  if (dim != p->size())
    p->resize(dim);

  copy((std::istream_iterator<double>(infile)), std::istream_iterator<double>(), (p->data()).begin());

  infile.close();
  return true;
}

const bool ioVector::write(const SiconosVector& m, const std::string& outputType) const
{
  std::ofstream outfile;
  if (Mode == "ascii")
    outfile.open(FileName.c_str(), std::ofstream::out);
  else if (Mode == "binary")
    outfile.open(FileName.c_str(), std::ofstream::binary);
  else
    SiconosVectorException::selfThrow("ioVector::write - Incorrect mode for writing");

  if (!outfile.good())
    SiconosVectorException::selfThrow("ioVector:: write error : Fail to open \"" + FileName + "\"");

  if (m.isBlock())
    SiconosVectorException::selfThrow(" ioVector::write : write BlockVector is not implemented");

  outfile.precision(15);

  if (outputType != "noDim")
    outfile << m.size() << std::endl;

  if (m.getNum() == 1)
  {
    DenseVect*  p = m.dense();
    std::copy(p->begin(), p->end(), std::ostream_iterator<double>(outfile, " "));
  }
  else if (m.getNum() == 4)
  {
    SparseVect* p = m.sparse();
    std::copy(p->begin(), p->end(), std::ostream_iterator<double>(outfile, " "));
  }

  outfile.close();
  return true;
}

