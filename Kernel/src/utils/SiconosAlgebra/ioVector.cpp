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
#include "ioVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorException.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include<fstream>

#include "SiconosAlgebra.hpp"

namespace ioVector
{
  bool read(const std::string& fileName, SiconosVector& m, const openmode& mode, const std::string& inputType)
  {
    // Note FP : .c_str() will be useless for std11
    std::ifstream infile(fileName.c_str(), mode);
    if (!infile.good())
      SiconosVectorException::selfThrow("ioVector::read error : Fail to open \"" + fileName + "\"");
    infile.precision(15);

    if(mode == BINARY_IN)
    {
      double * x =  m.getArray();
      if(inputType != "noDim")
      {
        unsigned int dim;
        infile.read((char*)(&dim), sizeof(m.size()));
      
      }
      infile.read((char*)(&x[0]), m.size() * sizeof (double));
    }
    else
    {
      DenseVect *p = m.dense();
      // Read the dimension of the vector in the first line of the input file
      // Just use to check that sizes are consistents.
      if(inputType != "noDim")
      {
        unsigned int dim;
        infile >> dim;
        if (dim != p->size())
          p->resize(dim);
      }
      copy((std::istream_iterator<double>(infile)), std::istream_iterator<double>(), (p->data()).begin());
    }
    infile.close();
    return true;
  }

  bool write(const std::string& fileName, const SiconosVector& m, const openmode& mode, const std::string& outputType)
  {
    std::ofstream outfile(fileName.c_str(), mode);
    if (!outfile.good())
      SiconosVectorException::selfThrow("ioVector:: write error : Fail to open \"" + fileName + "\"");
    outfile.precision(15);
    if(mode == BINARY_OUT)
    {
      double * x = m.getArray();
      if(outputType != "noDim")
      {
        unsigned int dim = m.size();
        outfile.write((char*)&dim, sizeof (dim));
      }
      outfile.write ((char*)(&x[0]), sizeof (double) * m.size());

    }
    else
    {
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
   }
    outfile.close();
    return true;
  }

}
