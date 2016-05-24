/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
  bool read(const std::string& fileName, 
            SiconosVector& m, 
            const openmode& mode, 
            int prec, 
            const std::string& inputType,
            const std::ios::fmtflags& flags
    )
  {
    // Note FP : .c_str() will be useless for std11
    std::ifstream infile(fileName.c_str(), mode);
    infile.flags(flags);
    if (!infile.good())
      SiconosVectorException::selfThrow("ioVector::read error : Fail to open \"" + fileName + "\"");
    if (infile.peek() == std::ifstream::traits_type::eof())
    {
      SiconosVectorException::selfThrow("ioVector::read : the given file is empty!");
    }

    infile.precision(prec);

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

  bool write(const std::string& fileName,
             const SiconosVector& m,
             const openmode& mode, 
             int prec, 
             const std::string& outputType,
             const std::ios::fmtflags& flags)
  {
    std::ofstream outfile(fileName.c_str(), mode);
    outfile.flags(flags);

    if (!outfile.good())
      SiconosVectorException::selfThrow("ioVector:: write error : Fail to open \"" + fileName + "\"");
    outfile.precision(prec);
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
    
      if (m.num() == 1)
      {
        DenseVect*  p = m.dense();
        std::copy(p->begin(), p->end(), std::ostream_iterator<double>(outfile, " "));
      }
      else if (m.num() == 4)
      {
        SparseVect* p = m.sparse();
        std::copy(p->begin(), p->end(), std::ostream_iterator<double>(outfile, " "));
      }
   }
    outfile.close();
    return true;
  }

}
