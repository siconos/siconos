/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include <iterator>
#include <iostream>
#include <fstream>
#include <limits>
#include "ioMatrix.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosException.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "SiconosAlgebra.hpp"

namespace ioMatrix
{
bool read(const std::string& fileName, const std::string& mode, SiconosMatrix& m)
{
  std::ifstream infile;
  if(mode == "ascii")
    infile.open(fileName.c_str(), std::ifstream::in);
  else if(mode == "binary")
    infile.open(fileName.c_str(), std::ifstream::binary);
  else
    THROW_EXCEPTION("incorrect mode for reading");

  if(!infile.good())
    THROW_EXCEPTION("");

  if(infile.peek() == std::ifstream::traits_type::eof())
  {
    THROW_EXCEPTION("the given file is empty!");
  }

  if(m.isBlock())
    THROW_EXCEPTION("not yet implemented for block matrix.");

  infile.precision(15);
  infile.setf(std::ios::scientific);
  DenseMat& p = *m.dense();

  // Dim of the matrix are given in the first line.
  // Just use to check that sizes are consistents.

  unsigned int s1, s2;
  infile >> s1;
  infile >> s2;

  if(s1 != p.size1() || s2 != p.size2())
    p.resize(s1, s2);

  // Note: using istream stl iterator seems to be 2-times faster than << with a loop over matrix data.
  //  copy((std::istream_iterator<double>(infile)), std::istream_iterator<double>(), (p->data()).begin());
  // But it fails with column-major saving ... (ok if user write its matrix in a column-major way)

  DenseMat::iterator1 it;
  DenseMat::iterator2 it2;
  //    std::cout.precision(15);
  //    std::cout.setf(std::ios::scientific);

  for(unsigned int i = 0; i < s1; i++)
  {
    for(unsigned int j = 0; j < s2; j++)
    {
      infile >> p(i, j);
      /* fail on ubuntu 14.04 assert(infile.good());*/
    }
  }

  // Old version: result in Boost format for ouptut
  //  DenseMat * p = m.dense();
  //  infile >> *p;

  infile.close();
  return true;
}

bool write(const std::string& fileName, const std::string& mode, const SiconosMatrix& m, const std::string& outputType)
{
  // Open file and various checks
  std::ofstream outfile;
  if(mode == "ascii")
    outfile.open(fileName.c_str(), std::ofstream::out);
  else if(mode == "binary")
    outfile.open(fileName.c_str(), std::ofstream::binary);
  else
    THROW_EXCEPTION("Incorrect mode for writing");

  if(!outfile.good())
    THROW_EXCEPTION("");

  if(m.isBlock())
    THROW_EXCEPTION("not yet implemented for BlockMatrix");

  outfile.precision(15);
  outfile.setf(std::ios::scientific);
  // Writing

  if(outputType != "noDim")
    outfile << m.size(0) << " " << m.size(1) << std::endl;

  if(m.num() == Siconos::DENSE)
  {
    // DenseMat * p = m.dense();
    DenseMat::iterator1 row;
    DenseMat::iterator2 col;
    double tmp;
    for(unsigned int i = 0; i < m.size(0); i++)
    {
      for(unsigned int j = 0; j < m.size(1); j++)
      {
        tmp = m(i, j);
        if(fabs(tmp) < std::numeric_limits<double>::min()) tmp = 0.0;
        outfile << tmp << " " ;
        assert(outfile.good());
      }
      outfile << std::endl;

    }
  }
  else if(m.num() == Siconos::TRIANGULAR)
  {
    TriangMat * p = m.triang();
    TriangMat::iterator1 row;
    for(row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else if(m.num() == Siconos::SYMMETRIC)
  {
    SymMat * p = m.sym();
    SymMat::iterator1 row;
    for(row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else if(m.num() == Siconos::SPARSE)
  {
    SparseMat * p = m.sparse();
    SparseMat::iterator1 row;
    for(row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else
  {
    BandedMat * p = m.banded();
    BandedMat::iterator1 row;
    for(row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }

  outfile.close();
  return true;
}

double compareRefFile(const SimpleMatrix& data, std::string filename, double epsilon,
                      Index index, SP::SimpleMatrix* ref, std::string mode,
                      bool verbose)
{
  SP::SimpleMatrix r;
  if(!ref) ref = &r;
  *ref = std::make_shared<SimpleMatrix>(data);
  (*ref)->zero();
  bool compare = false;

  try
  {
    compare = ioMatrix::read(filename, mode, **ref);
  }
  catch(...)
  {
    if(verbose)
      std::cout << "Warning: reference file " << filename
                << " not found, no comparison performed." << std::endl;
    Siconos::exception::process();
  }
  if(!compare)
    return -1.0;

  if(verbose)
    std::cout << "Comparison with reference file " << filename << std::endl;

  SP::SiconosVector err(new SiconosVector(data.size(1)));
  (data - **ref).normInfByColumn(err);

  if(verbose)
    err->display();

  if(index.size()==0)
    for(unsigned int i = 0; i < err->size(); ++i)
      index.push_back(i);

  /* Scalar error = max of columns */
  double error = 0.0;
  for(unsigned int i = 0; i < index.size(); ++i)
  {
    if(error < (*err)(index[i]))
      error = (*err)(index[i]);
  }

  if(verbose)
    std::cout << "Error = " << error << std::endl;
  if(error > epsilon)
  {
    if(verbose)
    {
      std::cout << "Warning. The results are rather different from the reference file." << std::endl;
      std::cout << "Error = "<< error << std::endl;
    }
  }

  return error;
}

} // namespace ioMatrix

// To be used later ... ?
//   template <class T, class iterator1> friend void write( const T& obj, iterator1 row, std::ofstream outfile)
//   {
//       for(row = obj.begin1(); row!=obj.end1() ; ++row)
//  {
//    std::copy(row.begin(),row.end(),std::ostream_iterator<double>(outfile," "));
//    outfile << std::endl;
//  }
//     }
