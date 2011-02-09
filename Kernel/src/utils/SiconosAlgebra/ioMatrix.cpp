/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

#include <iterator>
#include<iostream>
#include<fstream>
#include<limits>
#include "ioMatrix.hpp"
#include "SiconosMatrixException.hpp"
#include <boost/numeric/ublas/io.hpp>

template<> bool ioMatrix::read(SiconosMatrix& m) const
{
  std::ifstream infile;
  if (Mode == "ascii")
    infile.open(FileName.c_str(), std::ifstream::in);
  else if (Mode == "binary")
    infile.open(FileName.c_str(), std::ifstream::binary);
  else
    SiconosMatrixException::selfThrow("ioMatrix::read Incorrect mode for reading");

  if (!infile.good())
    SiconosMatrixException::selfThrow("ioMatrix::read error : Fail to open \"" + FileName + "\"");

  if (m.isBlock())
    SiconosMatrixException::selfThrow("ioMatrix::read not yet implemented for block matrix.");

  infile.precision(15);
  infile.setf(std::ios::scientific);
  DenseMat * p = m.dense();

  // Dim of the matrix are given in the first line.
  // Just use to check that sizes are consistents.

  unsigned int s1, s2;
  infile >> s1;
  infile >> s2;

  if (s1 != p->size1() || s2 != p->size2())
    p->resize(s1, s2);

  // Note: using istream stl iterator seems to be 2-times faster than << with a loop over matrix data.
  //  copy((std::istream_iterator<double>(infile)), std::istream_iterator<double>(), (p->data()).begin());
  // But it fails with column-major saving ... (ok if user write its matrix in a column-major way)

  DenseMat::iterator1 it;
  DenseMat::iterator2 it2;
  //   std::cout.precision(15);
  //   std::cout.setf(std::ios::scientific);

  for (it = p->begin1(); it != p->end1(); ++it)
  {


    //std::string sz;
    //double  tmp;
    //getline(infile,sz);
    //std::cout << sz;
    for (it2 = it.begin(); it2 != it.end(); ++it2)
    {
      //infile >> tmp;
      //std ::cout << tmp << " ";
      //(*it2)=tmp;
      //std ::cout << *it2 << " ";
      infile >> *it2;
      //        assert(!infile.eof());
      //        assert(!infile.fail());
      //        assert(!infile.bad());
      assert(infile.good());
    }
    // std::cout  << std:: endl;
  }

  // Old version: result in Boost format for ouptut
  //  DenseMat * p = m.dense();
  //  infile >> *p;

  infile.close();
  return true;
}

template<> bool ioMatrix::write(const SiconosMatrix& m, const std::string& outputType) const
{
  // Open file and various checks
  std::ofstream outfile;
  if (Mode == "ascii")
    outfile.open(FileName.c_str(), std::ofstream::out);
  else if (Mode == "binary")
    outfile.open(FileName.c_str(), std::ofstream::binary);
  else
    SiconosMatrixException::selfThrow("ioMatrix::write Incorrect mode for writing");

  if (!outfile.good())
    SiconosMatrixException::selfThrow("ioMatrix:: write error : Fail to open \"" + FileName + "\"");

  if (m.isBlock())
    SiconosMatrixException::selfThrow("ioMatrix:: write error : not yet implemented for BlockMatrix");

  outfile.precision(15);
  outfile.setf(std::ios::scientific);
  // Writing

  if (outputType != "noDim")
    outfile << m.size(0) << " " << m.size(1) << std::endl;

  if (m.getNum() == 1)
  {
    DenseMat * p = m.dense();
    DenseMat::iterator1 row;
    DenseMat::iterator2 col;
    double tmp;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      //std::copy(row.begin(),row.end(),std::ostream_iterator<double>(outfile," "));
      for (col = row.begin(); col != row.end() ; ++col)
      {
        tmp = *col;
        if (fabs(tmp) < std::numeric_limits<double>::min()) tmp = 0.0;
        outfile << tmp << " " ;
        assert(outfile.good());
      }
      outfile << std::endl;

    }
  }
  else if (m.getNum() == 2)
  {
    TriangMat * p = m.triang();
    TriangMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else if (m.getNum() == 3)
  {
    SymMat * p = m.sym();
    SymMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else if (m.getNum() == 4)
  {
    SparseMat * p = m.sparse();
    SparseMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else
  {
    BandedMat * p = m.banded();
    BandedMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }

  outfile.close();
  return true;
}


// To be used later ... ?
//   template <class T, class iterator1> friend void write( const T& obj, iterator1 row, std::ofstream outfile)
//   {
//       for(row = obj.begin1(); row!=obj.end1() ; ++row)
//  {
//    std::copy(row.begin(),row.end(),std::ostream_iterator<double>(outfile," "));
//    outfile << std::endl;
//  }
//     }
