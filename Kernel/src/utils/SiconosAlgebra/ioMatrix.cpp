#include "ioMatrix.h"
#include "SimpleMatrix.h"
#include "SiconosMatrixException.h"
#include <fstream>
#include <iterator>

// Default (private)
ioMatrix::ioMatrix(): ioObject() {}

ioMatrix::ioMatrix(const std::string& file, const std::string& mode): ioObject(file, mode) {}

const bool ioMatrix::read(SiconosMatrix& m) const
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
  DenseMat * p = m.getDensePtr();

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

  for (it = p->begin1(); it != p->end1(); ++it)
  {
    for (it2 = it.begin(); it2 != it.end(); ++it2)
      infile >> (*it2);
  }

  // Old version: result in Boost format for ouptut
  //  DenseMat * p = m.getDensePtr();
  //  infile >> *p;

  infile.close();
  return true;
}

const bool ioMatrix::write(const SiconosMatrix& m, const std::string& outputType) const
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

  // Writing

  if (outputType != "noDim")
    outfile << m.size(0) << " " << m.size(1) << std::endl;

  if (m.getNum() == 1)
  {
    DenseMat * p = m.getDensePtr();
    DenseMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else if (m.getNum() == 2)
  {
    TriangMat * p = m.getTriangPtr();
    TriangMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else if (m.getNum() == 3)
  {
    SymMat * p = m.getSymPtr();
    SymMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else if (m.getNum() == 4)
  {
    SparseMat * p = m.getSparsePtr();
    SparseMat::iterator1 row;
    for (row = p->begin1(); row != p->end1() ; ++row)
    {
      std::copy(row.begin(), row.end(), std::ostream_iterator<double>(outfile, " "));
      outfile << std::endl;
    }
  }
  else
  {
    BandedMat * p = m.getBandedPtr();
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
