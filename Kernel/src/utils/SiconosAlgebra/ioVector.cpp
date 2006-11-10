#include "ioVector.h"
#include "SimpleVector.h"
#include "SiconosVectorException.h"
#include <boost/numeric/ublas/io.hpp>
#include<fstream>
//SimpleVector *tmp = new SimpleVector();

SimpleVector* ioVector::temporary = new SimpleVector(1);
bool ioVector::writeSimpleBinary = false;

// Default private
ioVector::ioVector(): ioObject() {}

ioVector::ioVector(const std::string& file, const std::string& mode): ioObject(file, mode) {}

ioVector::~ioVector(void)
{
  if (writeSimpleBinary)
    delete(temporary);
}

const bool ioVector::read(SiconosVector& m) const
{

  if (Mode == "ascii")
  {

    std::ifstream infile(FileName.c_str(), std::ifstream::in);
    if (infile == NULL)
    {
      SiconosVectorException::selfThrow(" ioVector::read : Fail to open file \"" + FileName + "\"");
    }

    int col;
    infile >> col;

    if (m.isBlock())
      SiconosVectorException::selfThrow(" ioVector::read : read BlockVector is not implemented");
    DenseVect p(col);
    infile >> p;

    SimpleVector tmp(p);
    m = tmp;
    infile.close();
    return true;
  }
  else if (Mode == "binary")
  {

    FILE *infile = fopen(FileName.c_str(), "rb");

    if (infile == NULL)
    {
      SiconosVectorException::selfThrow(" ioVector::read : Fail to open file \"" + FileName + "\"");
    }

    int alpha;
    int err;
    fread((char*)&alpha, sizeof(int), 1, infile);

    if (m.isBlock())
      SiconosVectorException::selfThrow(" ioVector::read : read BlockVector is not implemented");

    DenseVect p(alpha);
    err = fread(&p, sizeof(DenseVect(alpha)), 1, infile);
    SimpleVector tmp(p);
    m = tmp;
    fclose(infile);
    return true;
  }
  else
  {
    SiconosVectorException::selfThrow(" ioVector::read : incorrect mode for reading vector");
    return false;
  }
}

const bool ioVector::write(const SiconosVector& m, const std::string&) const
{
  if (Mode == "ascii")
  {
    std::ofstream outfile(FileName.c_str());

    if (!outfile.is_open())
    {
      SiconosVectorException::selfThrow(" ioVector::write : Fail to open file \"" + FileName + "\"");
    }
    if (m.isBlock())
      SiconosVectorException::selfThrow(" ioVector::write : write BlockVector is not implemented");

    int col;
    col = m.size();
    outfile << col;
    outfile << '\n';

    DenseVect p;
    if (m.getNum() == 1)
      p = m.getDense();
    else if (m.getNum() == 2)
      p = m.getSparse();

    outfile << p;
    outfile.close();
    return true;
  }
  else if (Mode == "binary")
  {
    FILE *outfile = fopen(FileName.c_str(), "wb");
    if (outfile == NULL)
    {
      SiconosVectorException::selfThrow(" ioVector::write : Fail to open file \"" + FileName + "\"");
    }
    if (m.isBlock())
      SiconosVectorException::selfThrow(" ioVector::write : write BlockVector is not implemented");

    writeSimpleBinary = true;

    int col;
    col = m.size();

    temporary->resize(col, false);
    *temporary = m;
    fwrite((char*)&col, sizeof(int), 1, outfile);
    fwrite(temporary->getDensePtr(), sizeof(DenseVect(col)), 1, outfile);
    fclose(outfile);
    return true;
  }
  else
  {
    SiconosVectorException::selfThrow(" ioVector::write : incorrect mode for writing vector");
    return false;
  }
}

