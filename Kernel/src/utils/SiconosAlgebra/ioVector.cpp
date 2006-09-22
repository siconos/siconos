#include <fstream>
#include <iostream>
#include <stdio.h>

#include "ioVector.h"

//MySimpleVector *tmp = new MySimpleVector(DENSE);

MySimpleVector* ioVector::temporary = new MySimpleVector(DENSE);
bool ioVector::writeSimpleBinary = false;

// Default private
ioVector::ioVector(void)
{
  FileName = "";
  Mode = "ascii";
}

// Default public
ioVector::ioVector(const std::string& file, const std::string& mode)
{
  FileName = file;
  Mode = mode;
}

ioVector::~ioVector(void)
{
  if (writeSimpleBinary)
    delete(temporary);
}

bool ioVector::read(MySiconosVector& m)const
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

    if (m.isBlock() == false)
    {
      DenseVect p(col);
      infile >> p;

      MySimpleVector tmp(p);
      m = tmp;
    }
    else
    {
      SiconosVectorException::selfThrow(" ioVector::read : read BlockVector is not implemented");
    }
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

    if (m.isBlock() == false)
    {
      DenseVect p(alpha);
      err = fread(&p, sizeof(DenseVect(alpha)), 1, infile);
      MySimpleVector tmp(p);
      m = tmp;
    }
    else
    {
      SiconosVectorException::selfThrow(" ioVector::read : read BlockVector is not implemented");
    }
    fclose(infile);
    return true;
  }
  else
  {
    SiconosVectorException::selfThrow(" ioVector::read : incorrect mode for reading vector");

  }
}

bool ioVector::write(const MySiconosVector& m)
{

  if (Mode == "ascii")
  {
    std::ofstream outfile(FileName.c_str());

    if (!outfile.is_open())
    {
      SiconosVectorException::selfThrow(" ioVector::write : Fail to open file \"" + FileName + "\"");
    }
    if (m.isBlock() == false)
    {
      int col;
      col = m.size();
      outfile << col;
      outfile << '\n';

      DenseVect p;
      if (m.GetNum() == 1)
        p = m.GetDense();
      else if (m.GetNum() == 2)
        p = m.GetSparse();

      outfile << p;
    }
    else
    {
      SiconosVectorException::selfThrow(" ioVector::write : write BlockVector is not implemented");
    }
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
    if (m.isBlock() == false)
    {
      writeSimpleBinary = true;

      int col;
      col = m.size();

      temporary->resize(col, false);
      *temporary = m;
      fwrite((char*)&col, sizeof(int), 1, outfile);
      fwrite(temporary->GetDensePtr(), sizeof(DenseVect(col)), 1, outfile);

    }
    else
    {
      SiconosVectorException::selfThrow(" ioVector::write : write BlockVector is not implemented");
    }
    fclose(outfile);
    return true;
  }
  else
  {
    SiconosVectorException::selfThrow(" ioVector::write : incorrect mode for writing vector");

  }
}

bool ioVector::rawWrite(const MySiconosVector& m)
{

  if (Mode == "ascii")
  {
    std::ofstream outfile(FileName.c_str());

    if (!outfile.is_open())
    {
      SiconosVectorException::selfThrow(" ioVector::rawWrite : Fail to open file \"" + FileName + "\"");
    }
    if (m.isBlock() == false)
    {
      DenseVect p;
      if (m.GetNum() == 1)
        p = m.GetDense();
      else if (m.GetNum() == 2)
        p = m.GetSparse();

      outfile << p;
    }
    else
    {
      SiconosVectorException::selfThrow(" ioVector::rawWrite : rawWrite BlockVector is not implemented");
    }
    outfile.close();
    return true;
  }
  else if (Mode == "binary")
  {
    FILE *outfile = fopen(FileName.c_str(), "wb");
    if (outfile == NULL)
    {
      SiconosVectorException::selfThrow(" ioVector::rawWrite : Fail to open file \"" + FileName + "\"");
    }
    if (m.isBlock() == false)
    {
      int col;
      col = m.size();
      writeSimpleBinary = true;

      temporary->resize(col, false);
      *temporary = m;
      fwrite(temporary->GetDensePtr(), sizeof(DenseVect(col)), 1, outfile);

    }
    else
    {
      SiconosVectorException::selfThrow(" ioVector::rawWrite : rawWrite BlockVector is not implemented");
    }
    fclose(outfile);
    return true;
  }
  else
  {
    SiconosVectorException::selfThrow(" ioVector::rawWrite : incorrect mode for writing vector");
  }
}

bool ioVector::read(MySiconosMatrix& m)const
{
  SiconosVectorException::selfThrow(" ioVector::read(MySiconosMatrix&) is forbidden");
}
bool ioVector::write(const MySiconosMatrix& m)
{
  SiconosVectorException::selfThrow(" ioVector::write(MySiconosMatrix&) is forbidden");
}
bool ioVector::rawWrite(const MySiconosMatrix& m)
{
  SiconosVectorException::selfThrow(" ioVector::rawWrite(MySiconosMatrix&) is forbidden");
}

