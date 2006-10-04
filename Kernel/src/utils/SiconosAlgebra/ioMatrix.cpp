#include <fstream>
#include <iostream>
#include <stdio.h>

#include "ioMatrix.h"


MySimpleMatrix* ioMatrix::temporary = new MySimpleMatrix(1, 1);
bool ioMatrix::writeSimpleBinary = false;

// Default private
ioMatrix::ioMatrix(void)
{
  FileName = "";
  Mode = "ascii";
}

// Default public
ioMatrix::ioMatrix(const std::string& file, const std::string& mode)
{
  FileName = file;
  Mode = mode;
}

ioMatrix::~ioMatrix(void)
{
  if (writeSimpleBinary)
    delete(temporary);
}

bool ioMatrix::read(MySiconosMatrix& m)const
{

  if (Mode == "ascii")
  {

    std::ifstream infile(FileName.c_str(), std::ifstream::in);
    if (infile == NULL)
      SiconosMatrixException::selfThrow("function read error : Fail to open \"" + FileName + "\"");

    int row, col;
    infile >> row;
    infile >> col;

    if (m.isBlock() == false)
    {
      DenseMat p(row, col);
      infile >> p;
      MySimpleMatrix tmp(p);
      m = tmp;
    }
    else
    {
      int nbRow = row;
      int nbCol = col;
      DenseMat p ;
      mapped Mmap(nbRow, nbCol);
      for (int i = 0; i < nbRow; i++)
      {
        for (int j = 0; j < nbCol; j++)
        {
          infile >> row;
          infile >> col;
          p.resize(row, col, false);
          infile >> p;
          Mmap(i, j) = new MySimpleMatrix(p);
        }
      }
      MyBlockMatrix tmp(Mmap);
      m = tmp;

      for (int i = 0; i < nbRow; i++)
      {
        for (int j = 0; j < nbCol; j++)
        {
          delete(Mmap(i, j));
        }
      }
    }
    infile.close();
    return true;
  }
  else if (Mode == "binary")
  {

    FILE *infile = fopen(FileName.c_str(), "rb");

    if (infile == NULL)
    {
      SiconosMatrixException::selfThrow("function read error : Fail to open \"" + FileName + "\"");
    }

    int alpha;
    int beta;
    int err;
    fread((char*)&alpha, sizeof(int), 1, infile);
    fread((char*)&beta, sizeof(int), 1, infile);

    if (m.isBlock() == false)
    {
      DenseMat p(alpha, beta);
      err = fread(&p, sizeof(DenseMat(alpha, beta)), 1, infile);
      MySimpleMatrix tmp(p);
      m = tmp;
    }
    else
    {

      mapped mytest(alpha, beta);
      DenseMat p;

      int row, col;
      for (int i = 0; i < alpha; i++)
      {
        for (int j = 0; j < beta; j++)
        {
          fread((char*)&row, sizeof(int), 1, infile);
          fread((char*)&col, sizeof(int), 1, infile);
          p.resize(row, col, false);
          err = fread((char*)&p, sizeof(DenseMat(row, col)), 1, infile);
          mytest(i, j) = new MySimpleMatrix(p);
        }
      }
      MyBlockMatrix tmp(mytest);
      m = tmp;
      for (int i = 0; i < alpha; i++)
      {
        for (int j = 0; j < beta; j++)
        {
          delete mytest(i, j);
        }
      }
    }
    fclose(infile);
    return true;
  }
  else
  {
    SiconosMatrixException::selfThrow("Incorrect mode for reading");
    return false;
  }
}

bool ioMatrix::write(const MySiconosMatrix& m)
{

  if (Mode == "ascii")
  {
    std::ofstream outfile(FileName.c_str());

    if (!outfile.is_open())
    {
      SiconosMatrixException::selfThrow("function write error : Fail to open \"" + FileName + "\"");

    }
    if (m.isBlock() == false)
    {
      int row, col;
      row = m.size1();
      col = m.size2();
      outfile << row << ' ' << col;
      outfile << '\n';

      DenseMat p;
      if (m.getNum() == 1)
        p = m.getDense();
      else if (m.getNum() == 2)
        p = m.getTriang();
      else if (m.getNum() == 3)
        p = m.getSym();

      outfile << p;
    }
    else
    {
      mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).getMap();
      int col = Mmap.size2();
      int row = Mmap.size1();
      outfile << row << ' ' << col;
      outfile << '\n';

      int k, l;
      DenseMat p;
      if (STDMAP == 1)
      {
        mapped::array_type::iterator it;
        for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
        {
          if ((it->second)->getNum() == 1)
            p = (it->second)->getDense();
          else if ((it->second)->getNum() == 2)
            p = (it->second)->getTriang();
          else if ((it->second)->getNum() == 3)
            p = (it->second)->getSym();

          k = p.size1();
          l = p.size2();
          outfile << k << ' ';
          outfile << l << ' ';
          outfile << p;
          outfile << '\n';
        }
      }
      else
      {
        mapped::iterator1 it;
        mapped::iterator2 it2;

        for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
        {
          for (it2 = it.begin(); it2 != it.end(); it2 ++)
          {
            if ((**it2).getNum() == 1)
              p = (**it2).getDense();
            else if ((**it2).getNum() == 2)
              p = (**it2).getTriang();
            else if ((**it2).getNum() == 3)
              p = (**it2).getSym();

            k = p.size1();
            l = p.size2();
            outfile << k << ' ';
            outfile << l << ' ';
            outfile << p;
            outfile << '\n';
          }
        }
      }
    }
    outfile.close();
    return true;
  }
  else if (Mode == "binary")
  {
    FILE *outfile = fopen(FileName.c_str(), "wb");
    if (outfile == NULL)
    {
      SiconosMatrixException::selfThrow("function write error : Fail to open \"" + FileName + "\"");
    }
    if (m.isBlock() == false)
    {
      writeSimpleBinary = true;

      int row, col;
      row = m.size1();
      col = m.size2();

      temporary->resize(row, col, false);
      *temporary = m;
      fwrite((char*)&row, sizeof(int), 1, outfile);
      fwrite((char*)&col, sizeof(int), 1, outfile);
      fwrite(temporary->getDensePtr(), sizeof(DenseMat(row, col)), 1, outfile);

    }
    else
    {
      mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).getMap();
      int col = Mmap.size2();
      int row = Mmap.size1();
      fwrite((char*)&row, sizeof(int), 1, outfile);
      fwrite((char*)&col, sizeof(int), 1, outfile);
      int k, l, err;
      if (STDMAP == 1)
      {
        mapped::array_type::iterator it;
        for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
        {
          k = (it->second)->size1();
          l = (it->second)->size2();
          fwrite((char*)&k, sizeof(int), 1, outfile);
          fwrite((char*)&l, sizeof(int), 1, outfile);


          if ((it->second)->getNum() == 1)
          {
            const DenseMat *p = (it->second)->getDensePtr();
            err = fwrite(p, sizeof(DenseMat(k, l)), 1, outfile);
          }
          else if ((it->second)->getNum() == 2)
          {
            const TriangMat *p = (it->second)->getTriangPtr();
            err = fwrite(p, sizeof(TriangMat(k, l)), 1, outfile);
          }
          else if ((it->second)->getNum() == 3)
          {
            const SymMat *p = (it->second)->getSymPtr();
            err = fwrite(p, sizeof(SymMat(k, l)), 1, outfile);
          }
        }
      }
      else
      {

        mapped::iterator1 it;
        mapped::iterator2 it2;
        for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
        {
          for (it2 = it.begin(); it2 != it.end(); it2 ++)
          {
            k = (**it2).size1();
            l = (**it2).size2();
            fwrite((char*)&k, sizeof(int), 1, outfile);
            fwrite((char*)&l, sizeof(int), 1, outfile);
            err = fwrite((*it2), sizeof(MySimpleMatrix(k, l)), 1, outfile);
          }
        }
      }
    }
    fclose(outfile);
    return true;
  }
  else
  {
    SiconosMatrixException::selfThrow("Incorrect mode for writing");
    return false;
  }
}

bool ioMatrix::rawWrite(const MySiconosMatrix& m)
{

  if (Mode == "ascii")
  {
    std::ofstream outfile(FileName.c_str());

    if (!outfile.is_open())
    {
      SiconosMatrixException::selfThrow("function rawWrite error : Fail to open \"" + FileName + "\"");
    }
    if (m.isBlock() == false)
    {
      //int row, col;
      //row = m.size1 ();
      //col = m.size2 ();
      //outfile << row << ' '<<col;
      //outfile << '\n';
      DenseMat p;
      if (m.getNum() == 1)
        p = m.getDense();
      else if (m.getNum() == 2)
        p = m.getTriang();
      else if (m.getNum() == 3)
        p = m.getSym();

      outfile << p;
    }
    else
    {
      mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).getMap();
      //int col = Mmap.size2 ();
      //int row = Mmap.size1 ();
      //outfile << row << ' '<<col;
      //outfile << '\n';

      //int k, l;
      DenseMat p;
      if (STDMAP == 1)
      {
        mapped::array_type::iterator it;
        for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
        {
          if ((it->second)->getNum() == 1)
            p = (it->second)->getDense();
          else if ((it->second)->getNum() == 2)
            p = (it->second)->getTriang();
          else if ((it->second)->getNum() == 3)
            p = (it->second)->getSym();

          //k = p.size1 ();
          //l = p.size2 ();
          //outfile << k <<' ';
          //outfile << l <<' ';
          //outfile << p;
          outfile << '\n';
        }
      }
      else
      {
        mapped::iterator1 it;
        mapped::iterator2 it2;

        for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
        {
          for (it2 = it.begin(); it2 != it.end(); it2 ++)
          {
            if ((**it2).getNum() == 1)
              p = (**it2).getDense();
            else if ((**it2).getNum() == 2)
              p = (**it2).getTriang();
            else if ((**it2).getNum() == 3)
              p = (**it2).getSym();

            //k = p.size1 ();
            //l = p.size2 ();
            //outfile << ' '<<k<<' ';
            //outfile << l <<' ';
            outfile << p;
            outfile << '\n';
          }
        }
      }
    }
    outfile.close();
    return true;
  }
  else if (Mode == "binary")
  {
    FILE *outfile = fopen(FileName.c_str(), "wb");
    if (outfile == NULL)
    {
      SiconosMatrixException::selfThrow("function rawWrite error : Fail to open \"" + FileName + "\"");
    }
    if (m.isBlock() == false)
    {
      int row, col;
      row = m.size1();
      col = m.size2();
      writeSimpleBinary = true;

      temporary->resize(row, col, false);
      *temporary = m;
      //fwrite((char*)&row, sizeof(int), 1, outfile);
      //fwrite((char*)&col, sizeof(int), 1, outfile);
      fwrite(temporary->getDensePtr(), sizeof(DenseMat(row, col)), 1, outfile);

    }
    else
    {
      mapped Mmap = (dynamic_cast<const MyBlockMatrix&>(m)).getMap();
      //int col = Mmap.size2 ();
      //int row = Mmap.size1 ();
      //fwrite((char*)&row, sizeof(int), 1, outfile);
      //fwrite((char*)&col, sizeof(int), 1, outfile);
      int k, l, err;
      if (STDMAP == 1)
      {
        mapped::array_type::iterator it;
        for (it = (Mmap.data()).begin(); it != (Mmap.data()).end(); ++it)
        {
          k = (it->second)->size1();
          l = (it->second)->size2();
          //fwrite((char*)&k, sizeof(int), 1, outfile);
          //fwrite((char*)&l, sizeof(int), 1, outfile);


          if ((it->second)->getNum() == 1)
          {
            const DenseMat *p = (it->second)->getDensePtr();
            err = fwrite(p, sizeof(DenseMat(k, l)), 1, outfile);
          }
          else if ((it->second)->getNum() == 2)
          {
            const TriangMat *p = (it->second)->getTriangPtr();
            err = fwrite(p, sizeof(TriangMat(k, l)), 1, outfile);
          }
          else if ((it->second)->getNum() == 3)
          {
            const SymMat *p = (it->second)->getSymPtr();
            err = fwrite(p, sizeof(SymMat(k, l)), 1, outfile);
          }
        }
      }
      else
      {

        mapped::iterator1 it;
        mapped::iterator2 it2;
        for (it = Mmap.begin1(); it != Mmap.end1(); ++it)
        {
          for (it2 = it.begin(); it2 != it.end(); it2 ++)
          {
            k = (**it2).size1();
            l = (**it2).size2();
            //fwrite((char*)&k, sizeof(int), 1, outfile);
            //fwrite((char*)&l, sizeof(int), 1, outfile);
            err = fwrite((*it2), sizeof(MySimpleMatrix(k, l)), 1, outfile);
          }
        }
      }
    }
    fclose(outfile);
    return true;
  }
  else
  {
    SiconosMatrixException::selfThrow("Incorrect mode for writing");
    return false;
  }
}

bool ioMatrix::read(MySiconosVector& m)const
{
  SiconosMatrixException::selfThrow("ioMatrix::read(MySiconosVector&) is forbidden");
  return false;
}

bool ioMatrix::write(const MySiconosVector& m)
{
  SiconosMatrixException::selfThrow("ioMatrix::write(const MySiconosVector&) is forbidden");
  return false;
}

bool ioMatrix::rawWrite(const MySiconosVector& m)
{
  SiconosMatrixException::selfThrow("ioMatrix::rawWrite(const MySiconosVector&) is forbidden");
  return false;
}

