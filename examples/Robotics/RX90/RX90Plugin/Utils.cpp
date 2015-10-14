//
// Copyright (C) INRIA 1999-2008
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//%
// @file kernel/Utils.cpp
// @author Rémy MOZUL
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): mozul@inria.fr
//

#include "Utils.hpp"

/*********************************************************
* Matrix reduction           *
*********************************************************/

void reduceMatrix(matrix<double, column_major> & mat, vector<int> & index1, vector<int> & index2)
{
  if (index1.size() > 0)
  {
    if (mat.size1() < index1.size())
      std::cout << "Wrong size...\n";
    else
    {
      for (unsigned int i = 0; i < index1.size(); i++)
        row(mat, i) = row(mat, index1(i));
      mat.resize(index1.size(), mat.size2(), true);
    }
  }

  if (index2.size() > 0)
  {
    if (mat.size2() < index2.size())
      std::cout << "Wrong size...\n";
    else
    {
      for (unsigned int i = 0; i < index2.size(); i++)
        column(mat, i) = column(mat, index2(i));
      mat.resize(mat.size1(), index2.size(), true);
    }
  }
}

void reduceMatrix(matrix<double, column_major> & mat, vector<bool> & index1, vector<int> & index2)
{
  if (index1.size() > 0)
  {
    if (mat.size1() != index1.size())
      std::cout << "Wrong size...\n";
    else
    {
      matrix<double> tmp(index1.size(), mat.size2());
      int k = 0;
      for (unsigned int i = 0; i < index1.size(); i++)
        if (index1(i))
        {
          row(mat, k) = row(mat, i);
          k++;
        }
      mat.resize(k, mat.size2(), true);
    }
  }

  if (index2.size() > 0)
  {
    if (mat.size2() < index2.size())
      std::cout << "Wrong size...\n";
    else
    {
      for (unsigned int i = 0; i < index2.size(); i++)
        column(mat, i) = column(mat, index2(i));
      mat.resize(mat.size1(), index2.size(), true);
    }
  }
}

/*********************************************************
* Vector reduction           *
*********************************************************/

vector<double> reduceVector(vector<double> & vec, vector<int> & index)
{
  vector<double> ret(0);

  if (vec.size() < index.size())
  {
    std::cout << "Wrong size...\n";
    return ret;
  }

  if (index.size() > 0)
  {
    ret.resize(index.size());
    for (unsigned int i = 0; i < index.size(); i++)
      ret[i] = vec[index[i]];
    return ret;
  }

  ret.resize(vec.size());
  ret = vec;
  return ret;
}


vector<double> reduceVector(vector<double, array_adaptor<double> > & vec, vector<int> & index)
{
  vector<double> ret(0);

  if (vec.size() < index.size())
  {
    std::cout << "Wrong size...\n";
    return ret;
  }

  if (index.size() > 0)
  {
    ret.resize(index.size());
    for (unsigned int i = 0; i < index.size(); i++)
      ret[i] = vec[index[i]];
    return ret;
  }

  ret.resize(vec.size());
  ret = vec;
  return ret;
}

void reduceVector(vector<double> & vec, vector<bool> & index)
{
  if (vec.size() != index.size())
    std::cout << "Wrong size...\n";
  else
  {
    unsigned int k = 0;
    for (unsigned int i = 0; i < index.size(); i++)
    {
      if (index(i))
      {
        vec[k] = vec[i];
        k++;
      }
    }
    vec.resize(k);
  }
}

void reduceVector(vector<bool> & vec, vector<bool> & index)
{
  if (vec.size() != index.size())
    std::cout << "Wrong size...\n";
  else
  {
    unsigned int k = 0;
    for (unsigned int i = 0; i < index.size(); i++)
    {
      if (index(i))
      {
        vec[k] = vec[i];
        k++;
      }
    }
    vec.resize(k);
  }
}

/*********************************************************
* Vector expansion           *
*********************************************************/

void expandVector(double * tab, vector<double> & vector, vector<int> & index)
{
  //  if(vector2.size()>vector1.size() || vector2.size()!=index.size())
  //    std::cout << "Wrong size...\n";
  //  else
  if (index.size() > 0)
    for (unsigned int i = 0; i < vector.size(); i++)
      tab[index[i]] = vector[i];
  else
    for (unsigned int i = 0; i < vector.size(); i++)
      tab[i] = vector[i];
}

void expandVector(double * tab, vector<double> & vector, vector<bool> & index)
{
  //  if(vector.size()>index.size() && index.size()>0)
  //    std::cout << "Wrong size...\n";
  //  else
  if (index.size() > 0)
  {
    unsigned int k = 0;
    for (unsigned int i = 0; i < index.size(); i++)
      if (index[i])
      {
        tab[i] = vector[k];
        k++;
      }
  }
  else
  {
    for (unsigned int i = 0; i < vector.size(); i++)
      tab[i] = vector[i];
  }
}

