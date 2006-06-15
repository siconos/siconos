/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include "SiconosMatrix.h"
using namespace std;

// LAPACK F77 function not declared in lapackpp.h
Extern  "C"
{
  void F77NAME(dgetri)(integer * m, doublereal * A, integer * lda, integer * ipiv,
  doublereal * WORK, integer * LWORK, integer * info);
}

// Private functions (why?)

void SiconosMatrix::verbose(const string& msg)
{
  if (printVerbose)
    cout << msg << endl;
}

// --- CONSTRUCTORS ---

// Default (private)
SiconosMatrix::SiconosMatrix(const bool& isblock): isBlockMatrix(isblock)
{}

// --- DESTRUCTOR ---

SiconosMatrix::~SiconosMatrix()
{}


/****************** output stream ******************/
ostream& operator << (ostream &ostrm, SiconosMatrix& m)
{
  if (m.size(0) == 0)
  {
    cout << "Display SiconosMatrix : empty matrix" << endl;
  }
  //SiconosMatrixException::selfThrow("Try to display a 0-size matrix");

  else cout << m.getLaGenMatDouble() ;
  return ostrm;
}

/****************** input stream ******************/
istream& operator >> (istream &istrm, SiconosMatrix& m)
{
  for (unsigned int i = 0; i < m.size(0); i++)
    for (unsigned int j = 0; j < m.size(1); j++)
    {
      cout << '[' << i + 1 << ',' << j + 1 << "] = ";
      cin >> m(i, j);
      if (cin.fail())
        cout << "INPUT ERROR";
    }
  return istrm;
}


/*************************************************/
SiconosMatrix& SiconosMatrix::operator = (const SiconosMatrix& m)
{
  cout << "Warning, SiconosMatrix virtual class operator = call. This is strange ..." << endl;
  return *this;
}



