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
#include "LinearTIDS.h"
using namespace std;

// --- Constructors ---

// Default constructor
LinearTIDS::LinearTIDS():
  LinearDS()
{
  DSType = LITIDS;
}

// From xml file (newNsds is optional)
LinearTIDS::LinearTIDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  LinearDS(dsXML, newNsds)
{
  if (dsXML != NULL)
  {
    DSType = LITIDS;

    // Everything is done during the call to LinearDS constructor. In the present one, we just reject cases
    // where A or b are given as plug-in rather than as matrices.

    // pointer to xml
    LinearDSXML * ldsxml = (static_cast <LinearDSXML*>(dsxml));

    // reject case where A or b is a plug-in
    if (ldsxml->isAPlugin())
      RuntimeException::selfThrow("LinearTIDS - xml constructor, A is given as a plug-in, which should not be.");
    if (ldsxml->isBPlugin())
      RuntimeException::selfThrow("LinearTIDS - xml constructor, b is given as a plug-in, which should not be.");
  }
  else
    RuntimeException::selfThrow("LinearTIDS - xml constructor, xml file = NULL");

  isPlugin["A"] = false;
  isPlugin["b"] = false;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a plugin
LinearTIDS::LinearTIDS(const int newNumber, const SiconosVector& newX0, const SiconosMatrix& newA):
  LinearDS(newNumber, newX0, newA)
{
  DSType = LITIDS;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data, A from a plugin
LinearTIDS::LinearTIDS(const int newNumber, const SiconosVector& newX0,
                       const SiconosMatrix& newA, const SiconosVector& newB):
  LinearDS(newNumber, newX0, newA, newB)
{
  DSType = LITIDS;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// Copy constructor
LinearTIDS::LinearTIDS(const LinearTIDS & lds):
  LinearDS(lds)
{
  DSType = LITIDS;
  isPlugin["A"] = false;
  isPlugin["b"] = false;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

LinearTIDS::LinearTIDS(const DynamicalSystem & newDS):
  LinearDS(newDS)
{
  DSType = LITIDS;
  isPlugin["A"] = false;
  isPlugin["b"] = false;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

LinearTIDS::~LinearTIDS()
{}

// void LinearTIDS::initialize(const double time, const unsigned int sizeOfMemory)
// {
//   // reset x to x0, xFree and r to zero.
//   *x = *x0;
//   xFree->zero(); r->zero();

//   // Initialize memory vectors
//   initMemory(sizeOfMemory);

// }

bool LinearTIDS::checkDynamicalSystem()
{
  bool output = true;
  // n
  if (n == 0)
  {
    RuntimeException::selfThrow("LinearTIDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }
  // x0 != NULL
  if (x0 == NULL)
  {
    RuntimeException::selfThrow("LinearTIDS::checkDynamicalSystem - x0 not set.");
    output = false;
  }

  // A
  if (A == NULL)
  {
    RuntimeException::selfThrow("LinearTIDS::checkDynamicalSystem - A not set.");
    output = false;
  }

  return output;
}

void LinearTIDS::computeRhs(const double time, const bool)
{
  // compute right-hand side
  *rhs = *A * *x;

  // add b if required
  if (b != NULL)
    *rhs += *b;

  // compute and add Tu if required
  if (u != NULL)
  {
    if (isPlugin["u"]) // if u is a plug-in function
      computeU(time);
    if (T != NULL)
    {
      if (isPlugin["T"]) // if T is a plug-in function
        computeT();
      *rhs += *T ** u;
    }
    else
      *rhs += * u;
  }
}

void LinearTIDS::computeJacobianXRhs(const double time, const bool)
{}

void LinearTIDS::display() const
{
  DynamicalSystem::display();
  cout << "=== Linear Time-invariant system display ===" << endl;
  cout << "- A " << endl;
  if (A != NULL) A->display();
  else cout << "-> NULL" << endl;
  cout << "- b " << endl;
  if (b != NULL) b->display();
  else cout << "-> NULL" << endl;
  cout << "============================================" << endl;
}

LinearTIDS* LinearTIDS::convert(DynamicalSystem* ds)
{
  LinearTIDS* lsds = dynamic_cast<LinearTIDS*>(ds);
  return lsds;
}
