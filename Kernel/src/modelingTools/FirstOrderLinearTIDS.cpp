/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "FirstOrderLinearTIDS.h"
#include "FirstOrderLinearDSXML.h"

using namespace std;

// --- Constructors ---

// Default constructor
FirstOrderLinearTIDS::FirstOrderLinearTIDS():
  FirstOrderLinearDS()
{
  DSType = FOLTIDS;
}

// From xml file (newNsds is optional)
FirstOrderLinearTIDS::FirstOrderLinearTIDS(DynamicalSystemXMLSPtr dsXML, SP::NonSmoothDynamicalSystem newNsds):
  FirstOrderLinearDS(dsXML, newNsds)
{
  // Everything is done during the call to FirstOrderLinearDS constructor. In the present one, we just reject cases
  // where A or b are given as plug-in rather than as matrices.

  // pointer to xml

  FirstOrderLinearDSXMLSPtr ldsxml = (boost::static_pointer_cast <FirstOrderLinearDSXML>(dsxml));

  // reject case where A or b is a plug-in
  if (ldsxml->isAPlugin())
    RuntimeException::selfThrow("FirstOrderLinearTIDS - xml constructor, A is given as a plug-in, which should not be.");
  if (ldsxml->isBPlugin())
    RuntimeException::selfThrow("FirstOrderLinearTIDS - xml constructor, b is given as a plug-in, which should not be.");

  isPlugin["A"] = false;
  isPlugin["b"] = false;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a minimum set of data: A
FirstOrderLinearTIDS::FirstOrderLinearTIDS(int newNumber, const SiconosVector& newX0, const SiconosMatrix& newA):
  FirstOrderLinearDS(newNumber, newX0, newA)
{
  DSType = FOLTIDS;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// From a set of data: A and B
FirstOrderLinearTIDS::FirstOrderLinearTIDS(int newNumber, const SiconosVector& newX0,
    const SiconosMatrix& newA, const SiconosVector& newB):
  FirstOrderLinearDS(newNumber, newX0, newA, newB)
{
  DSType = FOLTIDS;
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

FirstOrderLinearTIDS::~FirstOrderLinearTIDS()
{}

void FirstOrderLinearTIDS::initRhs(double time)
{
  if (M)
  {
    if (! invM)
    {

#ifndef WithSmartPtr
      invM = new SimpleMatrix(*M);
      isAllocatedIn["invM"] = true;
#else
      invM.reset(new SimpleMatrix(*M));
#endif

    }
  }

  computeRhs(time);

  if (! jacobianXRhs)  // if not allocated with a set or anything else
  {
    if (A && ! M)  // if M is not defined, then A = jacobianXRhs, no memory allocation for that one.
      jacobianXRhs = A;
    else if (A && M)
    {

#ifndef WithSmartPtr
      jacobianXRhs = new SimpleMatrix(*A); // Copy A into jacobianXRhs
      isAllocatedIn["jacobianXRhs"] = true;
#else
      jacobianXRhs.reset(new SimpleMatrix(*A)); // Copy A into jacobianXRhs
#endif

      // Solve MjacobianXRhs = A
      invM->PLUForwardBackwardInPlace(*jacobianXRhs);
    }
    // else no allocation, jacobian is equal to 0.
  }
}

void FirstOrderLinearTIDS::computeRhs(const double time, const bool)
{

  *x[1] = * r; // Warning: r update is done in Interactions/Relations

  if (A != NULL)
    prod(*A, *x[0], *x[1], false);

  // compute and add b if required
  if (b != NULL)
    *x[1] += *b;

  if (M != NULL)
    invM->PLUForwardBackwardInPlace(*x[1]);
}

void FirstOrderLinearTIDS::computeJacobianXRhs(const double time, const bool)
{
  // Nothing to be done: jacobianXRhs is constant and computed during initialize. But this function is required to avoid call to base class function.
}

void FirstOrderLinearTIDS::display() const
{
  cout << "===> Linear Time-invariant First Order System display, " << number << ")." << endl;
  cout << "- A " << endl;
  if (A != NULL) A->display();
  else cout << "-> NULL" << endl;
  cout << "- b " << endl;
  if (b != NULL) b->display();
  else cout << "-> NULL" << endl;

  cout << "- M: " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;

  cout << "============================================" << endl;
}

FirstOrderLinearTIDS* FirstOrderLinearTIDS::convert(DynamicalSystem* ds)
{
  FirstOrderLinearTIDS* lsds = dynamic_cast<FirstOrderLinearTIDS*>(ds);
  return lsds;
}
