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
using namespace DS;
// --- Constructors ---

// From xml file
FirstOrderLinearTIDS::FirstOrderLinearTIDS(SP::DynamicalSystemXML dsXML): FirstOrderLinearDS(dsXML)
{
  // Everything is done during the call to FirstOrderLinearDS constructor. In the present one, we just reject cases
  // where A or b are given as plug-in rather than as matrices.

  // pointer to xml

  SP::FirstOrderLinearDSXML ldsxml = (boost::static_pointer_cast <FirstOrderLinearDSXML>(_dsxml));

  // reject case where A or b is a plug-in
  if (ldsxml->isAPlugin())
    RuntimeException::selfThrow("FirstOrderLinearTIDS - xml constructor, A is given as a plug-in, which should not be.");
  if (ldsxml->isBPlugin())
    RuntimeException::selfThrow("FirstOrderLinearTIDS - xml constructor, b is given as a plug-in, which should not be.");

  checkDynamicalSystem();
}

// From a minimum set of data: A
FirstOrderLinearTIDS::FirstOrderLinearTIDS(SP::SiconosVector newX0, SP::SiconosMatrix newA):
  FirstOrderLinearDS(newX0, newA)
{
  _DSType = FOLTIDS;
  checkDynamicalSystem();
}

// From a set of data: A and B
FirstOrderLinearTIDS::FirstOrderLinearTIDS(SP::SiconosVector newX0, SP::SiconosMatrix newA, SP::SiconosVector newB):
  FirstOrderLinearDS(newX0, newA, newB)
{
  _DSType = FOLTIDS;
  checkDynamicalSystem();
}

void FirstOrderLinearTIDS::initRhs(double time)
{
  if (_M && !_invM)
    _invM.reset(new SimpleMatrix(*_M));

  computeRhs(time);

  if (! _jacXRhs)  // if not allocated with a set or anything else
  {
    if (_A && ! _M)  // if M is not defined, then A = _jacXRhs, no memory allocation for that one.
      _jacXRhs = _A;
    else if (_A && _M)
    {
      _jacXRhs.reset(new SimpleMatrix(*_A)); // Copy A into _jacXRhs
      // Solve M_jacXRhs = A
      _invM->PLUForwardBackwardInPlace(*_jacXRhs);
    }
    // else no allocation, jacobian is equal to 0.
  }
}

void FirstOrderLinearTIDS::computeRhs(const double time, const bool)
{

  *_x[1] = * _r; // Warning: r update is done in Interactions/Relations

  if (_A)
    prod(*_A, *_x[0], *_x[1], false);

  // compute and add b if required
  if (_b)
    *_x[1] += *_b;

  if (_M)
    _invM->PLUForwardBackwardInPlace(*_x[1]);
}

void FirstOrderLinearTIDS::computeJacobianXRhs(const double time, const bool)
{
  // Nothing to be done: _jacXRhs is constant and computed during initialize. But this function is required to avoid call to base class function.
}

void FirstOrderLinearTIDS::display() const
{
  cout << "===> Linear Time-invariant First Order System display, " << _number << ")." << endl;
  cout << "- A " << endl;
  if (_A) _A->display();
  else cout << "-> NULL" << endl;
  cout << "- b " << endl;
  if (_b) _b->display();
  else cout << "-> NULL" << endl;

  cout << "- M: " << endl;
  if (_M) _M->display();
  else cout << "-> NULL" << endl;

  cout << "============================================" << endl;
}

FirstOrderLinearTIDS* FirstOrderLinearTIDS::convert(DynamicalSystem* ds)
{
  FirstOrderLinearTIDS* lsds = dynamic_cast<FirstOrderLinearTIDS*>(ds);
  return lsds;
}
