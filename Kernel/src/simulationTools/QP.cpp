/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

#include "QP.h"
#include "QPXML.h"

using namespace std;

QP::QP(OneStepNSProblemXML* osnspbxml, Simulation* newStrat):
  OneStepNSProblem("QP", osnspbxml, newStrat), Q(NULL), p(NULL),
  isQAllocatedIn(true), isPAllocatedIn(true)
{}

QP::~QP()
{}

void QP::compute(const double& time)
{
  RuntimeException::selfThrow("QP::compute not yet implemented");

}

void QP::display() const
{
  RuntimeException::selfThrow("QP::compute not yet implemented");
}

void QP::saveNSProblemToXML()
{
  RuntimeException::selfThrow("QP::compute not yet implemented");
}

void QP::savePToXML()
{
  RuntimeException::selfThrow("QP::compute not yet implemented");
}

void QP::saveQToXML()
{
  RuntimeException::selfThrow("QP::compute not yet implemented");
}

QP* QP::convert(OneStepNSProblem* osnsp)
{
  QP* qp = dynamic_cast<QP*>(osnsp);
  return qp;
}

// Default private constructor
QP::QP(): OneStepNSProblem(), Q(NULL), p(NULL),
  isQAllocatedIn(false), isPAllocatedIn(false)
{
  nspbType = "QP";
}
