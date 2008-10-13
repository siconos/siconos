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
#include "FirstOrderLinearTIR.h"
#include "Interaction.h"

using namespace std;
using namespace RELATION;

// Default (private) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(): FirstOrderLinearR()
{
  subType = LinearTIR;
}

// xml constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::RelationXML relxml):
  FirstOrderLinearR(relxml)
{
  subType = LinearTIR;
}

// Minimum data (C, B) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  FirstOrderLinearR()
{
  subType = LinearTIR;
  C.reset(new SimpleMatrix(newC));
  B.reset(new SimpleMatrix(newB));
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD,
    const SiconosMatrix& newF, const SimpleVector& newE, const SiconosMatrix& newB): FirstOrderLinearR()
{
  subType = LinearTIR;

  C.reset(new SimpleMatrix(newC));
  D.reset(new SimpleMatrix(newD));
  F.reset(new SimpleMatrix(newF));
  e.reset(new SimpleVector(newE));
  B.reset(new SimpleMatrix(newB));
  initPluginFlags(false);
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newB):
  FirstOrderLinearR(newC, newB)
{
  subType = LinearTIR;
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SimpleVector newE, SP::SiconosMatrix newB):
  FirstOrderLinearR(newC, newD, newF, newE, newB)
{
  subType = LinearTIR;
}

FirstOrderLinearTIR::~FirstOrderLinearTIR()
{}

void FirstOrderLinearTIR::computeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SP::SiconosVector y = interaction->getYPtr(0);
  SP::SiconosVector lambda = interaction->getLambdaPtr(0);

  // compute y
  if (C)
    prod(*C, *data["x"], *y);
  else
    y->zero();

  if (D)
    prod(*D, *lambda, *y, false);

  if (e)
    *y += *e;

  if (F)
    prod(*F, *data["z"], *y, false);
}

void FirstOrderLinearTIR::computeInput(double time, unsigned int level)
{
  if (B)
  {
    // We get lambda of the interaction (pointers)
    SP::SiconosVector lambda = interaction->getLambdaPtr(level);
    prod(*B, *lambda, *data["r"], false);
  }
}

FirstOrderLinearTIR* FirstOrderLinearTIR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearTIR*>(r);
}

