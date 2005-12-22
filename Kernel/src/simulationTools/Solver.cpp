/* Siconos version 1.0, Copyright INRIA 2005.
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
#include "Solver.h"
using namespace std;

Solver::Solver(const string& algoName, const string& solvingForm):
  solvingFormalisation(solvingForm), solverAlgorithmName(algoName), solvingMethod(NULL)
{
  solvingMethod = new method;
}

Solver::Solver(const Solver& newS):
  solvingFormalisation(newS.getSolvingFormalisation()), solverAlgorithmName(newS.getSolverAlgorithmName()), solvingMethod(NULL)
{
  solvingMethod = new method;
}

Solver::Solver(SolverXML* solvXml):
  solvingFormalisation(""), solverAlgorithmName(""), solvingMethod(NULL)
{
  if (solvXml != NULL)
  {
    solvingFormalisation = solvXml->getSolvingFormalisation();
    solverAlgorithmName  = solvXml->getSolverAlgorithmName();
  }
  else RuntimeException::selfThrow("Solver:: xml constructor, xml file=NULL");
  solvingMethod = new method;
}


Solver::~Solver()
{
  if (solvingMethod != NULL) delete solvingMethod;
  solvingMethod = NULL;
}

void Solver::display() const
{
  cout << "=== Solver data display ===" << endl;
  cout << " - Solving formalisation is: " << solvingFormalisation << endl;
  cout << " - Solver algorithm is: " << solverAlgorithmName << endl;
}

void Solver::saveSolverToXML()
{
  RuntimeException::selfThrow("saveSolverToXML: not yet implemented");
}
