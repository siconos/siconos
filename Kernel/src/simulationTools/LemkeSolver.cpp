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
#include "LemkeSolver.h"
using namespace std;

// default/from data constructor
LemkeSolver::LemkeSolver(const string& solvingForm, const unsigned int& iter):
  Solver("Lemke", solvingForm), maxIter(iter)
{
  // fill solvingMethod structure
  setSolvingMethod();
}

// copy constructor
LemkeSolver::LemkeSolver(const LemkeSolver& solv):
  Solver("Lemke", solv.getSolvingFormalisation()), maxIter(solv.getMaxIter())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

LemkeSolver::LemkeSolver(SolverXML* lemkeXML):
  Solver(lemkeXML), maxIter((static_cast<LemkeSolverXML*>(lemkeXML))->getMaxIter())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

LemkeSolver::~LemkeSolver()
{}

void LemkeSolver::display() const
{
  Solver::display();
  cout << " - Maximum number of iterations: " << maxIter << endl;
  cout << "===== End of Solver display =====" << endl;
}

void LemkeSolver::setSolvingMethod()
{
  if (solvingFormalisation == "LcpSolving")
  {
    strcpy(solvingMethod->lcp.name, solverAlgorithmName.c_str());
    solvingMethod->lcp.itermax = maxIter;
  }
  else
    RuntimeException::selfThrow("LemkeSolver constructor - solving method " + solvingFormalisation + " does not exist.");
}

//void LemkeSolver::saveNonSmoothLawToXML()
// {
// }

LemkeSolver* LemkeSolver::convert(Solver* solv)
{
  cout << "LemkeSolver::convert (Solver* nsl)" << endl;
  LemkeSolver* lemke = dynamic_cast<LemkeSolver*>(solv);
  return lemke;
}


