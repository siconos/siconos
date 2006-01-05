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
#include "LexicoLemkeSolver.h"
using namespace std;

// default/from data constructor
LexicoLemkeSolver::LexicoLemkeSolver(const string& solvingForm, const unsigned int& iter):
  Solver("LexicoLemke", solvingForm), maxIter(iter)
{
  // fill solvingMethod structure
  setSolvingMethod();
}

// copy constructor
LexicoLemkeSolver::LexicoLemkeSolver(const LexicoLemkeSolver& solv):
  Solver("LexicoLemke", solv.getSolvingFormalisation()), maxIter(solv.getMaxIter())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

LexicoLemkeSolver::LexicoLemkeSolver(SolverXML* lemkeXML):
  Solver(lemkeXML), maxIter((static_cast<LexicoLemkeSolverXML*>(lemkeXML))->getMaxIter())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

LexicoLemkeSolver::~LexicoLemkeSolver()
{}

void LexicoLemkeSolver::display() const
{
  Solver::display();
  cout << " - Maximum number of iterations: " << maxIter << endl;
  cout << "===== End of Solver display =====" << endl;
}

void LexicoLemkeSolver::setSolvingMethod()
{
  if (solvingFormalisation == "LcpSolving")
  {
    strcpy(solvingMethod->lcp.name, solverAlgorithmName.c_str());
    solvingMethod->lcp.itermax = maxIter;
    solvingMethod->lcp.chat = 1;
  }
  else
    RuntimeException::selfThrow("LexicoLemkeSolver constructor - solving method " + solvingFormalisation + " does not exist.");
}

LexicoLemkeSolver* LexicoLemkeSolver::convert(Solver* solv)
{
  cout << "LexicoLemkeSolver::convert (Solver* nsl)" << endl;
  LexicoLemkeSolver* lemke = dynamic_cast<LexicoLemkeSolver*>(solv);
  return lemke;
}


