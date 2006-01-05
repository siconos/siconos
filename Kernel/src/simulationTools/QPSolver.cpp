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
#include "QPSolver.h"
using namespace std;

// default/from data constructor
QPSolver::QPSolver(const string& solvingForm, const double& tol):
  Solver("QP", solvingForm), tolerance(tol)
{
  // fill solvingMethod structure
  setSolvingMethod();
}

// copy constructor
QPSolver::QPSolver(const QPSolver& solv):
  Solver("QP", solv.getSolvingFormalisation()), tolerance(solv.getTolerance())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

QPSolver::QPSolver(SolverXML* qpXML):
  Solver(qpXML), tolerance((static_cast<QPSolverXML*>(qpXML))->getTolerance())
{
  // fill solvingMethod structure
  setSolvingMethod();
}


QPSolver::~QPSolver()
{}

void QPSolver::display() const
{
  Solver::display();
  cout << " - Tolerance: " << tolerance << endl;
  cout << "===== End of Solver display =====" << endl;
}

void QPSolver::setSolvingMethod()
{
  if (solvingFormalisation == "LcpSolving")
  {
    strcpy(solvingMethod->lcp.name, solverAlgorithmName.c_str());
    solvingMethod->lcp.tol = tolerance;
  }
  else
    RuntimeException::selfThrow("QPSolver constructor - solving method " + solvingFormalisation + " not available for this solver.");
}

QPSolver* QPSolver::convert(Solver* solv)
{
  cout << "QPSolver::convert (Solver* nsl)" << endl;
  QPSolver* qp = dynamic_cast<QPSolver*>(solv);
  return qp;
}


