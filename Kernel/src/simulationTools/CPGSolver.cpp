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
#include "CPGSolver.h"
using namespace std;

// default/from data constructor
CPGSolver::CPGSolver(const string& solvingForm, const unsigned int& iter, const double& tol):
  Solver("CPG", solvingForm), maxIter(iter), tolerance(tol)
{
  // fill solvingMethod structure
  setSolvingMethod();
}

// copy constructor
CPGSolver::CPGSolver(const CPGSolver& solv):
  Solver("CPG", solv.getSolvingFormalisation()), maxIter(solv.getMaxIter()), tolerance(solv.getTolerance())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

CPGSolver::CPGSolver(SolverXML* cpgXML):
  Solver(cpgXML), maxIter((static_cast<CPGSolverXML*>(cpgXML))->getMaxIter()), tolerance((static_cast<CPGSolverXML*>(cpgXML))->getTolerance())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

CPGSolver::~CPGSolver()
{}

void CPGSolver::display() const
{
  Solver::display();
  cout << " - MaxIter: " << maxIter << endl;
  cout << " - Tolerance: " << tolerance << endl;
  cout << "===== End of Solver display =====" << endl;
}

void CPGSolver::setSolvingMethod()
{
  if (solvingFormalisation == "LcpSolving")
  {
    strcpy(solvingMethod->lcp.name, solverAlgorithmName.c_str());
    solvingMethod->lcp.itermax = maxIter;
    solvingMethod->lcp.tol = tolerance;
    solvingMethod->lcp.chat = 0;
  }
  else if (solvingFormalisation == "FrictionContact2DSolving")
  {
    strcpy(solvingMethod->pfc_2D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_2D.itermax = maxIter;
    solvingMethod->pfc_2D.tol = tolerance;
    solvingMethod->pfc_2D.chat = 0;
  }
  else if (solvingFormalisation == "FrictionContact3DSolving")
  {
    strcpy(solvingMethod->pfc_3D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_3D.itermax = maxIter;
    solvingMethod->pfc_3D.tol = tolerance;
    solvingMethod->pfc_3D.chat = 0;
  }
  else
    RuntimeException::selfThrow("CPGSolver constructor - solving method " + solvingFormalisation + " not available for this solver.");
}

CPGSolver* CPGSolver::convert(Solver* solv)
{
  cout << "CPGSolver::convert (Solver* nsl)" << endl;
  CPGSolver* nsqp = dynamic_cast<CPGSolver*>(solv);
  return nsqp;
}


