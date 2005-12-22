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
#include "LatinSolver.h"
using namespace std;

// default/from data constructor
LatinSolver::LatinSolver(const string& solvingForm, const unsigned int& iter, const double& tol, const double& sd):
  Solver("Latin", solvingForm), maxIter(iter), tolerance(tol), searchDirection(sd)
{
  // fill solvingMethod structure
  setSolvingMethod();
}

// copy constructor
LatinSolver::LatinSolver(const LatinSolver& solv):
  Solver("Latin", solv.getSolvingFormalisation()), maxIter(solv.getMaxIter()), tolerance(solv.getTolerance()), searchDirection(solv.getSearchDirection())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

LatinSolver::LatinSolver(SolverXML* latinXML):
  Solver(latinXML), maxIter((static_cast<LatinSolverXML*>(latinXML))->getMaxIter()),
  tolerance((static_cast<LatinSolverXML*>(latinXML))->getTolerance()),
  searchDirection((static_cast<LatinSolverXML*>(latinXML))->getSearchDirection())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

LatinSolver::~LatinSolver()
{}

void LatinSolver::display() const
{
  Solver::display();
  cout << " - MaxIter: " << maxIter << endl;
  cout << " - Tolerance: " << tolerance << endl;
  cout << " - Search direction: " << searchDirection << endl;
  cout << "===== End of Solver display =====" << endl;
}

void LatinSolver::setSolvingMethod()
{
  if (solvingFormalisation == "LcpSolving")
  {
    strcpy(solvingMethod->lcp.name, solverAlgorithmName.c_str());
    solvingMethod->lcp.itermax = maxIter;
    solvingMethod->lcp.tol = tolerance;
    solvingMethod->lcp.k_latin = searchDirection;
  }
  else if (solvingFormalisation == "PrimalRelaySolving")
  {
    strcpy(solvingMethod->pr.name, solverAlgorithmName.c_str());
    solvingMethod->pr.itermax = maxIter;
    solvingMethod->pr.tol = tolerance;
    solvingMethod->pr.k_latin = searchDirection;
  }
  else if (solvingFormalisation == "DualRelaySolving")
  {
    strcpy(solvingMethod->dr.name, solverAlgorithmName.c_str());
    solvingMethod->dr.itermax = maxIter;
    solvingMethod->dr.tol = tolerance;
    solvingMethod->dr.k_latin = searchDirection;
  }
  else if (solvingFormalisation == "FrictionContact2DSolving")
  {
    strcpy(solvingMethod->pfc_2D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_2D.itermax = maxIter;
    solvingMethod->pfc_2D.tol = tolerance;
    solvingMethod->pfc_3D.k_latin = searchDirection;
  }
  else if (solvingFormalisation == "FrictionContact3DSolving")
  {
    strcpy(solvingMethod->pfc_3D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_3D.itermax = maxIter;
    solvingMethod->pfc_3D.tol = tolerance;
    solvingMethod->pfc_3D.k_latin = searchDirection;
  }
  else
    RuntimeException::selfThrow("LatinSolver constructor - solving method " + solvingFormalisation + " does not exist.");
}

LatinSolver* LatinSolver::convert(Solver* solv)
{
  cout << "LatinSolver::convert (Solver* nsl)" << endl;
  LatinSolver* latin = dynamic_cast<LatinSolver*>(solv);
  return latin;
}


