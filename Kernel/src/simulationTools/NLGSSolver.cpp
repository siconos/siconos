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
#include "NLGSSolver.h"
using namespace std;

// default/from data constructor
NLGSSolver::NLGSSolver(const string& solvingForm, const unsigned int& iter, const double& tol):
  Solver("NLGS", solvingForm), maxIter(iter), tolerance(tol)
{
  // fill solvingMethod structure
  setSolvingMethod();
}

// copy constructor
NLGSSolver::NLGSSolver(const NLGSSolver& solv):
  Solver("NLGS", solv.getSolvingFormalisation()), maxIter(solv.getMaxIter()), tolerance(solv.getTolerance())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

NLGSSolver::NLGSSolver(SolverXML* nlgsXML):
  Solver(nlgsXML), maxIter((static_cast<NLGSSolverXML*>(nlgsXML))->getMaxIter()), tolerance((static_cast<NLGSSolverXML*>(nlgsXML))->getTolerance())
{
  // fill solvingMethod structure
  setSolvingMethod();
}

NLGSSolver::~NLGSSolver()
{}

void NLGSSolver::display() const
{
  Solver::display();
  cout << " - MaxIter: " << maxIter << endl;
  cout << " - Tolerance: " << tolerance << endl;
  cout << "===== End of Solver display =====" << endl;
}

void NLGSSolver::setSolvingMethod()
{
  if (solvingFormalisation == "LcpSolving")
  {
    strcpy(solvingMethod->lcp.name, solverAlgorithmName.c_str());
    solvingMethod->lcp.itermax = maxIter;
    solvingMethod->lcp.tol = tolerance;
  }
  else if (solvingFormalisation == "PrimalRelaySolving")
  {
    strcpy(solvingMethod->pr.name, solverAlgorithmName.c_str());
    solvingMethod->pr.itermax = maxIter;
    solvingMethod->pr.tol = tolerance;
  }
  else if (solvingFormalisation == "DualRelaySolving")
  {
    strcpy(solvingMethod->dr.name, solverAlgorithmName.c_str());
    solvingMethod->dr.itermax = maxIter;
    solvingMethod->dr.tol = tolerance;
  }
  else if (solvingFormalisation == "FrictionContact2DSolving")
  {
    strcpy(solvingMethod->pfc_2D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_2D.itermax = maxIter;
    solvingMethod->pfc_2D.tol = tolerance;
  }
  else if (solvingFormalisation == "FrictionContact3DSolving")
  {
    strcpy(solvingMethod->pfc_3D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_3D.itermax = maxIter;
    solvingMethod->pfc_3D.tol = tolerance;
  }
  else
    RuntimeException::selfThrow("NLGSSolver constructor - solving method " + solvingFormalisation + " does not exist.");
}

NLGSSolver* NLGSSolver::convert(Solver* solv)
{
  cout << "NLGSSolver::convert (Solver* nsl)" << endl;
  NLGSSolver* nsqp = dynamic_cast<NLGSSolver*>(solv);
  return nsqp;
}


