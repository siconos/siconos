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

Solver::Solver(): nonSmoothPbType("undefined"), solverAlgorithmName(DEFAULT_SOLVER), solvingMethod(NULL), maxIter(DEFAULT_ITER), tolerance(DEFAULT_TOL),
  normType(DEFAULT_NORMTYPE), searchDirection(DEFAULT_SEARCHDIR)
{}

Solver::Solver(const string& nspType, const string& algoName, const unsigned int & iter, const double & tol,
               const string & norm, const double & searchDir):
  nonSmoothPbType(nspType), solverAlgorithmName(algoName), solvingMethod(NULL), maxIter(iter), tolerance(tol),
  normType(norm), searchDirection(searchDir)
{
  solvingMethod = new method();
  setSolvingMethod();
}

Solver::Solver(const Solver& newS):
  nonSmoothPbType(newS.getNonSmoothPbType()), solverAlgorithmName(newS.getSolverAlgorithmName()), solvingMethod(NULL),
  maxIter(newS.getMaxIter()), tolerance(newS.getTolerance()), normType(newS.getNormType()), searchDirection(newS.getSearchDirection())
{
  solvingMethod = new method();
  setSolvingMethod();
}

Solver::Solver(SolverXML* solvXml, const string& nspbType):
  nonSmoothPbType(nspbType), solverAlgorithmName("undefined"), solvingMethod(NULL), maxIter(DEFAULT_ITER), tolerance(DEFAULT_TOL),
  normType(DEFAULT_NORMTYPE), searchDirection(DEFAULT_SEARCHDIR)
{
  if (solvXml != NULL)
  {
    solverAlgorithmName = solvXml->getType();
    maxIter = solvXml->getMaxIter();
    tolerance = solvXml->getTolerance();
    normType = solvXml->getNormType();
    searchDirection = solvXml->getSearchDirection();
  }
  else RuntimeException::selfThrow("Solver:: xml constructor, xml file=NULL");
  solvingMethod = new method();
  setSolvingMethod();
}

Solver::~Solver()
{
  if (solvingMethod != NULL) delete solvingMethod;
  solvingMethod = NULL;
}

void Solver::display() const
{
  cout << "=== Solver data display ===" << endl;
  cout << " - Solver algorithm is: " << solverAlgorithmName << endl;
  cout << " - MaxIter: " << maxIter << endl;
  cout << " - Tolerance: " << tolerance << endl;
  cout << " - Search direction: " << searchDirection << endl;
  // cout<<" - Norm Type: " << normType <<endl; never used at the time.
  cout << "Warning: depending on the solver type, some of the data above may be useless. See documentation for more details." << endl;
  cout << "===== End of Solver display =====" << endl;
}

void Solver::saveSolverToXML()
{
  RuntimeException::selfThrow("saveSolverToXML: not yet implemented");
}

void Solver::setSolvingMethod()
{
  if (nonSmoothPbType == "LCP")
  {
    strcpy(solvingMethod->lcp.name, solverAlgorithmName.c_str());
    solvingMethod->lcp.itermax = maxIter;
    solvingMethod->lcp.tol = tolerance;
    solvingMethod->lcp.k_latin = searchDirection;
  }
  else if (nonSmoothPbType == "Relay")
  {
    strcpy(solvingMethod->pr.name, solverAlgorithmName.c_str());
    solvingMethod->pr.itermax = maxIter;
    solvingMethod->pr.tol = tolerance;
    solvingMethod->pr.k_latin = searchDirection;
  }
  else if (nonSmoothPbType == "DualRelay")
  {
    strcpy(solvingMethod->dr.name, solverAlgorithmName.c_str());
    solvingMethod->dr.itermax = maxIter;
    solvingMethod->dr.tol = tolerance;
    solvingMethod->dr.k_latin = searchDirection;
  }
  else if (nonSmoothPbType == "FrictionContact2D")
  {
    strcpy(solvingMethod->pfc_2D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_2D.itermax = maxIter;
    solvingMethod->pfc_2D.tol = tolerance;
    solvingMethod->pfc_2D.k_latin = searchDirection;
  }
  else if (nonSmoothPbType == "FrictionContact3D")
  {
    strcpy(solvingMethod->pfc_3D.name, solverAlgorithmName.c_str());
    solvingMethod->pfc_3D.itermax = maxIter;
    solvingMethod->pfc_3D.tol = tolerance;
    solvingMethod->pfc_3D.k_latin = searchDirection;
  }
  else if (nonSmoothPbType == "dualFrictionContact2D")
  {
    strcpy(solvingMethod->dfc_2D.name, solverAlgorithmName.c_str());
    solvingMethod->dfc_2D.itermax = maxIter;
    solvingMethod->dfc_2D.tol = tolerance;
    solvingMethod->dfc_2D.k_latin = searchDirection;
  }
  else
    RuntimeException::selfThrow("Solver constructor - non smooth formalisation of type " + nonSmoothPbType + " not available.");
}
