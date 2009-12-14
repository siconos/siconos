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
#include "NonSmoothSolver.hpp"
#include "RuntimeException.hpp"
#include <iterator>
#include <iostream>

using namespace std;

// Private function used to fill in Solver_Options structure numerics_solver_options with current object data,
// only if isSet = true. Else set numerics_solver_options->isSet to false, which means that parameters
// will be read in the default_parameters file during Numerics driver call.
void NonSmoothSolver::fillSolverOptions()
{
  numerics_solver_options.reset(new Solver_Options);
  numerics_solver_options->filterOn = 0;
  if (!isSet)
    numerics_solver_options->isSet = 0;

  else
  {
    numerics_solver_options->isSet = 1;
    strcpy(numerics_solver_options->solverName, _name.c_str());
    // No memory allocation for iparam and dparam of numerics_solver_options, just pointer links
  }
  // Link is required for iparam and dparam even when isSet == 0, since we need to recover output parameters (error ...)
  numerics_solver_options->iSize = int_parameters->size();
  numerics_solver_options->dSize = double_parameters->size();
  numerics_solver_options->iparam = &((*int_parameters)[0]);
  numerics_solver_options->dparam = &((*double_parameters)[0]);
}


// Default constructor: build an empty Solver_Options structure
// Parameters will be read in the default input parameters file during Numerics driver call.
NonSmoothSolver::NonSmoothSolver(): _name("undefined"), isSet(false)
{
  int_parameters.reset(new IntParameters(NB_PARAM));
  double_parameters.reset(new DoubleParameters(NB_PARAM));
  fillSolverOptions();
}

// Copy constructor
NonSmoothSolver::NonSmoothSolver(const NonSmoothSolver& newS):
  _name(newS.name()), isSet(newS.isSolverSet())
{
  int_parameters.reset(new IntParameters(*newS.intParameters()));
  double_parameters.reset(new DoubleParameters(*newS.doubleParameters()));
  fillSolverOptions();
}

// Constructor from a set of data
NonSmoothSolver::NonSmoothSolver(const std::string& newName, IntParameters& iparam, DoubleParameters& dparam):
  _name(newName), isSet(true)
{
  int_parameters.reset(new IntParameters(iparam));
  double_parameters.reset(new DoubleParameters(dparam));
  fillSolverOptions();
}
NonSmoothSolver::NonSmoothSolver(const std::string& newName, IntParameters& iparam, DoubleParameters& dparam, double * dWork, int * iWork):
  _name(newName), isSet(true)
{
  int_parameters.reset(new IntParameters(iparam));
  double_parameters.reset(new DoubleParameters(dparam));
  fillSolverOptions();
  numerics_solver_options->dWork = dWork;
  numerics_solver_options->iWork = iWork;
}

// Construction using XML object
NonSmoothSolver::NonSmoothSolver(SP::NonSmoothSolverXML solvXML):
  _name("undefined"), isSet(true)
{
  if (! solvXML)
    RuntimeException::selfThrow("NonSmoothSolver, XML constructor, NULL input");

  // Read name
  _name = solvXML->name();
  // Built and read int_parameters and double_parameters
  int_parameters.reset(new IntParameters);
  double_parameters.reset(new DoubleParameters);
  if (!solvXML->hasIparam() || !solvXML->hasDparam())
    RuntimeException::selfThrow("NonSmoothSolver, XML constructor, missing input in XML file (int or double vector)");
  solvXML->getIParam(*int_parameters);
  solvXML->getDParam(*double_parameters);
  //   if(int_parameters->size()>NB_PARAM || double_parameters->size()>NB_PARAM)
  //     {
  //       std::cout << "NonSmoothSolver xml constructor warning: too large number of provided int and/or double parameters. "<< std::endl;
  //       std::cout << "Some of them might be ignored. Check in the solver documentation to know what are the required parameters." << std::endl;
  //     }
  fillSolverOptions();

}

// Construction by reading in a file (XML)
NonSmoothSolver::NonSmoothSolver(const string& inputFile):
  _name("undefined"), isSet(false)
{}

NonSmoothSolver::~NonSmoothSolver()
{}

void NonSmoothSolver::display() const
{
  cout << "=== Non Smooth Solver based on algorithm of type " << _name << " with:" << endl;
  if (isSet)
  {
    cout << "int parameters:  " ;
    copy(int_parameters->begin(), int_parameters->end(), ostream_iterator<int>(cout, " "));
    cout << endl;
    cout << "double parameters:  " ;
    copy(double_parameters->begin(), double_parameters->end(), ostream_iterator<double>(cout, " "));
    cout << endl;
  }
  else
    cout << "no user input parameters" << endl;
  cout << "===== End of NonSmoothSolver display =====" << endl;
}

void NonSmoothSolver::saveNonSmoothSolverToXML()
{
  RuntimeException::selfThrow("NonSmoothSolver, saveNonSmoothSolverToXML: not yet implemented");
}

