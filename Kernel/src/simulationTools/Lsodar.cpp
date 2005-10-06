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

#include "Lsodar.h"
using namespace std;

Lsodar::Lsodar(OneStepIntegratorXML* osiXML): OneStepIntegrator(osiXML)
{
  integratorType = LSODAR_INTEGRATOR;
  RuntimeException::selfThrow("Lsodar::Lsodar() - xml constructor - not yet implemented");
}

Lsodar::Lsodar(TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(td, ds)
{
  integratorType = LSODAR_INTEGRATOR;
}

Lsodar::~Lsodar()
{}

void Lsodar::computeFreeState()
{
  IN("Lsodar::computeFreeState\n");
  integrate();
  OUT("Lsodar::computeFreeState\n");
}

void Lsodar::integrate()
{
  IN("Lsodar::integrate\n");
  //tryfunction(ds->getVectorFieldPtr());
  OUT("Lsodar::integrate\n");
}

void tryfunction(fctPtr f)
{
  printf("in testFunctionInParameter \n") ;

  //typedef void (* fctPtr)(int *sizeOfX, double *time, double *x, double *xdot);

  //  double x = [2];
  //  double time = 0;
  //  double xDot[2] = 0;
  //  int sizeOfX = 2;

  double *x;
  double *time;
  double *xDot;
  int *sizeOfX;

  x = (double*)malloc(sizeof(double) * 2);
  xDot = (double*)malloc(sizeof(double) * 2);
  time = (double*)malloc(sizeof(double));
  sizeOfX = (int*)malloc(sizeof(int));

  *sizeOfX = 2;
  x[0] = 1;
  x[1] = -1;
  printf("sizeOfX = %i -- %d\n", *sizeOfX, *sizeOfX);

  f(sizeOfX, time, x, xDot);
  printf("out testFunctionInParameter \n")  ;

  free(x);
  free(xDot);
  free(time);
  free(sizeOfX);
}

Lsodar* Lsodar::convert(OneStepIntegrator* osi)
{
  cout << "Lsodar::convert (OneStepIntegrator* osi)" << endl;
  Lsodar* lsodar = dynamic_cast<Lsodar*>(osi);
  return lsodar;
}


Lsodar::Lsodar(): OneStepIntegrator()
{
  integratorType = LSODAR_INTEGRATOR;
}
