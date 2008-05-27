/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
#include "Numerics_Options.h"

/* Default value for verbose mode: turned to off
Warning: global variable
*/
int verbose = 0;

void setNumericsOptions(Numerics_Options* opt)
{
  verbose = opt->verboseMode;
}

void numericsError(char * functionName, char* message)
{
  char output[200] = "Numerics error - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, output);
  exit(EXIT_FAILURE);
}

void numericsWarning(char * functionName, char* message)
{
  char output[200] = "Numerics warning - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, output);
  exit(EXIT_FAILURE);
}
