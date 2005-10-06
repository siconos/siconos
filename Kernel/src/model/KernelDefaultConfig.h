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
#ifndef __KERNELDEFAULTCONFIG__
#define __KERNELDEFAULTCONFIG__

#include <string>

/* A matrix is saved in the XML output file if his size is not higher then MatrixMaxSize */
/*const*/ unsigned int    MATRIX_MAX_SIZE = 10;
/*const*/
unsigned int    VECTOR_MAX_SIZE = 10;
/*const*/
std::string   FILE_STORAGE = "ascii"; // N_ASCII or N_BINARY

/*const*/
std::string   XML_SCHEMA = "/share/SICONOS/SiconosModelSchema-V1.2.xsd";

std::string   DefaultSolver = "default";
std::string   DefaultAlgoName = "default";
std::string   DefaultAlgoNormType = "default";
double  DefaultAlgoTolerance = -1.0;
int   DefaultAlgoMaxIter = -1;
double  DefaultAlgoSearchDirection = -1.0;

std::string  DefaultComputeInput = "BasicPlugin:computeInput";
std::string  DefaultComputeOutput = "BasicPlugin:computeOutput";

#endif

