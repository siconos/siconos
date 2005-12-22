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
 *
 * This file contains global variables defined in Siconos/Kernel, that are
 * used to set default values.
 *
 */
#ifndef __KERNELDEFAULTCONFIG__
#define __KERNELDEFAULTCONFIG__

#include <string>

const unsigned int MATRIX_MAX_SIZE = 10;
const unsigned int VECTOR_MAX_SIZE = 10;
const std::string FILE_STORAGE = "ascii";

const std::string XML_SCHEMA = "/share/SICONOS/SiconosModelSchema-V1.2.xsd";

// Default values for non smooth problem solver
const std::string DefaultSolvingForm = "LcpSolving";
const std::string DefaultAlgoName = "NSQP";
const std::string DefaultAlgoNormType = "default";
const double  DefaultAlgoTolerance = 0.01;
const unsigned int DefaultAlgoMaxIter = 10;
const double  DefaultAlgoSearchDirection = 1.0;

// Defaults values for input/output plug-in in DSIO
const std::string  DefaultComputeInput = "DefaultPlugin:computeInput";
const std::string  DefaultComputeOutput = "DefaultPlugin:computeOutput";

#endif

