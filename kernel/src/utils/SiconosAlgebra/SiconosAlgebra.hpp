/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/


/** ! \file SiconosAlgebra.hpp
    \brief Header file for Siconos Algebra objects

    This file provides typedef for matrix and vector objects, const values and so on ...
*/

#ifndef SiconosAlgebra
#define SiconosAlgebra

//#include "SiconosConfig.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/banded.hpp>

/** Those cannot be forward declare  */

/** iterator (type1) through a sparse mat */
//typedef SparseMat::iterator1 SpMatIt1;

/** iterator (type2) through a sparse mat */
//typedef SparseMat::iterator2 SpMatIt2;


#endif
