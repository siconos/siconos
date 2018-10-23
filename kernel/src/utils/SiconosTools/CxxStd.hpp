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
/*! \file CxxStd.hpp
  \brief Management of different c++ standards and compiler
*/



// Proper definition of isnan
#ifndef SICONOS_ISNAN
#define SICONOS_ISNAN
#if __cplusplus >= 201103L
#include <cmath>
#else
#if ((!defined(_MSC_VER)) && (!defined( __SUNPRO_CC)))
#include <cmath>
using std::isnan;
using std::isinf;
#else
#include <math.h>
#endif
#endif
#endif
