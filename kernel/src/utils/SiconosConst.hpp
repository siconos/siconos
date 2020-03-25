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
/*! \file SiconosConst.hpp
\brief General constants for Siconos Kernel.
*/
#ifndef __SICONOSCONST__
#define __SICONOSCONST__

/**
   Internal bound max levels for time integrators.
   This value may be checked to see if initialization has occured.
 */
#define LEVELMAX 999

/** double precision machine */
#define MACHINE_PREC std::numeric_limits<double>::epsilon()

// #ifndef nullptr
// const int nullptr = 0;
// #endif

#endif

