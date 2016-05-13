/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*! \file PluginTypes.hpp
  \brief list of typedef for pointers to functions used in plugin mechanism.
*/

#ifndef PLUGINTYPES_HPP
#define PLUGINTYPES_HPP

/** Pointer to function used for plug-in for matrix-type operators that depends only on time */
typedef void (*MatrixFunctionOfTime)(double, unsigned int, unsigned int, double*, unsigned int, double*);

/** Pointer to function used for plug-in for vector-type operators that depends only on time */
typedef void (*VectorFunctionOfTime)(double, unsigned int, double*, unsigned int, double*);

/** */
typedef void (*FPtr1)(double, unsigned int, double*, double*, unsigned int, double*);

/** */
typedef void (*FPtr2)(unsigned int, double*, unsigned int, double*, double*, unsigned int, double*);

/** */
typedef void (*FPtr3)(unsigned int, double*, unsigned int, double*, unsigned int, double*);

typedef void (*FPtr4bis)(unsigned int, double*, unsigned int, double*, unsigned int, double*, unsigned int, double*);

/** */
typedef void (*FPtr4)(unsigned int, double*, double, unsigned int, double*, unsigned int, double*);

/** */
typedef void (*FPtr5)(unsigned int, double*, double*, double*, unsigned int, double*);

typedef void (*FPtr5bis)(unsigned int, double*, unsigned int, double*, unsigned int, double*, unsigned int, double*);

/** */
typedef void (*FPtr6)(double, unsigned int, double*, double*, double*, unsigned int, double*);

/** */
typedef void (*FPtr7)(unsigned int, double*, double*, unsigned int, double*);

typedef void (*OutPtr)(unsigned int, double*, double, unsigned int, double*, double*, unsigned int, double*);

typedef void (*InPtr)(unsigned int, double*, double, unsigned int, double*, unsigned int, double*);

#endif
