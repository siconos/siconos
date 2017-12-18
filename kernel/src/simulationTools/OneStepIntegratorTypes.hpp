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

/*! \file OneStepIntegratorTypes.hpp
  \brief enum of the available types for one-step time integrators.
*/

#ifndef OSITYPES_HPP
#define OSITYPES_HPP
#include "SiconosFwd.hpp"

/** Namespace for one-step integrators. */
namespace OSI
{

/** List of OneStepIntegrator types/ids */
enum TYPES
{
  /** Euler-Moreau scheme*/
  EULERMOREAUOSI,
  /** Moreau-Jean scheme */
  MOREAUJEANOSI,
  /** ?? */
  MOREAUJEANGOSI,
  /** LSodar (ode solver with rootfinding process) */
  LSODAROSI,
  /** odepack HEM5 (Hairer) solver */
  HEM5OSI,
  /** Moreau-Jean with direct projection */
  MOREAUDIRECTPROJECTIONOSI,
  /** Moreau-Jean with combined projection */
  MOREAUCOMBINEDPROJECTIONOSI,
  /** Intregrator based on discontinuous Galerkin methods*/
  D1MINUSLINEAROSI,
  /** Schatzman-Paoli scheme */
  SCHATZMANPAOLIOSI,
  /** Zero-order osi (?) */
  ZOHOSI,
  /** Newmark-like scheme*/
  NEWMARKALPHAOSI,
  /** Moreau-Jean-Bilbao */
  MOREAUJEANBILBAOOSI
};

}
#endif
