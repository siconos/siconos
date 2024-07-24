/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef OCC_UTILS
#define OCC_UTILS

#include <Standard_TypeDef.hxx>

#include "MechanicsFwd.hpp"
#include "SiconosFwd.hpp"

struct OccContactFace;
struct OccContactEdge;

void occ_move(TopoDS_Shape& shape, const SiconosVector& pos);

void occ_distanceFaceFace(const OccContactFace& csh1, const OccContactFace& csh2,
                          Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                          Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                          Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                          Standard_Real& MinDist);

void occ_distanceFaceEdge(const OccContactFace& csh1, const OccContactEdge& csh2,
                          Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                          Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                          Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                          Standard_Real& MinDist);

#endif
