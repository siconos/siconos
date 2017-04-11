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
/*! \file ModelingTools.hpp
Include files related to modelingTools
Note that not all files from the current location are listed below, since some of them are already included inside the ones below.
*/

#include "EqualityConditionNSL.hpp"
#include "ComplementarityConditionNSL.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "RelayNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NormalConeNSL.hpp"

#include "Interaction.hpp"

#include "RelationTypes.hpp"
#include "FirstOrderType1R.hpp"
#include "FirstOrderType2R.hpp"

#include "LagrangianLinearTIDS.hpp"
#include "LagrangianLinearDiagonalDS.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonEulerFrom1DLocalFrameR.hpp"
#include "NewtonEulerFrom3DLocalFrameR.hpp"

#include "BoundaryCondition.hpp"
#include "HarmonicBC.hpp"
#include "FixedBC.hpp"



#include "NonSmoothDynamicalSystem.hpp"

