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
#include "FirstOrderR.hpp"

void FirstOrderR::initialize(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workVInter, VectorOfSMatrices& workMInter)
{
  workVInter.resize(FirstOrderR::workVecSize);
  workMInter.resize(FirstOrderR::mat_workMatSize);
  initComponents(inter, DSlink, workVInter, workMInter);
  unsigned int sizeY = inter.getSizeOfY();
  
  if (!workVInter[FirstOrderR::h_alpha])
        workVInter[FirstOrderR::h_alpha].reset(new SiconosVector(sizeY));
  if (!workVInter[FirstOrderR::vec_residuY])
        workVInter[FirstOrderR::vec_residuY].reset(new SiconosVector(sizeY));

  if (requireResidu())
    {
      unsigned int sizeOfDS = inter.getSizeOfDS();
      if (!workVInter[FirstOrderR::g_alpha])
        workVInter[FirstOrderR::g_alpha].reset(new SiconosVector(sizeOfDS));
      if (!workVInter[FirstOrderR::vec_residuR])
        workVInter[FirstOrderR::vec_residuR].reset(new SiconosVector(sizeOfDS));
    }
}


