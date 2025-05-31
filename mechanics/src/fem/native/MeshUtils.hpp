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

/*! \file MeshUtils.hpp

 */
#ifndef MESHUTILS_H
#define MESHUTILS_H

#include "SiconosVector.hpp"
#include <memory>
#include <string>


namespace siconos::mechanics::fem {

class Mesh;
class FiniteElementModel;

std::shared_ptr<Mesh> create2dMesh2x1();

std::shared_ptr<Mesh> create2dMeshnxm(int n, int m, double Lx, double Ly);

std::shared_ptr<Mesh> createMeshFromGMSH2(std::string gmsh_filename);

void writeMeshforPython(std::shared_ptr<Mesh> mesh);

std::string prepareWriteDisplacementforPython(std::string basename);

void writeDisplacementforPython(std::shared_ptr<Mesh> mesh,
                                std::shared_ptr<FiniteElementModel> femodel,
                                std::shared_ptr<siconos::algebra::SiconosVector> x,
                                std::string filename);

}  // namespace siconos::mechanics::fem
#endif
