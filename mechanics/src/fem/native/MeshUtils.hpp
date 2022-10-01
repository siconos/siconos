/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include "Mesh.hpp"
#include "FiniteElementModel.hpp"
#include <SiconosKernel.hpp>
#include <string.h>

namespace siconos::mechanics::fem::native
{

Mesh* create2dMesh2x1();

Mesh* create2dMeshnxm(int n, int m, double Lx, double Ly);

Mesh* createMeshFromGMSH2(std::string gmsh_filename);

void  writeMeshforPython(std::shared_ptr<Mesh>  mesh);

std::string prepareWriteDisplacementforPython(std::string basename);

void  writeDisplacementforPython(std::shared_ptr<Mesh>  mesh,
                                 std::shared_ptr<FiniteElementModel> femodel,
                                 std::shared_ptr<SiconosVector> x, std::string filename);

}
#endif
