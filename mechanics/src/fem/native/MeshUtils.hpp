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

#include <string>

#include "SiconosVector.hpp"

namespace siconos::mechanics::fem {

class Mesh;
class FiniteElementModel;

std::shared_ptr<Mesh> create2dMesh2x1();

std::shared_ptr<Mesh> create2dMeshnxm(int n, int m, double Lx, double Ly);

/** @brief read a mesh from a gmsh file and convert it into a siconos mesh
 *  @param fname gmsh file name
 *  @return a Siconos Mesh
 *
 */
std::shared_ptr<Mesh> createMeshFromGMSH2(std::string fname);

/** @brief write a python file describing the mesh
 *  @param[in] mesh input mesh
 *  @param[in] fname outputfile name (.py)
 */
void writeMeshforPython(const Mesh& mesh, std::string);

std::string prepareWriteDisplacementforPython(std::string basename);

void writeDisplacementforPython(const Mesh& mesh, const FiniteElementModel& femodel,
                                const siconos::algebra::SiconosVector& x,
                                std::string filename);

}  // namespace siconos::mechanics::fem
#endif
