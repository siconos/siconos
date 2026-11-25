// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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

std::shared_ptr<Mesh> create2dMeshnxm(size_t n, size_t m, double Lx, double Ly);

/** @brief read a mesh from a gmsh file and convert it into a siconos mesh
 *  @param fname gmsh file name
 *  @return a Siconos Mesh
 *
 */
std::shared_ptr<Mesh> createMeshFromGMSH2(const std::string& input, bool is_filename = true);

std::shared_ptr<Mesh> createBeamMesh(const siconos::algebra::SiconosVector3& coords_start,
                                     const siconos::algebra::SiconosVector3& coords_end,
                                     size_t nb_elements, int dim);

/** @brief write a python file describing the mesh
 *  @param[in] mesh input mesh
 *  @param[in] fname outputfile name (.py)
 */
void writeMeshforPython(const Mesh& mesh, std::string);

std::string prepareWriteForPython(std::string basename, std::string fieldname,
                                  const std::vector<std::string>& fields);

std::string prepareWriteDisplacementforPython(std::string basename);

void writeDisplacementforPython(const Mesh& mesh, const FiniteElementModel& femodel,
                                const siconos::algebra::SiconosVector& x,
                                std::string filename);

std::string prepareWriteTensorforPython(std::string basename, std::string tensorName);
std::string prepareWriteBeamTensorforPython(std::string basename, std::string tensorName);

void writeTensorforPython(const FiniteElementModel& femodel,
                          const siconos::algebra::SiconosVector& x, std::string filename,
                          std::string tensorName);

void prepareWriteBeamPositionforSOFA(std::string filename);
void prepareWriteBlockPositionforSOFA(std::string filename);

void writeBeamPositionforSOFA(const Mesh& mesh, const FiniteElementModel& femodel,
                              const siconos::algebra::SiconosVector& positions,
                              std::string filename, double t);
void writeBlockPositionforSOFA(const siconos::algebra::SiconosVector& positions,
                               std::string filename, double t);

}  // namespace siconos::mechanics::fem
#endif
