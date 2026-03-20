/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

#include "FemTools.hpp"

#include <FiniteElementModel.hpp>
#include <Material.hpp>
#include <MeshUtils.hpp>
#include <NodeFem1d2DR.hpp>
#include <NodeFem2d2DR.hpp>
#include <SiconosKernel.hpp>
#include <SiconosVector.hpp>
#include <Simulation.hpp>
#include <TimeStepping.hpp>
#include <Tools.hpp>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <memory>
#include <optional>

std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS>
siconos::mechanics::fem::build_dynamicalsystem_from_gmsh(
    std::string gmsh_filename, Tags& tags, const Material& material,
    const siconos::algebra::SiconosVector& nodal_forces,
    const std::vector<int>& bc_dof_index) {
  // Read mesh from input file
  auto mesh = siconos::mechanics::fem::createMeshFromGMSH2(gmsh_filename);

  std::filesystem::path fname(gmsh_filename);
  auto basename = fname.stem();
  siconos::mechanics::fem::writeMeshforPython(*mesh, basename);

  std::map<int, const siconos::mechanics::fem::Material> materials = {
      {tags[MeshTags::bulk_material], material}};

  auto start = std::chrono::system_clock::now();
  auto FEsolid =
      std::make_shared<siconos::mechanics::fem::FiniteElementLinearTIDS>(mesh, materials);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "FiniteElementLinearTIDS - Assembly time: " << elapsed << " ms\n";

  //  Applied forces
  if (tags.find(MeshTags::applied_forces) != tags.end())
    FEsolid->applyNodalForces(tags[MeshTags::applied_forces], nodal_forces);

  // Boundary Conditions  */
  /* This part should be hidden in a new BC function for a node number
   * and a dof index. */
  if (tags.find(MeshTags::boundary_conditions) != tags.end()) {
    FEsolid->applyDirichletBoundaryConditions(tags[MeshTags::boundary_conditions],
                                              bc_dof_index);
    //    FEsolid->boundaryConditions()->display();
  }
  return FEsolid;  // RVO
}

siconos::mechanics::fem::ContactDetection::ContactDetection(
    double initial_gap, const siconos::algebra::SiconosVector& normal,
    std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw,
    std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> fesolid,
    const siconos::mechanics::fem::ContactConditionFunction& func,
    const std::optional<int>& contact_condition_tag)
    : InteractionManager(),
      initial_gap_{initial_gap},
      normal_{normal},
      nslaw_{nslaw},
      fesolid_(fesolid),
      contact_condition_{func} {
  femodel_ = fesolid_->FEModel();

  auto nslaw_with_friction =
      std::dynamic_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(nslaw);
  if (nslaw_with_friction) contact_frame_dimension_ = 2;

  // Set the list of nodes potentially concerned by a contact
  if (contact_condition_tag)
    contacting_node_zone_ = femodel_->contactingNodes(*contact_condition_tag);
  else  // all nodes
  {
    contacting_node_zone_.assign(femodel_->nodes().begin(), femodel_->nodes().end());
  }

  if (contact_frame_dimension_ == 2)  // create tangent vertor
  {
    tangent_.resize(2);
    tangent_(0) = -normal_(1);
    tangent_(1) = normal_(0);
  }
}

std::shared_ptr<siconos::mechanics::fem::ContactDetection::ContactingNode>
siconos::mechanics::fem::ContactDetection::find_contacting_node(size_t node_index) {
  auto it = std::find_if(contacting_nodes_.begin(), contacting_nodes_.end(),
                         [node_index](std::shared_ptr<ContactingNode> node) {
                           return node_index == node->node_index;
                         });
  if (it != contacting_nodes_.end())
    return *it;
  else
    return nullptr;
}

void siconos::mechanics::fem::ContactDetection::updateInteractions(
    std::shared_ptr<siconos::simulation::Simulation> simulation) {
  // std::cout<< "\nCall to updateInteractions in ContactDetection" <<
  // std::endl;

  auto displacement = fesolid_->q_read();
  auto solid = simulation->nonSmoothDynamicalSystem();

  // Update the list of contacting node by brute contact detection
  // contacting_nodes_.clear();
  for (auto node : contacting_node_zone_) {
    if (contact_condition_(*node)) {
      auto dof_index = node->global_dof_index()[0];
      if (!find_contacting_node(dof_index)) {
        auto cn = std::make_shared<ContactingNode>(node, node->global_dof_index()[0]);
        contacting_nodes_.insert(cn);
      }
    } else {  // No contact -> remove not if present in the vector of contacting nodes
      auto cn = find_contacting_node(node->global_dof_index()[0]);
      if (cn) contacting_nodes_.erase(cn);
    }
  }

  // update/create Interaction for contacting points
  for (auto cn : contacting_nodes_) {
    auto node = cn->node_;
    auto node_idx = node->global_dof_index()[0];
    siconos::algebra::SiconosVector2 pc2;
    pc2(0) = cn->node_->x() + displacement(node_idx);
    pc2(1) = -initial_gap_;
    if (cn->inter_)  // update interaction->relation
    {
      // std::cout << "update interaction" << std::endl;

      if (contact_frame_dimension_ == 2) {
        auto r = std::dynamic_pointer_cast<siconos::mechanics::fem::NodeFem2d2DR>(
            cn->inter_->relation());
        assert(r);
        r->updateContactPoint(pc2, normal_, tangent_);
      } else {
        auto r = std::dynamic_pointer_cast<siconos::mechanics::fem::NodeFem1d2DR>(
            cn->inter_->relation());
        assert(r);
        r->updateContactPoint(pc2, normal_);
      }

    } else  // create an interaction and link
    {
      // std::cout << "create interaction" << std::endl;
      std::shared_ptr<siconos::modeling::Relation> relation;
      if (contact_frame_dimension_ == 2) {
        relation = std::make_shared<siconos::mechanics::fem::NodeFem2d2DR>(cn->node_, pc2,
                                                                           normal_, tangent_);
      } else {
        relation =
            std::make_shared<siconos::mechanics::fem::NodeFem1d2DR>(cn->node_, pc2, normal_);
      }
      auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw_, relation);
      cn->inter_ = inter;
      // link the interaction and the dynamical system
      solid->link(inter, fesolid_);
    }
  }
}
