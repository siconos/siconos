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
#ifndef FEMTOOLS
#define FEMTOOLS

#include <FENode.hpp>
#include <FiniteElementLinearTIDS.hpp>
#include <SiconosKernel.hpp>
#include <functional>
#include <optional>
#include <set>

/** @file FemTools.hpp toolbox of functions to handle finite element model with Siconos
 */

namespace siconos::mechanics::fem {

/** Names to identify properly gmsh tags*/
enum class MeshTags : int {
  bulk_material,
  boundary_conditions,
  applied_forces,
};

/** Map between explicit names (enum) and gmsh tag value (ing)*/
using Tags = std::map<MeshTags, int>;

/** Function type used to activate (or not) contact at a given node */
using ContactConditionFunction = std::function<bool(const siconos::mechanics::fem::FENode&)>;

/** @brief Build a dynamical system (Lagrangian linear) from a gmsh file
 *
 *  @param[in] gmsh_filename input file name (must be gmsh 2)
 *  @param[in] tags tags map (must fit with MeshTags allowed values).
 *  @param[in] material description of material properties of the elements
 *   (those tagged with tags[bulk_material])
 *  @param nodal_forces a vector of forces values to be applied on nodes listed with tag
 *   MeshTags::applied_forces. usefull only if applied_forces key is present in tags.
 *   Use {} if  not.
 *  @param bc_dof_index a vector of dof for which dirichlet BC must be applied (nodes with tag
 *   MeshTags::boundary_conditions) usefull only if boundary_conditions key is present in tags.
 *   Use {} if not.
 */
std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS>
build_dynamicalsystem_from_gmsh(std::string gmsh_filename, Tags& tags,
                                const Material& material,
                                const siconos::algebra::SiconosVector& nodal_forces,
                                const std::vector<int>& bc_dof_index);

/** @brief A class dedicated to automatic contact detection
 *   for a pre-defined FiniteElementLinearTIDS and its associated FiniteElementModel
 *
 */
class ContactDetection : public siconos::simulation::InteractionManager {
 protected:
  double initial_gap_{0.};

  /** Normal vector at contact */
  std::shared_ptr<siconos::algebra::SiconosVector> normal_{nullptr};

  /** Tangent vector at contact */
  std::shared_ptr<siconos::algebra::SiconosVector> tangent_{nullptr};

  /** Non smooth law  */
  std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw_{nullptr};

  /** Dimension */
  unsigned int contact_frame_dimension_ = 1;

  /** Associated finite element model */
  std::shared_ptr<siconos::mechanics::fem::FiniteElementModel> femodel_{nullptr};

  /** Associated dynamical system */
  std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> fesolid_{nullptr};

  /** Set of nodes (from finite element model) that may be in contact */
  std::vector<std::shared_ptr<siconos::mechanics::fem::FENode>> contacting_node_zone_ = {};

  /** Function used to activate or not a contact */
  ContactConditionFunction contact_condition_{nullptr};

  /** Internal structure used to define an active contact node */
  struct ContactingNode {
    std::shared_ptr<siconos::mechanics::fem::FENode> node_{nullptr};
    std::shared_ptr<siconos::modeling::Interaction> inter_{nullptr};
    std::size_t node_index{0};

    ContactingNode(std::shared_ptr<siconos::mechanics::fem::FENode> node, size_t index)
        : node_(node), node_index(index) {}
  };

  /** Set of nodes where contact is currently active */
  std::set<std::shared_ptr<ContactingNode>> contacting_nodes_ = {};

 public:
  /** @brief constructor
   *  @param[in] initial_gap
   *  @param[in] normal normal vector
   *  @param[in] nslaw nonsmooth law
   *  @param[in] fesolid dynamical system
   *  @param[in] func lambda function used to activate (or not) the contact for a node
   *  @param[in, optional] contact_condition_tag tag (in gmsh file) to identify nodes where
   *   contact could occur If not set, use all nodes of the FE model.
   */
  ContactDetection(double initial_gap, std::shared_ptr<siconos::algebra::SiconosVector> normal,
                   std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw,
                   std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> fesolid,
                   const ContactConditionFunction& func,
                   const std::optional<int>& contact_condition_tag = std::nullopt);

  virtual ~ContactDetection() noexcept = default;

  // shoud be done with find and a lambda function
  std::shared_ptr<ContactingNode> find_contacting_node(std::size_t node_index);

  /** Called by Simulation after updating positions prior to starting
   * the Newton loop. */
  void updateInteractions(std::shared_ptr<siconos::simulation::Simulation> simulation);
};

// Run simulation, work for all examples with T3 FEM
int run_T3_simulation(
    std::shared_ptr<siconos::simulation::TimeStepping> simulation,
    std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> FEsolid,
    std::string basename, std::string reference_file_name);

}  // namespace siconos::mechanics::fem
#endif