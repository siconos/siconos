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

/*! \file SolidLinearTIDS.hpp

 */
#ifndef SOLIDLAGRANGIANTIDS_H
#define SOLIDLAGRANGIANTIDS_H

#include "FiniteElementLinearTIDS.hpp"

/** Finite Element discretization of an elastic solids that inherits from FEM with time invariant coefficients and exhibits sigma
 *  to enable easy plasticty implementation using a stress relation.
 * - \f$M\dot v + K u +  B^T \sigma = F_{ext}(t,z) + p
 * v = \dot u
 * \dot \sigma = B v + \dot \epsilon_p \f$
 */
namespace siconos::mechanics::fem
{
class SolidLinearTIDS : public FiniteElementLinearTIDS
{

protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);
//  using FiniteElementLinearTIDS::_mesh;
//  using FiniteElementLinearTIDS::_materials;
//  using FiniteElementLinearTIDS::_storageType;
//  /** M Matrix */
//  std::shared_ptr<siconos::algebra::SimpleMatrix> _M;
  /** Elasticity Matrix */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _S;
  /** B matrix from FEM */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _B;

  /** stress of the system */
  std::shared_ptr<siconos::algebra::SiconosVector> _sigma{nullptr};
  /** Initial stress of the system */
  std::shared_ptr<siconos::algebra::SiconosVector> _sigma0{nullptr};

  std::shared_ptr<siconos::algebra::SiconosVector> _plasticRate = {nullptr};

  unsigned int _dimStress{0};
  /** memory of previous coordinates of the system */
  siconos::algebra::SiconosMemory _stressMemory;
  siconos::algebra::SiconosMemory _plasticRateMemory;

public:

  /** constructor from initial state and all matrix operators.
   * \param mesh the mesh that defined the spatial discretization
   * \param material
   */
  SolidLinearTIDS(std::shared_ptr<Mesh> mesh,
                          std::map<unsigned int, std::shared_ptr<Material> > materials,
                          siconos::algebra::UblasType storageType=siconos::algebra::UblasType::DENSE);


  /** destructor */
  ~SolidLinearTIDS() noexcept = default;


  /** Compute kinetic energy
   */
  double kineticEnergy() const;
  double elasticPotentialEnergy() const;

  void display(bool brief) const override;

  inline std::shared_ptr<siconos::algebra::SiconosVector> stress() const
  {
    return _sigma;
  }


  /** get \dot \epsilon_p
   *
   *  \param level unsigned int, required level for p, default = 2
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> plasticRate() const
  {
    return _plasticRate;
  }


  /** get S matrix (pointer link)
   *
   *  \return pointer on a SiconosMatrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> S() const { return _S; }


  /** get B matrix (pointer link)
   *
   *  \return pointer on a SiconosMatrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> B() const { return _B; }

  inline unsigned int stressDimension()  { return _dimStress; }
  inline unsigned int dimension() const override { return _n; }
  inline unsigned int velocityDimension() const { return _ndof; }


  /** initialize the siconos::algebra::SiconosMemory objects with a positive size.
   *
   *  \param size the size of the siconos::algebra::SiconosMemory. must be >= 0
   */
  void initMemory(unsigned int size) override;

  inline const siconos::algebra::SiconosMemory &stressMemory()
  {
    return _stressMemory;
  }

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   *  \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory() override;

  modeling::Type acceptType(types::FindType &ft) const override { return ft.visit(*this); }


};

} // namespace siconos::mechanics::fem
#endif // SOLIDLAGRANGIANTIDS_H
