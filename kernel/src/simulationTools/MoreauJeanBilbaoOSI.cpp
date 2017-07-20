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
#include "CxxStd.hpp"

#include "MoreauJeanBilbaoOSI.hpp"
#include "LagrangianLinearDiagonalDS.hpp"
#include "BlockVector.hpp"
#include "OneStepNSProblem.hpp"
#include "TypeName.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include <debug.h>
#include "LagrangianR.hpp"
#include "NewtonImpactNSL.hpp"

using namespace RELATION;

// --- constructor from a set of data ---
MoreauJeanBilbaoOSI::MoreauJeanBilbaoOSI():
  OneStepIntegrator(OSI::MOREAUJEANBILBAOOSI)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =1;
  _levelMinForInput =0;
  _levelMaxForInput =1;
  _steps=1;
}

void MoreauJeanBilbaoOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  // Get work buffers from the graph
  // const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
  //if(!(checkOSI(dsv)))
  // RuntimeException::selfThrow("MoreauJeanBilbaoOSI::initializeDynamicalSystem(m,t,ds) - ds does not belong to the OSI.");
  VectorOfVectors& work_ds = *_initializeDSWorkVectors(ds);
  //VectorOfVectors& work_ds = *_dynamicalSystemsGraph->properties(dsv).workVectors;
  Type::Siconos dsType = Type::value(*ds);
  if(dsType != Type::LagrangianLinearDiagonalDS)
    RuntimeException::selfThrow("MoreauJeanBilbaoOSI::initializeDynamicalSystem - not yet implemented for Dynamical system of type : " + Type::name(*ds));

  work_ds.resize(MoreauJeanBilbaoOSI::WORK_LENGTH);
  // - work buffer, used to save vfree.
  work_ds[MoreauJeanBilbaoOSI::VFREE].reset(new SiconosVector(ds->dimension()));
  // Initialize the iteration matrix and other parameters of the osi (all ds-dependent)
  _initialize_iteration_matrix(ds);
  LagrangianLinearDiagonalDS& lldds = static_cast<LagrangianLinearDiagonalDS&> (*ds);
  // Update dynamical system components (for memory swap).
  lldds.computeForces(t, lldds.q(), lldds.velocity());
  lldds.swapInMemory();
}

void MoreauJeanBilbaoOSI::fillDSLinks(Interaction &inter, InteractionProperties& interaction_properties,
				  DynamicalSystemsGraph & DSG)
{
  // Get the dynamical systems linked by inter
  SP::DynamicalSystem ds1= interaction_properties.source;
  SP::DynamicalSystem ds2= interaction_properties.target;
  assert(ds1);
  assert(ds2);

  VectorOfVectors& workV = *interaction_properties.workVectors;
  workV.resize(MoreauJeanBilbaoOSI::WORK_INTERACTION_LENGTH);
  workV[MoreauJeanBilbaoOSI::OSNSP_RHS].reset(new SiconosVector(inter.getSizeOfY()));

  // work vector of the interaction (from interaction_properties)
  VectorOfBlockVectors& DSlink = *interaction_properties.DSlink;

  // // Relation type.
  Relation &relation =  *inter.relation();
  RELATION::TYPES relationType = relation.getType();
  RELATION::SUBTYPES relationSubType = inter.relation()->getSubType();

  if(relationType != Lagrangian || relationSubType != LinearTIR)
    RuntimeException::selfThrow("MoreauJeanBilbaoOSI::computeFreeOutput only Lagrangian Linear Relations are allowed.");

  // // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

  // allocate and set work vectors for the osi

  if(checkOSI(DSG.descriptor(ds1)))
  {
    VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
    if (!DSlink[LagrangianR::xfree])
    {
      DSlink[LagrangianR::xfree].reset(new BlockVector());
      DSlink[LagrangianR::xfree]->insertPtr(workVds1[MoreauJeanBilbaoOSI::VFREE]);
    }
    else
    {
      DSlink[LagrangianR::xfree]->setVectorPtr(0,workVds1[MoreauJeanBilbaoOSI::VFREE]);
    }
  }
  if (ds1 != ds2)
  {
    if(checkOSI(DSG.descriptor(ds2)))
    {
      VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      if (!DSlink[LagrangianR::xfree])
      {
        DSlink[LagrangianR::xfree].reset(new BlockVector());
        //dummy insertion to reserve first vector for ds1
        DSlink[LagrangianR::xfree]->insertPtr(workVds2[MoreauJeanBilbaoOSI::VFREE]);
        DSlink[LagrangianR::xfree]->insertPtr(workVds2[MoreauJeanBilbaoOSI::VFREE]);
      }
      else
      {
        DSlink[LagrangianR::xfree]->insertPtr(workVds2[MoreauJeanBilbaoOSI::VFREE]);
      }
    }
  }
}

void MoreauJeanBilbaoOSI::initialize_nonsmooth_problems()
{
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->setIndexSetLevel(1);
  ((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY])->setInputOutputLevel(1);
}

void MoreauJeanBilbaoOSI::_initialize_iteration_matrix(SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical system, depending on its type.
  //
  // W = mass + time_step**2/2(I - Theta)Stiffness + time_step * Sigma*

  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
  VectorOfVectors& work = *_dynamicalSystemsGraph->properties(dsv).workVectors;

  if(work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]
     || work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR]
     || _dynamicalSystemsGraph->properties(dsv).W)
    RuntimeException::selfThrow("MoreauJeanBilbaoOSI::_initialize_iteration_matrix(t,ds) - W has already been initialized by another osi");

  Type::Siconos dsType = Type::value(*ds);
  if(dsType != Type::LagrangianLinearDiagonalDS)
    RuntimeException::selfThrow("MoreauJeanBilbaoOSI::initialize_iteration_matrix - not yet implemented for Dynamical system of type : " + Type::name(*ds));
  LagrangianLinearDiagonalDS& lldds = static_cast<LagrangianLinearDiagonalDS&> (*ds);
  unsigned int ndof = lldds.dimension();
  // Allocate work buffers for:
  // - Iteration matrix
  _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(ndof, ndof, Siconos::BANDED, 0, 0));

  // - I - theta
  work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA].reset(new SiconosVector(ndof));
  // - dt * sigma*
  work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR].reset(new SiconosVector(ndof));
  SimpleMatrix& iteration_matrix = *_dynamicalSystemsGraph->properties(dsv).W;
  // Todo : check perf of two different methods to compute W and friends
  _compute_osi_parameters_v1(lldds, work, iteration_matrix);
}

// void MoreauJeanBilbaoOSI::_compute_osi_parameters_v0(LagrangianLinearDiagonalDS& lldds, VectorOfVectors& work, SimpleMatrix& iteration_matrix_w)
// {
//   // notations:
//   // - stiffness(k) = omega
//   // - damping(k) = sigma
//   SiconosVector& stiffness = *lldds.stiffness();
//   SiconosVector & damping = *lldds.damping();
//   struct op { double operator() (double d) const { return std::exp(d); } };

//   double time_step = _simulation->timeStep();

//   DenseVect& omega2 = *stiffness.dense();
//   DenseVect& sigma = *damping.dense();
//   DenseVect& buffer = *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]->dense();
//   DenseVect& sigma_star = *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR]->dense();
//   // sigma**2, saved in sigma_star. Remind that damping == 2.*sigma
//   std::transform(sigma.begin(), sigma.end(),
//                 sigma.begin(), sigma_star.begin(), std::multiplies<double>() );
//   *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR] *= (0.25);
//   // sqrt(sigma**2 - omega**2) saved in buffer
//   sub(*work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR], omega2, *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]);
//   std::transform(buffer.begin(), buffer.end(), buffer.begin(),
//                  static_cast<std::complex<double> (*)(std::complex<double>)>(std::sqrt));
//   *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA] *= time_step;
//   // exp( work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]) = exp(sqrt(sigma**2 - omega**2), in sigma_star
//   std::transform(buffer.begin(), buffer.end(), sigma_star.begin(), exp);
//   // exp(-work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]) = exp(-sqrt(sigma**2 - omega**2) in buffer
//   *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA] *= -1.;
//   std::transform(buffer.begin(), buffer.end(), buffer.begin(), exp);
//   // sum of exp in buffer
//   add(*work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA], *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR], *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]);
//   // exp(-sigma*time_step) in work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR]
//   *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR] = -time_step * sigma;
//   std::transform(sigma_star.begin(), sigma_star.end(), sigma_star.begin(), exp);
//   // Ak = exp(-sigma*time_step)*(exp(-sqrt(sigma**2 - omega**2)) + exp(sqrt(sigma**2 - omega**2))
//   // i.e. Ak = sigma_star * buffer (component wise), result in buffer
//   std::transform(sigma_star.begin(), sigma_star.end(),
//                 buffer.begin(), buffer.begin(), std::multiplies<double>() );
//   // ek = exp(-2sigma*time_step) = sigma_star**2 in sigma_star
//   std::transform(sigma_star.begin(), sigma_star.end(),
//                 sigma_star.begin(), sigma_star.begin(), std::multiplies<double>() );
//   // 1 + ek in iteration_mat (temp buffer in VFREE)
//   work[MoreauJeanBilbaoOSI::VFREE]->fill(1.);
//   DenseVect& iteration_matrix = *work[MoreauJeanBilbaoOSI::VFREE]->dense();
//   add(*work[MoreauJeanBilbaoOSI::VFREE], *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR],
//       *work[MoreauJeanBilbaoOSI::VFREE]);
//   // 1 - ek in ek (sigma_star)
//   *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR] *= -1.;
//   transform(sigma_star.begin(), sigma_star.end(), sigma_star.begin(),
//             bind2nd(std::plus<double>(), 1.0));
//   // 1 - ek/1 +ek in sigma_star
//   std::transform(sigma_star.begin(), sigma_star.end(),
//                 iteration_matrix.begin(), sigma_star.begin(), std::divides<double>() );
//   // 1 + ek - Ak (i.e. iteration_matrix - buffer) in iteration_matrix
//   sub(*work[MoreauJeanBilbaoOSI::VFREE], *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA],
//       *work[MoreauJeanBilbaoOSI::VFREE]);
//   // Ak / (1 + ek - Ak) in buffer (i.e. buffer/iteration_matrix in buffer)
//   std::transform(buffer.begin(), buffer.end(),
//                 iteration_matrix.begin(), buffer.begin(), std::divides<double>() );
//   // 2./(omega_k**2 * time_step ** 2) in iteration_matrix
//   double coeff = 2. / (time_step * time_step);
//   work[MoreauJeanBilbaoOSI::VFREE]->fill(coeff);
//   std::transform(iteration_matrix.begin(), iteration_matrix.end(),
//                 omega2.begin(), iteration_matrix.begin(), std::divides<double>() );
//   std::transform(iteration_matrix.begin(), iteration_matrix.end(),
//                 omega2.begin(), iteration_matrix.begin(), std::divides<double>() );
//   // 1 - theta_k = 1 - 2/w**2*dt**2 + (Ak / (1+ek -Ak) (i.e. 1. - iteration_matrix + buffer) in buffer
//   sub(*work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA], *work[MoreauJeanBilbaoOSI::VFREE],
//       *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]);
//   transform(buffer.begin(), buffer.end(), buffer.begin(),
//             bind2nd(std::plus<double>(), 1.0));
//   // sigma* * dt= (1 + (1-theta_k)*omega**2*dt**2/2) * (1-ek)/(1+ek) in sigma_star
//   // ie sigma_star = (1 + buffer / iteration_matrix ) * sigma_star
//   std::transform(buffer.begin(), buffer.end(),
//                  iteration_matrix.begin(), iteration_matrix.begin(), std::divides<double>() );
//   transform(iteration_matrix.begin(), iteration_matrix.end(), iteration_matrix.begin(),
//             bind2nd(std::plus<double>(), 1.0));
//   std::transform(iteration_matrix.begin(), iteration_matrix.end(),
//                  sigma_star.begin(), sigma_star.begin(), std::multiplies<double>() );
//   // SUMMARY:
//   // sigma_star = time_step * sigma*
//   // buffer = 1 - theta

//   // now we can compute the iteration matrix ...
//   // W = M + K(1-theta)*dt**2/2 + dt*sigma^*
//   // i.e W = mass + omega*buffer*dt**2*0.5 + sigma_star
//   std::transform(omega2.begin(), omega2.end(),
//                  buffer.begin(), iteration_matrix.begin(), std::multiplies<double>() );
//   SiconosMatrix& mass = *lldds.mass();
//   for(unsigned int i=0;i<lldds.dimension();++i)
//   {
//     iteration_matrix_w(i, i) =  mass(i, i) + 0.5 * time_step * time_step * (*work[MoreauJeanBilbaoOSI::VFREE])(i);
//     iteration_matrix_w(i, i) += (*work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR])(i);
//     iteration_matrix_w(i, i) = 1. / iteration_matrix_w(i, i);
//   }
//   // save 2. * dt * sigma_* in sigma_star
//   double two = 2.;
//   *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR] *= two;
// }

void MoreauJeanBilbaoOSI::_compute_osi_parameters_v1(LagrangianLinearDiagonalDS& lldds, VectorOfVectors& work, SimpleMatrix& iteration_matrix)
{
  SP::SiconosVector omega2 = lldds.stiffness();
  SP::SiconosVector damp = lldds.damping();
  SP::SiconosMatrix mass = lldds.mass();

  double one_minus_theta, dt_sigma_star;
  double time_step = _simulation->timeStep();
  double coeff = 0.5 * time_step * time_step;
  unsigned int ndof = lldds.dimension();
  double one = 1.;
  double two = 2.;
  double omega2k, sigmak, massk;
  for(unsigned int k=0;k<ndof;++k)
  {
    massk = 1.;
    sigmak = 0.;
    omega2k = 0.;
    if(mass)
      massk = (*mass)(k, k);
    if(damp)
      sigmak = 0.5 * (*damp)(k);
    if(omega2)
      omega2k = (*omega2)(k);
    compute_parameters(time_step, omega2k, sigmak, one_minus_theta, dt_sigma_star);
    iteration_matrix(k, k) = one / (massk + coeff * one_minus_theta * omega2k + dt_sigma_star);
    (*work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA])(k) = one_minus_theta;
    (*work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR])(k) = two * dt_sigma_star;
  }
}

void MoreauJeanBilbaoOSI::compute_parameters(double time_step, double omega2, double sigma, double& one_minus_theta, double& dt_sigma_star)
{
  // Computes:
  // 1 - theta
  // sigma_star = dt * sigma^*
  double ek = std::exp(-sigma * time_step);
  std::complex<double> buff = std::sqrt(std::complex<double>(sigma*sigma -omega2)) * time_step;
  std::complex<double> cAk = ek * (std::exp(buff) + std::exp(-buff));
  assert(std::imag(cAk) < 1e-12);
  double Ak = std::real(cAk);
  double one = 1.;
  double two = 2.;
  double one_over_two = 0.5;
  ek = ek * ek;
  double res = omega2 * time_step * time_step;
  one_minus_theta = one - two / res + Ak / (one + ek - Ak);
  dt_sigma_star = (one - ek) / (one + ek) * (one + one_over_two * res * one_minus_theta);
}

// void MoreauJeanBilbaoOSI::_compute_osi_parameters_v2(LagrangianLinearDiagonalDS& lldds, VectorOfVectors& work, SimpleMatrix& iteration_matrix_w)
// {
//   // notations:
//   // - stiffness(k) = omega
//   // - damping(k) = sigma
//   SiconosVector& stiffness = *lldds.stiffness(); // equal to omega**2
//   SiconosVector & damping = *lldds.damping(); // equal to 2*sigma
//   struct op { double operator() (double d) const { return std::exp(d); } };
//   unsigned int ndof = lldds.dimension();
//   double time_step = _simulation->timeStep();

//   DenseVect& omega2 = *stiffness.dense();
//   DenseVect& sigma = *damping.dense();
//   DenseVect& buffer = *work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]->dense();
//   DenseVect& sigma_star = *work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR]->dense();
//   vector<std::complex<double> > buffer_cplx(ndof);
//   // sigma**2, saved in buffer
//   noalias(buffer) = element_prod(sigma, sigma);
//   buffer_cplx = 0.25 * buffer;
//   // dt*sqrt(sigma**2 - omega**2) saved in buffer
//   buffer_cplx -= omega2;
//   std::transform(buffer_cplx.begin(), buffer_cplx.end(), buffer_cplx.begin(),
//                  (std::complex<double>(*)(std::complex<double>)) sqrt);
//   buffer_cplx *= time_step;
//   // exp( work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]) = exp(sqrt(sigma**2 - omega**2), in sigma_star
//   std::transform(buffer_cplx.begin(), buffer_cplx.end(), sigma_star.begin(), exp);
//   // exp(-work[MoreauJeanBilbaoOSI::ONE_MINUS_THETA]) = exp(-sqrt(sigma**2 - omega**2) in buffer
//   buffer_cplx *= -1.;
//   std::transform(buffer_cplx.begin(), buffer_cplx.end(), buffer_cplx.begin(), exp);
//   // sum of exp in buffer
//   buffer_cplx += sigma_star;
//   // exp(-sigma*time_step) in work[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR]
//   sigma_star = -time_step * sigma;
//   std::transform(sigma_star.begin(), sigma_star.end(), sigma_star.begin(), exp);
//   // Ak = exp(-sigma*time_step)*(exp(-sqrt(sigma**2 - omega**2)) + exp(sqrt(sigma**2 - omega**2))
//   // i.e. Ak = sigma_star * buffer (component wise), result in buffer
//   buffer = element_prod(sigma_star, buffer);
//   // ek = exp(-2sigma*time_step) = sigma_star**2 in sigma_star
//   sigma_star = element_prod(sigma_star, sigma_star);
//   // 1 + ek in iteration_mat
//   DenseVect& iteration_matrix = *work[MoreauJeanBilbaoOSI::VFREE]->dense();
//   work[MoreauJeanBilbaoOSI::VFREE]->fill(1.);
//   iteration_matrix += sigma_star;
//   // 1 - ek in ek (sigma_star)
//   sigma_star *= -1.;
//   transform(sigma_star.begin(), sigma_star.end(), sigma_star.begin(),
//             bind2nd(std::plus<double>(), 1.0));
//   // 1 - ek/1 +ek in sigma_star
//   sigma_star = element_div(sigma_star, iteration_matrix);
//   // 1 + ek - Ak (i.e. iteration_matrix - buffer) in iteration_matrix
//   iteration_matrix -= buffer;
//   // Ak / (1 + ek - Ak) in buffer (i.e. buffer/iteration_matrix in buffer)
//   buffer = element_div(buffer, iteration_matrix);
//   // 2./(omega_k**2 * time_step ** 2) in iteration_matrix
//   double coeff = 2. / (time_step * time_step);
//   work[MoreauJeanBilbaoOSI::VFREE]->fill(coeff);
//   iteration_matrix = element_div(iteration_matrix, omega2);
//   // 1 - theta_k = 1 - 2/w**2*dt**2 + (Ak / (1+ek -Ak) (i.e. 1. - iteration_matrix + buffer) in buffer
//   buffer -= iteration_matrix;
//   transform(buffer.begin(), buffer.end(), buffer.begin(),
//             bind2nd(std::plus<double>(), 1.0));
//   // sigma* * dt= (1 + (1-theta_k)*omega**2*dt**2/2) * (1-ek)/(1+ek) in sigma_star
//   // ie sigma_star = (1 + buffer / iteration_matrix ) * sigma_star
//   iteration_matrix = element_div(buffer, iteration_matrix);
//   transform(iteration_matrix.begin(), iteration_matrix.end(), iteration_matrix.begin(),
//             bind2nd(std::plus<double>(), 1.0));
//   sigma_star = element_prod(iteration_matrix, sigma_star);
//   // SUMMARY:
//   // sigma_star = time_step * sigma*
//   // buffer = 1 - theta

//   // now we can compute the iteration matrix ...
//   // W = M + K(1-theta)*dt**2/2 + dt*sigma^*
//   // i.e W = mass + omega**2*buffer*dt**2*0.5 + sigma_star
//   noalias(iteration_matrix) = element_prod(omega2, buffer);
//   SiconosMatrix& mass = *lldds.mass();
//   for(unsigned int i=0;i<lldds.dimension();++i)
//   {
//     iteration_matrix_w(i, i) = mass(i, i) + iteration_matrix(i) * 0.5 * time_step * time_step;
//     iteration_matrix_w(i, i) += sigma_star(i);
//     iteration_matrix_w(i, i) = 1. / iteration_matrix_w(i, i);
//   }
//   // save 2. * dt * sigma_* in sigma_star
//   sigma_star *= 2.;
// }

double MoreauJeanBilbaoOSI::computeResidu()
{
  return 0.;
}

void MoreauJeanBilbaoOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  double time_step = _simulation->timeStep();
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    // Nothing to be done if the osi is not linked to the ds
    if(!checkOSI(dsi)) continue;
    //
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    LagrangianLinearDiagonalDS& lldds = static_cast<LagrangianLinearDiagonalDS&> (*ds);
    VectorOfVectors& work_ds = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // Get velocity computed at the beginning of the time step.
    SiconosVector& v_i = *lldds.velocityMemory()->getSiconosVector(0);
    SiconosVector& q_i = *lldds.qMemory()->getSiconosVector(0);
    SiconosVector& stiffness = *lldds.stiffness();
    // Get iteration matrix
    const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
    SimpleMatrix& inv_iteration_matrix = *_dynamicalSystemsGraph->properties(dsv).W;
    // Get 2.*dt*sigma^*
    SiconosVector& two_dt_sigma_star = *work_ds[MoreauJeanBilbaoOSI::TWO_DT_SIGMA_STAR];
    // buffer for vfree
    SiconosVector& vfree = *work_ds[MoreauJeanBilbaoOSI::VFREE];
    // Compute vfree
    unsigned int dimension = lldds.dimension();
    for(unsigned int k=0; k<dimension;++k)
      vfree(k) = v_i(k) - inv_iteration_matrix(k, k) * (time_step * stiffness(k) * q_i(k) + two_dt_sigma_star(k) * v_i(k));
  }

}

struct MoreauJeanBilbaoOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  Interaction& _inter;
  InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(OneStepNSProblem *p, Interaction& inter, InteractionProperties& interProp) :
    _osnsp(p), _inter(inter), _interProp(interProp) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter.nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[MoreauJeanBilbaoOSI::OSNSP_RHS];
    subscal(e, *_inter.y_k(_osnsp->inputOutputLevel()), osnsp_rhs, subCoord, false);
  }

};

void MoreauJeanBilbaoOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  InteractionsGraph& indexSet = *osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  assert(indexSet.bundle(vertex_inter));
  VectorOfBlockVectors& DSlink = *indexSet.properties(vertex_inter).DSlink;
  // current interaction
  Interaction& inter = *indexSet.bundle(vertex_inter);
  assert(inter.relation());
  // Get relation and non smooth law types
  // RELATION::TYPES relationType = inter.relation()->getType();
  // RELATION::SUBTYPES relationSubType = inter.relation()->getSubType();
  // check relation type: done in fillDSLinks.
  // if(relationType != Lagrangian || relationSubType != LinearR)
  //   RuntimeException::selfThrow("MoreauJeanBilbaoOSI::computeFreeOutput only Lagrangian Linear Relations are allowed.");

  unsigned int sizeY = inter.nonSmoothLaw()->size();
  unsigned int relativePosition = 0;

  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  // buffer used to save output
  BlockVector& x_free = *DSlink[LagrangianR::xfree];
  SiconosVector& osnsp_rhs = *(*indexSet.properties(vertex_inter).workVectors)[MoreauJeanBilbaoOSI::OSNSP_RHS];


  if(inter.relation()->C())
  {
    SiconosMatrix&  C = *inter.relation()->C() ;
    coord[3] = C.size(1);
    coord[5] = C.size(1);
    // osnsp_rhs[coord] = C.x_free
    subprod(C, x_free, osnsp_rhs, coord, true);
  }
  _NSLEffectOnFreeOutput nslEffectOnFreeOutput = _NSLEffectOnFreeOutput(osnsp, inter, indexSet.properties(vertex_inter));
  inter.nonSmoothLaw()->accept(nslEffectOnFreeOutput);

}

void MoreauJeanBilbaoOSI::integrate(double& tinit, double& tend, double& tout, int& notUsed)
{
  RuntimeException::selfThrow("MoreauJeanBilbaoOSI::integrate - Not yet implemented for MoreauJeanBilbaoOSI.");
}

void MoreauJeanBilbaoOSI::updatePosition(DynamicalSystem& ds)
{
  // --  Update current position for the last computed velocities --
  // q(i+1) = q(i) + dt v(i+1)
  double time_step = _simulation->timeStep();
  // get dynamical system
  LagrangianLinearDiagonalDS& d = static_cast<LagrangianLinearDiagonalDS&> (ds);
  // Compute q
  SiconosVector& v = *d.velocity();
  SiconosVector& q = *d.q();
  //  -> get positions at the beginning of the time step
  SiconosVector& qold = *d.qMemory()->getSiconosVector(0);
  // update positions
  scal(time_step, v, q, true) ;
  q += qold;
}


void MoreauJeanBilbaoOSI::updateState(const unsigned int )
{
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if(useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SiconosMatrix& inv_iteration_matrix = *_dynamicalSystemsGraph->properties(*dsi).W;
    // get dynamical system and work vector
    LagrangianLinearDiagonalDS& lldds = static_cast<LagrangianLinearDiagonalDS&> (*_dynamicalSystemsGraph->bundle(*dsi));
    VectorOfVectors& work_ds = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    SiconosVector& vfree = *work_ds[MoreauJeanBilbaoOSI::VFREE];

    SiconosVector &v = *lldds.velocity();
    if(lldds.p(_levelMaxForInput) && lldds.p(_levelMaxForInput)->size() > 0)
    {
      v = *lldds.p(_levelMaxForInput); // v = p
      if(lldds.boundaryConditions())
        for(std::vector<unsigned int>::iterator
              itindex = lldds.boundaryConditions()->velocityIndices()->begin() ;
            itindex != lldds.boundaryConditions()->velocityIndices()->end();
            ++itindex)
          v.setValue(*itindex, 0.0);
      unsigned int ndof = lldds.dimension();
      for(unsigned int k=0;k<ndof;++k)
        v(k) = vfree(k) + v(k) * inv_iteration_matrix(k, k);
    }
    else
      v =  vfree;
    // Update positions with the last computed velocities.
    updatePosition(lldds);
  }
}

void MoreauJeanBilbaoOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== MoreauJeanBilbaoOSI OSI display ======" <<std::endl;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  if(_dynamicalSystemsGraph)
  {
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

      std::cout << "--------------------------------" <<std::endl;
      std::cout << "--> W of dynamical system number " << ds->number() << ": " <<std::endl;
      if(_dynamicalSystemsGraph->properties(*dsi).W) _dynamicalSystemsGraph->properties(*dsi).W->display();
      else std::cout << "-> NULL" <<std::endl;
    }
  }
  std::cout << "================================" <<std::endl;
}

void MoreauJeanBilbaoOSI::prepareNewtonIteration(double time)
{
}


bool MoreauJeanBilbaoOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  DEBUG_PRINT("addInteractionInIndexSet(SP::Interaction inter, unsigned int i)\n");

  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity

  double gamma = 1.0 / 2.0;
  y += gamma * h * yDot;
  assert(!isnan(y));
  DEBUG_EXPR(
    if(y <= 0)
      DEBUG_PRINT("MoreauJeanBilbaoOSI::addInteractionInIndexSet ACTIVATE.\n");
    );
  return (y <= 0.0);
}


bool MoreauJeanBilbaoOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  DEBUG_PRINTF("MoreauJeanBilbaoOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!isnan(y));

  DEBUG_EXPR(
    if(y > 0)
      DEBUG_PRINT("MoreauJeanOSI::removeInteractionInIndexSet DEACTIVATE.\n");
    );
  return (y > 0.0);
}



