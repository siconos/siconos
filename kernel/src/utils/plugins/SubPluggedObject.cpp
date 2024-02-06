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

#include "SubPluggedObject.hpp"

#include "SimpleMatrix.hpp"

siconos::plugins::SubPluggedObject::SubPluggedObject(const PluggedObject& PO,
                                                     const unsigned int n,
                                                     const unsigned int p,
                                                     const unsigned int indx)
  : PluggedObject{"Sub" + PO.pluginName()}, _indx(indx), _p(p)
{
  _tmpMat = std::make_shared<siconos::algebra::SimpleMatrix>(n, p);
#if (__GNUG__ && !(__clang__ || __INTEL_COMPILER || __APPLE__) && \
     (((__GNUC__ > 5) && (__GNUC_MINOR__ > 0))))
#pragma GCC diagnostic ignored "-Wpmf-conversions"
  fPtr = (void*)&SubPluggedObject::computeAndExtract;
  _parentfPtr = PO.fPtr;
#else
  THROW_EXCEPTION("SubPluggedObject must be compiled with GCC !");
#endif
};

// SubPluggedObject(const SubPluggedObject& SPO)
//   : PluggedObject(SPO), _indx(SPO.getIndex()), _p(SPO.getp())
// {
//   _parentfPtr = SPO.getParentfPtr();
//   _tmpMat = std::make_shared<siconos::algebra::SimpleMatrix>(SPO.getTmpMat());
// }

void siconos::plugins::SubPluggedObject::computeAndExtract(double time, unsigned int n,
                                                           double* M, unsigned int sizez,
                                                           double* z)
{
  ((matrixPlugin)_parentfPtr)(time, n, _p, &(*_tmpMat)(0, 0), sizez, z);
  for (unsigned int i = 0; i < n; i++) {
    M[i] = (*_tmpMat)(i, _indx);
  }
};
