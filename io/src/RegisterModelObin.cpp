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

#include "SiconosConfig.h"
#ifdef WITH_SERIALIZATION
#include "SiconosFull.hpp"

#include "RegisterModel.hpp"

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_sparse.hpp>
#include <boost/numeric/bindings/ublas/matrix_sparse.hpp>

#include <boost/archive/binary_oarchive.hpp>

void RegisterModelObin(std::ofstream& ofs, SP::Model& model)
{
  boost::archive::binary_oarchive ar(ofs);
  siconos_io_register_Numerics(ar);
  siconos_io_register_Kernel(ar);
  siconos_io_register_Mechanics(ar);
  siconos_io_register_Control(ar);
  ar << NVP(model);
}
#endif
