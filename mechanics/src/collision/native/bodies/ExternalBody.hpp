/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file ExternalBody.hpp
  \brief External visitable lagrangian class.
*/
#ifndef ExternalBody_hpp
#define ExternalBody_hpp

#include <LagrangianDS.hpp>

namespace siconos::collision::native {
class SpaceFilter;
}

namespace siconos::collision::native::bodies {

class ExternalBody : public siconos::modeling::LagrangianDS,
                     public std::enable_shared_from_this<ExternalBody> {
 public:
  virtual void selfHash(siconos::collision::native::SpaceFilter&) = 0;

  virtual void selfFindInteractions(
      std::shared_ptr<siconos::collision::native::SpaceFilter>) = 0;

 protected:
  ACCEPT_SERIALIZATION(ExternalBody);
};
}  // namespace siconos::collision::native::bodies
#endif
