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
class OccBody;

#include <SiconosVisitables.hpp>
#undef SICONOS_VISITABLES

#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()                              \
  REGISTER(OccBody)



#include "OccTimeStepping.hpp"
#include "OccBody.hpp"

#include <NonSmoothDynamicalSystem.hpp>

#include <SiconosVisitor.hpp>

#define VISITOR_CLASSES()                       \
  REGISTER(OccBody)

#include <VisitorMaker.hpp>



using namespace Experimental;

struct UpdateShapes : public SiconosVisitor
{
  using SiconosVisitor::visit;

  template<typename T>
  void operator()(const T& ds)
  {
    const_cast<T&>(ds).updateShapes();
    const_cast<T&>(ds).updateContactShapes();
  }
};


void OccTimeStepping::updateWorldFromDS()
{
  DynamicalSystemsGraph& dsg = *_nsds->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg.vertices();

  Visitor< Classes < OccBody >, UpdateShapes >::Make up;

  for(; dsi != dsiend; ++dsi)
  {
    dsg.bundle(*dsi)->accept(up);
  }

}
