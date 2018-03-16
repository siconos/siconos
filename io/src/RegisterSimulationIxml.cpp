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

#include "RegisterSimulation.hpp"

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_sparse.hpp>
#include <boost/numeric/bindings/ublas/matrix_sparse.hpp>

#include <boost/archive/xml_iarchive.hpp>

// Work-around for issue reading inf/nan double values
#include <boost/archive/basic_text_iprimitive.hpp>
namespace boost { namespace archive {
template<> template<>
void basic_text_iprimitive<std::istream>::load<double>( double& t )
{
  char s[32]={}, *b, *a = b = s;
  int i=0;
  for (int i=0; is.peek()!='<' && i < 31; i++)
    s[i] = is.get();
  errno = 0;
  t = std::strtod(s, &b);
  if (errno == ERANGE || errno == EINVAL || a==b)
    boost::serialization::throw_exception(
      archive_exception(archive_exception::input_stream_error));
}
} } // close namespaces

void RegisterSimulationIxml(std::ifstream& ifs, SP::Simulation& sim)
{
  boost::archive::xml_iarchive ar(ifs);
  siconos_io_register_Numerics(ar);
  siconos_io_register_Kernel(ar);
  siconos_io_register_Mechanics(ar);
  siconos_io_register_Control(ar);
  ar >> NVP(sim);
}
#endif
