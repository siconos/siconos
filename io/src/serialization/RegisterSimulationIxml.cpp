// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifdef WITH_SERIALIZATION
#include <boost/archive/xml_iarchive.hpp>

#include "RegisterSimulation.hpp"
#include "SiconosFull.hpp"

// Work-around for issue reading inf/nan double values
#include <boost/archive/basic_text_iprimitive.hpp>
namespace boost {
namespace archive {
template <>
template <>
void basic_text_iprimitive<std::istream>::load<double>(double& t) {
  char s[32] = {}, *b, *a = b = s;
  int i = 0;
  for (int i = 0; is.peek() != '<' && i < 31; i++) s[i] = is.get();
  errno = 0;
  t = std::strtod(s, &b);
  if (errno == ERANGE || errno == EINVAL || a == b)
    boost::serialization::throw_exception(
        archive_exception(archive_exception::input_stream_error));
}
}  // namespace archive
}  // namespace boost

void RegisterSimulationIxml(std::ifstream& ifs,
                            std::shared_ptr<siconos::simulation::Simulation>& sim) {
  boost::archive::xml_iarchive ar(ifs);
  siconos_io_register_Numerics(ar);
  siconos_io_register_Kernel(ar);
  siconos_io_register_Mechanics(ar);
  siconos_io_register_Control(ar);
  ar >> NVP(sim);
}
#endif
