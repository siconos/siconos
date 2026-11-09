// The one place where the disks world is constructed, hence where its columns
// are registered: kept in its own translation unit, compiled in parallel with
// the binding parts (see data_t in py_head.hpp).
#include "py_head.hpp"

namespace siconos::python::disks {
data_t::data_t() : _data(idata_t{}) {}
}  // namespace siconos::python::disks
