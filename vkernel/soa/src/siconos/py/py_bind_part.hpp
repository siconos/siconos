#pragma once

// Body of one binding part.  Not compiled on its own: CMake generates
// bind_part_<k>.cpp files that define NONOS_PART and include this header, and
// NONOS_PARTS (the number of parts) is a compile definition of _nonos.

#include "py_bind.hpp"

#if !defined(NONOS_PART) || !defined(NONOS_PARTS)
#error "NONOS_PART and NONOS_PARTS must be defined"
#endif

namespace siconos::python::disks {

template <std::size_t Part>
void bind_part(py::module_& disks) {
  if constexpr (Part == 0) {
    disks.def("make_storage", &make_storage, py::return_value_policy::reference,
              R"pbdoc(
        Create a new data object for 2D disks simulation
    )pbdoc");

    py::class_<data_t>(disks, "data_t");
  }

  bind_items<Part, NONOS_PARTS>(disks);
}

// the one part of this translation unit; main.cpp calls them all
template void bind_part<NONOS_PART>(py::module_& disks);

}  // namespace siconos::python::disks
