// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef IO_TOOLS

/*!\file io_tools.h
  \brief utilities to deal with files reading/writing in numerics.
*/

#if defined(__cplusplus)
extern "C" {
#endif

/** return 1 if filename is hdf5

  Warning: we only check the extension (h5 or hdf5)

*/
int check_hdf5_file(const char* filename);

#if defined(__cplusplus)
}
#endif
#endif
