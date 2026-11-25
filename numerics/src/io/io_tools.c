/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "io_tools.h"

#include <string.h>  // for strcmp, strrchr

#include "numerics_errors.h"

int check_hdf5_file(const char* filename) {
  // get last occurence of dot in the file name
  // and shift one char forward to get the extension
  const char* ext = strrchr(filename, '.') + 1;
  // compare ext with standard ext for hdf5 files.
  int res = strcmp(ext, "hdf5");
  if (res != 0) res = strcmp(ext + 1, "h5");

  return !res;
}
