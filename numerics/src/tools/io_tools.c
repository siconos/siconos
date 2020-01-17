/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
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


#include "io_tools.h"
#include <string.h>  // for strcmp, strrchr

int check_hdf5_file(const char * filename)
{
  // get last occurence of dot in the file name
  // and shift one char forward to get the extension
  const char * ext = strrchr(filename, '.') + 1;
  // compare ext with standard ext for hdf5 files.
  int res = strcmp(ext, "hdf5");
  if(res != 0)
    res = strcmp(ext+1, "h5");

  return !res;
}


