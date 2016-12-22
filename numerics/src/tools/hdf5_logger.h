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





#include <stdbool.h>
#include <stdint.h>

#include "NumericsMatrix.h"

#ifdef WITH_HDF5
#include <hdf5.h>

#else

#define hid_t void*

#endif /* WITH_HDF5 */

/** \struct SN_logh5
 HDF5 logger (for full debug purposes) */
typedef struct {
  char* itername; /**< name of the group for the current iteration */
  unsigned itername_len; /**< maximum length of itername */
  hid_t file; /**< handle to the HDF5 file */
  hid_t group; /**< Handle to the current (top-level) group (e.g. "/foo") */
} SN_logh5;


bool SN_logh5_check_gzip(void);

SN_logh5*  SN_logh5_init(const char* filename, const unsigned iter_max);

bool SN_logh5_end(SN_logh5* logger);

bool SN_logh5_new_iter(unsigned iter, SN_logh5* logger);

bool SN_logh5_end_iter(SN_logh5* logger);

bool SN_logh5_scalar_double(double val, const char* name, hid_t loc_id);

bool SN_logh5_scalar_integer(ptrdiff_t val, const char* name, hid_t loc_id);

bool SN_logh5_scalar_uinteger(size_t val, const char* name, hid_t loc_id);

bool SN_logh5_attr_uinteger(size_t val, const char* name, hid_t loc_id);

bool SN_logh5_vec_double(size_t size, double* vec, const char* name, hid_t loc_id);

bool SN_logh5_NM(NumericsMatrix* mat, const char* name, SN_logh5* logger);

bool SN_logh5_mat_dense(size_t size0, size_t size1, double* mat, const char* name, hid_t loc_id);

bool SN_logh5_csparse(CSparseMatrix* cs, const char* name, hid_t loc_id);

bool SN_logh5_vec_int32(size_t size, int32_t* vec, const char* name, hid_t loc_id);

bool SN_logh5_vec_int64(size_t size, int64_t* vec, const char* name, hid_t loc_id);



