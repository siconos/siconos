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

#include <stdio.h>
#include <math.h>

#include "SparseMatrix_internal.h"
#include "hdf5_logger.h"
#include "NumericsSparseMatrix.h"

#ifdef WITH_HDF5

bool SN_logh5_check_gzip(void)
{

  unsigned filter_info;

  htri_t avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE);
  if (!avail) {
    printf ("gzip filter not available.\n");
    return false;
  }

  herr_t status = H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
  if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
      !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
    printf ("gzip filter not available for encoding and decoding.\n");
    return false;
  }

  return true;

}

SN_logh5* SN_logh5_init(const char* filename, const unsigned iter_max)
{
  herr_t status;

  assert(filename);
  assert(iter_max > 0);

  SN_logh5* logger = (SN_logh5*) malloc(sizeof(SN_logh5));
  if (!logger)
  {
    fprintf(stderr, "SN_logh5_init :: could not allocate a logger struct!\n");
    exit(EXIT_FAILURE);
  }

  logger->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  /* formula : 2 ("i-") + 1 (end '\0') + number of char to represent the number
   * of iteration*/
  logger->itername_len = 4 + (unsigned)ceil(log10(iter_max));

  logger->itername = (char*)calloc(logger->itername_len, sizeof(char));

  logger->group = 0;

  return logger;
}

bool SN_logh5_end(SN_logh5* logger)
{
  assert(logger);
  herr_t status = 0;

  if (logger->group)
  {
    fprintf(stderr, "SN_logh5_end :: group pointer was not 0! Most likely the last opened group was not properly closed\n");
    SN_logh5_end_iter(logger);
  }
  status = H5Fclose(logger->file);
  logger->file = 0;

  free(logger->itername);
  logger->itername = NULL;

  free(logger);

  return !status ? true : false;
}

bool SN_logh5_new_iter(unsigned iter, SN_logh5* logger)
{
  assert(logger);
  herr_t status = 0;

  if (logger->group)
  {
    fprintf(stderr, "SN_logh5_new_iter :: group pointer was not 0! Most likely the last opened group was not properly closed\n");
    status = H5Gclose (logger->group);
  }

  snprintf(logger->itername, logger->itername_len, "i-%d", iter);

  logger->group = H5Gcreate (logger->file, logger->itername, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  return !status ? true : false;
}

bool SN_logh5_end_iter(SN_logh5* logger)
{
  assert(logger);
  assert(logger->group);

  if (!logger->group)
  {
    fprintf(stderr, "SN_logh5_end_iter :: group pointer is NULL !!\n");
    return false;
  }

  herr_t status = H5Gclose (logger->group);
  logger->group = 0;

  return status;
}

static bool SN_logh5_write_dset(hid_t loc_id, const char* name, hid_t type, hid_t space, void* val)
{
  herr_t status;
  hid_t dset = H5Dcreate (loc_id, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
  status = H5Dclose (dset);

  return !status ? true : false;
}

static bool SN_logh5_write_attr(hid_t loc_id, const char* name, hid_t type, hid_t space, void* val)
{
  herr_t status;
  hid_t dset = H5Acreate (loc_id, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite (dset, type, val);
  status = H5Aclose (dset);

  return !status ? true : false;
}

bool SN_logh5_scalar_double(double val, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(name);

  herr_t status;
  hid_t space = H5Screate(H5S_SCALAR);

  status = SN_logh5_write_dset(loc_id, name, H5T_NATIVE_DOUBLE, space, &val);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_scalar_integer(ptrdiff_t val, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(name);

  hid_t type;
  herr_t status;

  switch (sizeof(ptrdiff_t))
  {
    case sizeof(int32_t):
    {
      type = H5T_NATIVE_INT_LEAST32;
      break;
    }
    case sizeof(int64_t):
    {
      type = H5T_NATIVE_INT_LEAST64;
      break;
    }
    default:
      fprintf(stderr, "logh5 :: unsupported integer length %lu\n", sizeof(ptrdiff_t));
      return false;
  }

  hid_t space = H5Screate(H5S_SCALAR);
  status = SN_logh5_write_dset(loc_id, name, type, space, &val);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_scalar_uinteger(size_t val, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(name);

  hid_t type;
  herr_t status;

  switch (sizeof(size_t))
  {
    case sizeof(int32_t):
    {
      type = H5T_NATIVE_UINT_LEAST32;
    break;
    }
    case sizeof(int64_t):
    {
      type = H5T_NATIVE_UINT_LEAST64;
      break;
    }
    default:
      fprintf(stderr, "logh5 :: unsupported unsigned integer length %lu\n", sizeof(size_t));
      return false;
  }

  hid_t space = H5Screate(H5S_SCALAR);
  status = SN_logh5_write_dset(loc_id, name, type, space, &val);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_attr_uinteger(size_t val, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(name);

  hid_t type;
  herr_t status;

  switch (sizeof(size_t))
  {
    case sizeof(int32_t):
    {
      type = H5T_NATIVE_UINT_LEAST32;
    break;
    }
    case sizeof(int64_t):
    {
      type = H5T_NATIVE_UINT_LEAST64;
      break;
    }
    default:
      fprintf(stderr, "logh5 :: unsupported unsigned integer length %lu\n", sizeof(size_t));
      return false;
  }

  hid_t space = H5Screate(H5S_SCALAR);
  status = SN_logh5_write_attr(loc_id, name, type, space, &val);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_vec_double(size_t size, double* vec, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(vec);
  assert(name);

  hsize_t dims[1] = {size};
  hid_t type = H5T_NATIVE_DOUBLE;
  herr_t status;

  hid_t space =  H5Screate_simple(1, dims, NULL);
  status = SN_logh5_write_dset(loc_id, name, type, space, vec);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_csparse(CSparseMatrix* cs, const char* name, hid_t loc_id)
{
  assert(cs);
  assert(loc_id);
  hid_t mat_group;

  bool result = false;

  if (name)
  {
    /* Create subgroup here */
    mat_group = H5Gcreate(loc_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  else
  {
    /* We are already in the matrix subgroup */
    mat_group = loc_id;
  }

  result = SN_logh5_scalar_integer(cs->nzmax, "nzmax", mat_group);
  result = SN_logh5_scalar_integer(cs->m, "m", mat_group);
  result = SN_logh5_scalar_integer(cs->n, "n", mat_group);
  result = SN_logh5_scalar_integer(cs->nz, "nz", mat_group);

  if (sizeof(CS_INT) == sizeof(int64_t))
  {
    result = SN_logh5_vec_int64(cs->n+1, cs->p, "p", mat_group);
    result = SN_logh5_vec_int64(cs->p[cs->n], cs->i, "i", mat_group);
  }
/*   else if (sizeof(CS_INT) == sizeof(int32_t))
  {
    result = SN_logh5_vec_int32(cs->n+1, cs->p, "p", mat_group);
    result = SN_logh5_vec_int32(cs->nzmax, cs->i, "i", mat_group);
  }*/
  else
  {
    fprintf(stderr, "SN_logh5_NM :: unknown pointer size %lu\n", sizeof(CS_INT));
    result = false;
  }
  result = SN_logh5_vec_double(cs->p[cs->n], cs->x, "x", mat_group);

  if (name)
  {
    H5Gclose (mat_group);
  }
  return result;
}

bool SN_logh5_NM(NumericsMatrix* mat, const char* name, SN_logh5* logger)
{
  assert(mat);
  assert(name);

  bool result = false;
  herr_t status;

  unsigned storageType = mat->storageType;

  /* Create subgroup based on the name of the matrix*/
  hid_t mat_group;
  if (logger->group)
  {
    mat_group = H5Gcreate(logger->group, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  else
  {
    mat_group = H5Gcreate(logger->file, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  /* Add the storagetype as an attribute  */
  SN_logh5_attr_uinteger(storageType, "storageType", mat_group);
  SN_logh5_attr_uinteger(mat->size0, "size0", mat_group);
  SN_logh5_attr_uinteger(mat->size1, "size1", mat_group);

  switch (storageType)
  {
    case NM_DENSE:
    {
      result = SN_logh5_mat_dense(mat->size0, mat->size1, mat->matrix0, "values", mat_group);
      break;
    }
    case NM_SPARSE_BLOCK:
    {
      if (!(mat->matrix2 &&  mat->matrix2->csc))
      {
        fprintf(stderr, "SN_logh5_NM :: sparse block matrix are currently not supported\n");
        break;
      }
    }
    case NM_SPARSE:
    {
      CSparseMatrix* cs = mat->matrix2->csc;
      assert(cs && "SN_logh5_NM only csc matrix are supported");
      SN_logh5_csparse(cs, NULL, mat_group);
      break;
    }
    default:
    {
      fprintf(stderr, "SN_logh5_NM :: unknown matrix type %d\n", mat->storageType);
      return false;
    }
  }

  status = H5Gclose (mat_group);
  return true;
}

bool SN_logh5_mat_dense(size_t size0, size_t size1, double* mat, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(mat);
  assert(name);

  hsize_t dims[2] = {size0, size1};
  hid_t type = H5T_NATIVE_DOUBLE;
  herr_t status;

  hid_t space =  H5Screate_simple(2, dims, NULL);
  status = SN_logh5_write_dset(loc_id, name, type, space, mat);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_vec_int32(size_t size, int32_t* vec, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(vec);
  assert(name);

  hsize_t dims[1] = {size};
  hid_t type = H5T_NATIVE_INT_LEAST32;
  herr_t status;

  hid_t space =  H5Screate_simple(1, dims, NULL);
  status = SN_logh5_write_dset(loc_id, name, type, space, vec);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_vec_int64(size_t size, int64_t* vec, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(vec);
  assert(name);

  hsize_t dims[1] = {size};
  hid_t type = H5T_NATIVE_INT_LEAST64;
  herr_t status;

  hid_t space =  H5Screate_simple(1, dims, NULL);
  status = SN_logh5_write_dset(loc_id, name, type, space, vec);
  status = H5Sclose(space);

  return !status ? true : false;
}

bool SN_logh5_vec_uint64(size_t size, uint64_t* vec, const char* name, hid_t loc_id)
{
  assert(loc_id);
  assert(vec);
  assert(name);

  hsize_t dims[1] = {size};
  hid_t type = H5T_NATIVE_UINT_LEAST64;
  herr_t status;

  hid_t space =  H5Screate_simple(1, dims, NULL);
  status = SN_logh5_write_dset(loc_id, name, type, space, vec);
  status = H5Sclose(space);

  return !status ? true : false;
}
#else /* WITH_HDF5 */

bool SN_logh5_check_gzip(void)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

SN_logh5* SN_logh5_init(const char* filename, const unsigned iter_max)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return NULL;
}

bool SN_logh5_end(SN_logh5* logger)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_new_iter(unsigned iter, SN_logh5* logger)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_end_iter(SN_logh5* logger)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_scalar_double(double val, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_scalar_integer(ptrdiff_t val, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_scalar_uinteger(size_t val, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_vec_double(size_t size, double* vec, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_NM(NumericsMatrix* mat, const char* name, SN_logh5* logger)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_mat_dense(size_t size0, size_t size1, double* mat, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_vec_int32(size_t size, int32_t* vec, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_vec_int64(size_t size, int64_t* vec, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

bool SN_logh5_vec_uint64(size_t size, uint64_t* vec, const char* name, hid_t loc_id)
{
  fprintf(stderr, "SN_logh5 :: Siconos/Numerics has been compiled with no HDF5 support!\n");
  return false;
}

#endif /*  WITH_HDF5 */
