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

#ifndef SiconosFullNumerics_hpp
#define SiconosFullNumerics_hpp

#include "SiconosConfig.h"
#ifdef WITH_SERIALIZATION
#include "Register.hpp"



template <class Archive>
void siconos_io(Archive& ar, Callback&v, unsigned int version)
{
}
REGISTER_BOOST_SERIALIZATION(Callback);

template <class Archive>
void siconos_io(Archive& ar, SolverOptions&v, unsigned int version)
{
  SERIALIZE(v, (solverId)(isSet)(iSize)(dSize)(filterOn)(numberOfInternalSolvers), ar);

  if (Archive::is_loading::value)
  {
    solver_options_nullify(&v);
    v.iparam = (int *) malloc(v.iSize * sizeof(int));
    v.dparam = (double *) malloc(v.dSize * sizeof(double));
    v.internalSolvers = (SolverOptions *) malloc(v.numberOfInternalSolvers * sizeof(SolverOptions));
    v.callback = (Callback *) malloc(sizeof(Callback));
  }
  SERIALIZE(v, (callback), ar);

  SERIALIZE_C_ARRAY(v.iSize, v, iparam, ar);
  SERIALIZE_C_ARRAY(v.dSize, v, dparam, ar);
  SERIALIZE_C_ARRAY(v.numberOfInternalSolvers, v, internalSolvers, ar);
}
REGISTER_BOOST_SERIALIZATION(SolverOptions);

template <class Archive>
void siconos_io(Archive& ar, LinearComplementarityProblem& v, unsigned int version)
{
  SERIALIZE(v, (size), ar);

  if(Archive::is_loading::value)
  {
    v.q = (double *) malloc(v.size * sizeof(double));
    v.M = NM_new();
  }
  SERIALIZE(v, (M), ar);
  SERIALIZE_C_ARRAY(v.size, v, q, ar);
}
REGISTER_BOOST_SERIALIZATION(LinearComplementarityProblem);

template <class Archive>
void siconos_io(Archive& ar, FrictionContactProblem& p, const unsigned int file_version)
{
  SERIALIZE(p, (dimension)(numberOfContacts), ar);

  if (Archive::is_loading::value)
  {
    p.q = (double *) malloc(p.dimension * p.numberOfContacts * sizeof(double));
    p.mu = (double *) malloc(p.numberOfContacts * sizeof(double));
    p.M = NM_new();
  }

  SERIALIZE(p, (M), ar);
  SERIALIZE_C_ARRAY(p.dimension * p.numberOfContacts, p, q, ar);
  SERIALIZE_C_ARRAY(p.dimension, p, mu, ar);
}
REGISTER_BOOST_SERIALIZATION(FrictionContactProblem);

template <class Archive>
void siconos_io(Archive& ar, SparseBlockStructuredMatrix& v, unsigned int version)
{
  SERIALIZE(v, (nbblocks)(blocknumber0)(blocknumber1)(filled1)(filled2), ar);
  if (Archive::is_loading::value)
  {
    v.block = (double **) malloc(v.nbblocks * sizeof(double *));
    v.blocksize1 = (unsigned int *) malloc (v.blocknumber1* sizeof(unsigned int));
    v.blocksize0 = (unsigned int *) malloc (v.blocknumber0* sizeof(unsigned int));
    SERIALIZE_C_ARRAY(v.blocknumber1, v, blocksize1, ar);
    SERIALIZE_C_ARRAY(v.blocknumber0, v, blocksize0, ar);
    int diagonalblocknumber  = v.blocknumber1 +
      ((v.blocknumber0 - v.blocknumber1) & -(v.blocknumber0 < v.blocknumber1));
    for (unsigned int i=0; i< diagonalblocknumber; ++i)
    {
      unsigned int size0 = v.blocksize0[i];
      if (i != 0) size0 -= v.blocksize0[i - 1];
      unsigned int size1 = v.blocksize1[i];
      if (i != 0) size1 -= v.blocksize1[i - 1];
      v.block[i] = (double*) malloc(size0 * size1 * sizeof(double));
    }
    v.index1_data = (size_t *) malloc (v.filled1 * sizeof(size_t));
    v.index2_data = (size_t *) malloc (v.filled2 * sizeof(size_t));
  }
  else
  {
    SERIALIZE_C_ARRAY(v.blocknumber1, v, blocksize1, ar);
    SERIALIZE_C_ARRAY(v.blocknumber0, v, blocksize0, ar);
  }

  int diagonalblocknumber  = v.blocknumber1 +
      ((v.blocknumber0 - v.blocknumber1) & -(v.blocknumber0 < v.blocknumber1));

  for (unsigned int i=0; i< v.nbblocks; ++i)
  {
    ar & ::boost::serialization::make_nvp("block", (long&) v.block[i]);
  }

  for (unsigned int i=0; i< diagonalblocknumber; ++i)
  {
    unsigned int size0 = v.blocksize0[i];
    if (i != 0) size0 -= v.blocksize0[i - 1];
    unsigned int size1 = v.blocksize1[i];
    if (i != 0) size1 -= v.blocksize1[i - 1];
    for (unsigned int k=0; k<size0 * size1; ++k)
    {
      ar & ::boost::serialization::make_nvp("item", v.block[i][k]);
    }
  }

  SERIALIZE_C_ARRAY(v.filled1,  v, index1_data, ar);
  SERIALIZE_C_ARRAY(v.filled2,  v, index2_data, ar);
}
REGISTER_BOOST_SERIALIZATION(SparseBlockStructuredMatrix);

template <class Archive>
void siconos_io(Archive&ar, NumericsMatrix& v, unsigned int version)
{
  SERIALIZE(v, (storageType)(size0)(size1), ar);
  if (v.storageType == 0)
  {
    if (Archive::is_loading::value)
    {
      v.matrix0 = (double *) malloc(v.size0 * v.size1 * sizeof(double));
      v.matrix1 = NULL;
      v.matrix2 = NULL;
      v.internalData = NULL;
    }
    SERIALIZE_C_ARRAY(v.size0 * v.size1, v, matrix0, ar);
  }
  else
  {
    {
      if (Archive::is_loading::value)
      {
        v.matrix0 = NULL;
        v.matrix1 = (SparseBlockStructuredMatrix*) malloc(sizeof(SparseBlockStructuredMatrix));
        v.matrix2 = NULL;
        v.internalData = NULL;
      }
      SERIALIZE(v, (matrix1), ar);
    }
  }
}
REGISTER_BOOST_SERIALIZATION(NumericsMatrix);

template <class Archive>
void siconos_io_register_Numerics(Archive& ar)
{
  ar.register_type(static_cast<SolverOptions*>(NULL));
  ar.register_type(static_cast<LinearComplementarityProblem*>(NULL));
  ar.register_type(static_cast<NumericsMatrix*>(NULL));
  ar.register_type(static_cast<SparseBlockStructuredMatrix*>(NULL));
  ar.register_type(static_cast<FrictionContactProblem*>(NULL));

}
#endif
#endif
