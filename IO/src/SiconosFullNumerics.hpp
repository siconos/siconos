/* Siconos-IO, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef SiconosFullNumerics_hpp
#define SiconosFullNumerics_hpp

#include "IOConfig.h"
#ifdef WITH_SERIALIZATION
#include "Register.hpp"


SICONOS_IO_REGISTER(NumericsOptions, (verboseMode));


template <class Archive>
void siconos_io(Archive& ar, Callback&v, unsigned int version)
{
}
REGISTER_BOOST_SERIALIZATION(Callback);

template <class Archive>
void siconos_io(Archive& ar, _SolverOptions&v, unsigned int version)
{
  SERIALIZE(v, (solverId)(isSet)(iSize)(dSize)(filterOn)(numberOfInternalSolvers), ar);

  if (Archive::is_loading::value)
  {
    v.iparam = (int *) malloc(v.iSize * sizeof(int));
    v.dparam = (double *) malloc(v.dSize * sizeof(double));
    v.internalSolvers = (SolverOptions *) malloc(v.numberOfInternalSolvers * sizeof(SolverOptions));
    v.numericsOptions = (NumericsOptions *) malloc(sizeof(NumericsOptions));
    v.callback = (Callback *) malloc(sizeof(Callback));
    v.iWork = NULL;
    v.dWork = NULL;
  }
  SERIALIZE(v, (numericsOptions)(callback), ar);

  SERIALIZE_C_ARRAY(v.iSize, v, iparam, ar);
  SERIALIZE_C_ARRAY(v.dSize, v, dparam, ar);
  SERIALIZE_C_ARRAY(v.numberOfInternalSolvers, v, internalSolvers, ar);
}
REGISTER_BOOST_SERIALIZATION(_SolverOptions);

template <class Archive>
void siconos_io(Archive& ar, LinearComplementarityProblem& v, unsigned int version)
{
  SERIALIZE(v, (size), ar);

  if(Archive::is_loading::value)
  {
    v.q = (double *) malloc(v.size * sizeof(double));
    v.M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));
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
    p.M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));
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
      v.matrix3 = NULL;
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
        v.matrix3 = NULL;
      }
      SERIALIZE(v, (matrix1), ar);
    }
  }
}
REGISTER_BOOST_SERIALIZATION(NumericsMatrix);

template <class Archive>
void siconos_io_register_Numerics(Archive& ar)
{
  ar.register_type(static_cast<_SolverOptions*>(NULL));
  ar.register_type(static_cast<LinearComplementarityProblem*>(NULL));
  ar.register_type(static_cast<NumericsMatrix*>(NULL));
  ar.register_type(static_cast<SparseBlockStructuredMatrix*>(NULL));
  ar.register_type(static_cast<FrictionContactProblem*>(NULL));

}
#endif
#endif
