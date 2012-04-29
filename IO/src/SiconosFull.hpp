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

#ifndef SiconosFull_hpp
#define SiconosFull_hpp

#include "Register.hpp"

#include <boost/typeof/typeof.hpp>

#define NVP(X) BOOST_SERIALIZATION_NVP(X)

#include <SiconosKernel.hpp>
#include <Circle.hpp>
#include <CircleCircleR.hpp>
#include <CircularDS.hpp>
#include <Disk.hpp>
#include <DiskDiskR.hpp>
#include <DiskMovingPlanR.hpp>
#include <DiskPlanR.hpp>
#include <SphereLDS.hpp>
#include <SphereLDSPlanR.hpp>
#include <SphereLDSSphereLDSR.hpp>
#include <SphereNEDS.hpp>
#include <SphereNEDSPlanR.hpp>
#include <SphereNEDSSphereNEDSR.hpp>
#include <SiconosBodies.hpp>
#include <SpaceFilter.hpp>
#include <ExternalBody.hpp>

#include "SiconosFullGenerated.hpp"

/* hand written */

SICONOS_IO_REGISTER(NumericsOptions, (verboseMode));

SICONOS_IO_REGISTER(RelationProperties, (block)(source)(target));

SICONOS_IO_REGISTER(SystemProperties, (upper_block)(lower_block));

SICONOS_IO_REGISTER(SiconosVector,)

BOOST_TYPEOF_REGISTER_TYPE(_SolverOptions);

BOOST_TYPEOF_REGISTER_TYPE(LinearComplementarityProblem);

\

template <class Archive>
void siconos_io(Archive& ar, InteractionsSet& v, unsigned int)
{
  if (Archive::is_loading::value)
  {
    v.fpt = &Interaction::getSort;
  }
  ar & boost::serialization::make_nvp("setOfT", v.setOfT);
}



template <class Archive>
void siconos_io(Archive& ar, DynamicalSystemsGraph& v, unsigned int version)
{

  ar & boost::serialization::make_nvp("g", v.g);

  if (Archive::is_loading::value)
  {
    DynamicalSystemsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = v.vertices(); vi != viend; ++vi)
    {
      v.vertex_descriptor[v.bundle(*vi)] = *vi;
    }
  }

}

template <class Archive>
void siconos_io(Archive& ar, UnitaryRelationsGraph& v, unsigned int version)
{

  ar & boost::serialization::make_nvp("g", v.g);

  if (Archive::is_loading::value)
  {
    DynamicalSystemsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = v.vertices(); vi != viend; ++vi)
    {
      v.vertex_descriptor[v.bundle(*vi)] = *vi;
    }
  }

}





template <class Archive>
void siconos_io(Archive& ar, std::basic_ofstream<char>&v , unsigned int version)
{
  // do nothing
}

template <class Archive>
void siconos_io(Archive& ar, FrictionContact &v, unsigned int version)
{
  SERIALIZE(v, (_contactProblemDim)(_mu)(_numerics_solver_options)(_numerics_solver_id), ar);

  if (Archive::is_loading::value)
  {
    if (v._contactProblemDim == 2)
      v._frictionContact_driver = &frictionContact2D_driver;
    else
      v._frictionContact_driver = &frictionContact3D_driver;
  }

  ar & boost::serialization::make_nvp("LinearOSNS",
                                      boost::serialization::base_object<LinearOSNS>(v));

}



template <class Archive>
void siconos_io(Archive& ar, __mpz_struct& v, unsigned int version)
{
  SERIALIZE(v, (_mp_alloc)(_mp_size), ar);
  SERIALIZE_C_ARRAY(v._mp_alloc, v, _mp_d, ar);
}


template <class Archive>
void siconos_io(Archive& ar, _SolverOptions&v, unsigned int version)
{
  SERIALIZE(v, (solverId)(isSet)(iSize)(dSize)(filterOn)(numberOfInternalSolvers), ar);

  SERIALIZE_C_ARRAY(v.iSize, v, iparam, ar);
  SERIALIZE_C_ARRAY(v.dSize, v, dparam, ar);
  SERIALIZE_C_ARRAY(v.numberOfInternalSolvers, v, internalSolvers, ar);
}

template <class Archive>
void siconos_io(Archive& ar, LinearComplementarityProblem& v, unsigned int version)
{
  SERIALIZE(v, (size)(M), ar);
  SERIALIZE_C_ARRAY(v.size, v, q, ar);
}

template <class Archive>
void siconos_io(Archive& ar, SparseBlockStructuredMatrix& v, unsigned int version)
{
  SERIALIZE(v, (nbblocks)(blocknumber0)(blocknumber1)(filled1)(filled2), ar);
  SERIALIZE_C_ARRAY(v.filled1, v, index1_data, ar);
  SERIALIZE_C_ARRAY(v.filled2, v, index2_data, ar);
}

template <class Archive>
void siconos_io(Archive&ar, NumericsMatrix& v, unsigned int version)
{
  SERIALIZE(v, (storageType)(size0)(size1)(matrix1), ar);
  SERIALIZE_C_ARRAY(v.size0 * v.size1, v, matrix0, ar);
}

template <class Archive>
void siconos_io(Archive& ar, DynamicalSystemsSet& v, unsigned int version)
{
  ar &  boost::serialization::make_nvp("ThisShouldNotBeASetAnyMore",
                                       boost::serialization::base_object< std::vector<SP::DynamicalSystem> >(v));
}

template <class Archive>
void siconos_io(Archive & ar, SimpleVector & v, unsigned int version)
{
  ar & boost::serialization::make_nvp("_dense", v._dense);
  if (v._dense)
  {
    ar & boost::serialization::make_nvp("vect", v.vect.Dense);
  }
  if (!v._dense)
  {
    ar & boost::serialization::make_nvp("vect", v.vect.Sparse);
  }
  ar &  boost::serialization::make_nvp("SiconosVector",
                                       boost::serialization::base_object<SiconosVector>(v));
}

template <class Archive>
void siconos_io(Archive & ar, SimpleMatrix & m, unsigned int version)
{
  ar & boost::serialization::make_nvp("num", m.num);
  ar & boost::serialization::make_nvp("dimRow", m.dimRow);
  ar & boost::serialization::make_nvp("dimCol", m.dimCol);
  ar & boost::serialization::make_nvp("ipiv", m.ipiv);
  ar & boost::serialization::make_nvp("isPLUFactorized", m.isPLUFactorized);
  ar & boost::serialization::make_nvp("isPLUInversed", m.isPLUInversed);
  switch (m.num)
  {
  case 1:
  {
    ar & boost::serialization::make_nvp("mat", m.mat.Dense);
    break;
  }
  case 2:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Triang);
    break;
  }
  case 3:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Sym);
    break;
  }
  case 4:
  {
    ar & boost::serialization::make_nvp("mat", m.mat.Sparse);
    break;
  }
  case 5:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Banded);
    break;
  }
  case 6:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Zero);
    break;
  }
  case 7:
  {
    //      *ar & boost::serialization::make_nvp("mat", c.mat.Identity);
    break;
  }
  }
  ar &  boost::serialization::make_nvp("SiconosMatrix",
                                       boost::serialization::base_object<SiconosMatrix>(m));
}


namespace boost
{
namespace serialization
{

template <class Archive>
void serialize(Archive& ar, FrictionContact& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, __mpz_struct& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, InteractionsSet& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, PluginHandle& v, unsigned int version)
{

}

template <class Archive>
void serialize(Archive& ar, UnitaryRelationsGraph& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, DynamicalSystemsGraph& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, std::basic_ofstream<char>& v, unsigned int version)
{
  siconos_io(ar, v, version);
}



template <class Archive>
void serialize(Archive& ar, NumericsMatrix& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, SparseBlockStructuredMatrix& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, LinearComplementarityProblem& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, _SolverOptions& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, DynamicalSystemsSet& v, unsigned int version)
{
  siconos_io(ar, v, version);
}


template <class Archive>
void serialize(Archive& ar, SimpleVector& v, unsigned int version)
{
  siconos_io(ar, v, version);
}

template <class Archive>
void serialize(Archive& ar, SimpleMatrix& m, unsigned int version)
{
  siconos_io(ar, m, version);
}

}
}


template <class Archive>
void siconos_io_register(Archive& ar)
{
  siconos_io_register_generated(ar);
  ar.register_type(static_cast<SimpleMatrix*>(NULL));
  ar.register_type(static_cast<SimpleVector*>(NULL));
  ar.register_type(static_cast<DynamicalSystemsSet*>(NULL));
  //  ar.register_type(static_cast<_SolverOptions*>(NULL));
  ar.register_type(static_cast<LinearComplementarityProblem*>(NULL));
  ar.register_type(static_cast<SparseBlockStructuredMatrix*>(NULL));
  ar.register_type(static_cast<NumericsMatrix*>(NULL));
  ar.register_type(static_cast<std::basic_ofstream<char>*>(NULL));
  ar.register_type(static_cast<DynamicalSystemsGraph*>(NULL));
  ar.register_type(static_cast<UnitaryRelationsGraph*>(NULL));
  //  ar.register_type(static_cast<PluginHandle*>(NULL));
  ar.register_type(static_cast<InteractionsSet*>(NULL));
  ar.register_type(static_cast<__mpz_struct*>(NULL));
  ar.register_type(static_cast<FrictionContact*>(NULL));
}

#endif
