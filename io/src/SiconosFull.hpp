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

#ifndef SiconosFull_hpp
#define SiconosFull_hpp
#include "SiconosConfig.h"
#ifdef WITH_SERIALIZATION
#include "Register.hpp"

#define NVP(X) BOOST_SERIALIZATION_NVP(X)

#include <SparseBlockMatrix.h>
#ifdef KERNEL_ONLY
#include "SiconosFullKernelGenerated.hpp"
#else
//only kernel+mechanics here! <- need dispatch for kernel+control, etc.
#include "SiconosFullGenerated.hpp"
#endif
#include "SiconosFullNumerics.hpp"

#include <SimpleMatrix.hpp>
#include <SiconosVector.hpp>
#include <fc2d_Solvers.h>
#include <fc3d_Solvers.h>
/* hand written */


//BOOST_TYPEOF_REGISTER_TYPE(_SolverOptions);

//BOOST_TYPEOF_REGISTER_TYPE(LinearComplementarityProblem);


template <class Archive>
void siconos_io(Archive& ar, _DynamicalSystemsGraph& v, unsigned int version)
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
REGISTER_BOOST_SERIALIZATION(_DynamicalSystemsGraph);

template <class Archive>
void siconos_io(Archive& ar, _InteractionsGraph& v, unsigned int version)
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
REGISTER_BOOST_SERIALIZATION(_InteractionsGraph);




template <class Archive>
void siconos_io(Archive& ar, std::basic_ofstream<char>&v , unsigned int version)
{
  // do nothing
}
REGISTER_BOOST_SERIALIZATION(std::basic_ofstream<char>);

template <class Archive>
void siconos_io(Archive& ar, FrictionContact &v, unsigned int version)
{
  SERIALIZE(v, (_contactProblemDim)(_mu)(_numerics_solver_options)(_numerics_solver_id), ar);

  if (Archive::is_loading::value)
  {
    if (v._contactProblemDim == 2)
      v._frictionContact_driver = &fc2d_driver;
    else
      v._frictionContact_driver = &fc3d_driver;
  }

  ar & boost::serialization::make_nvp("LinearOSNS",
                                      boost::serialization::base_object<LinearOSNS>(v));

}
REGISTER_BOOST_SERIALIZATION(FrictionContact);


template <class Archive>
void siconos_io(Archive& ar, __mpz_struct& v, unsigned int version)
{
  SERIALIZE(v, (_mp_alloc)(_mp_size), ar);
  SERIALIZE_C_ARRAY(v._mp_alloc, v, _mp_d, ar);
}
REGISTER_BOOST_SERIALIZATION(__mpz_struct);

template <class Archive>
void siconos_io(Archive& ar, __mpf_struct& v, unsigned int version)
{
  SERIALIZE(v, (_mp_prec)(_mp_size)(_mp_exp), ar);
  SERIALIZE_C_ARRAY(abs(v._mp_size), v, _mp_d, ar);
}
REGISTER_BOOST_SERIALIZATION(__mpf_struct);






template <class Archive>
void siconos_io(Archive & ar, SiconosVector & v, unsigned int version)
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
}
REGISTER_BOOST_SERIALIZATION(SiconosVector);

template <class Archive>
void siconos_io(Archive & ar, SimpleMatrix & m, unsigned int version)
{
  ar & boost::serialization::make_nvp("num", m._num);
  ar & boost::serialization::make_nvp("_ipiv", m._ipiv);
  ar & boost::serialization::make_nvp("_isPLUFactorized", m._isPLUFactorized);
  ar & boost::serialization::make_nvp("_isPLUInversed", m._isPLUInversed);
  switch (m._num)
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
REGISTER_BOOST_SERIALIZATION(SimpleMatrix);

#include <f2c.h>
template<typename Archive>
void siconos_io(Archive& ar, LsodarOSI& osi, unsigned int version)
{
  ar & boost::serialization::make_nvp("_intData", osi._intData);

  if (Archive::is_loading::value)
  {
    osi.rtol.reset(new doublereal[osi._intData[0]]);
    osi.atol.reset(new doublereal[osi._intData[0]]);
    osi.iwork.reset(new integer[osi._intData[7]]);
    osi.rwork.reset(new doublereal[osi._intData[6]]);
    osi.jroot.reset(new integer[osi._intData[1]]);
  }
  {
    boost_ser_array<doublereal>
      wrapper = boost::serialization::make_array(osi.rtol.get(),osi._intData[0]); 
    ar & boost::serialization::make_nvp("rtol",wrapper); 
  }
  {
    boost_ser_array<doublereal>
      wrapper = boost::serialization::make_array(osi.atol.get(),osi._intData[0]); 
    ar & boost::serialization::make_nvp("atol",wrapper); 
  }
  {
    boost_ser_array<integer>
      wrapper = boost::serialization::make_array(osi.iwork.get(),osi._intData[7]); 
    ar & boost::serialization::make_nvp("iwork",wrapper); 
  }
  {
    boost_ser_array<doublereal>
      wrapper = boost::serialization::make_array(osi.rwork.get(),osi._intData[6]); 
    ar & boost::serialization::make_nvp("rwork",wrapper); 
  }
  {
    boost_ser_array<integer>
      wrapper = boost::serialization::make_array(osi.jroot.get(),osi._intData[1]); 
    ar & boost::serialization::make_nvp("jroot",wrapper); 
  }
  
  ar & boost::serialization::make_nvp("OneStepIntegrator", 
                                      boost::serialization::base_object<OneStepIntegrator>(osi));
}
REGISTER_BOOST_SERIALIZATION(LsodarOSI);


template<typename Archive, typename P>
void siconos_property_io(Archive& ar, P& p)
{

  typename P::Access::iterator vi, viend;
  for (boost::tie(vi, viend) = p.access.elements(p._g); vi != viend; ++vi)
  {
    ar & boost::serialization::make_nvp("property", (*p._store)[*vi]);
  }

}


#define MAKE_SICONOS_IO_PROPERTIES(CLASS)                               \
  template<class Archive>                                               \
  void siconos_io(Archive& ar, Siconos::VertexProperties<CLASS, _DynamicalSystemsGraph>& p, unsigned int version) \
  {                                                                     \
    SERIALIZE(p, (_g)(_stamp), ar);                                     \
    siconos_property_io(ar, p);                                         \
  }                                                                     \
  template<class Archive>                                               \
  void siconos_io(Archive& ar, Siconos::VertexProperties<CLASS, _InteractionsGraph>& p, unsigned int version) \
  {                                                                     \
    SERIALIZE(p, (_g)(_stamp), ar);                                     \
    siconos_property_io(ar, p);                                         \
  }                                                                     \
  template<class Archive>                                               \
  void siconos_io(Archive& ar, Siconos::EdgeProperties<CLASS, _DynamicalSystemsGraph>& p, unsigned int version) \
  {                                                                     \
    SERIALIZE(p, (_g)(_stamp), ar);                                     \
    siconos_property_io(ar, p);                                         \
  }                                                                     \
  template<class Archive>                                               \
  void siconos_io(Archive& ar, Siconos::EdgeProperties<CLASS, _InteractionsGraph>& p, unsigned int version) \
  {                                                                     \
    SERIALIZE(p, (_g)(_stamp), ar);                                     \
    siconos_property_io(ar, p);                                         \
  }                                                                     \

#define MAKE_SICONOS_IO_SP_PROPERTIES(CLASS)                            \
  template<class Archive>                                               \
  void siconos_io(Archive& ar, Siconos::VertexSPProperties<CLASS, _DynamicalSystemsGraph>& p, unsigned int version) \
  {                                                                     \
    SERIALIZE(p, (_g)(_stamp), ar);                                     \
    siconos_property_io(ar, p);                                         \
  }                                                                     \
  template<class Archive>                                               \
  void siconos_io(Archive& ar, Siconos::VertexSPProperties<CLASS, _InteractionsGraph>& p, unsigned int version) \
  {                                                                     \
    SERIALIZE(p, (_g)(_stamp), ar);                                     \
    siconos_property_io(ar, p);                                         \
  }                                                                     \
 
namespace Siconos
{
  MAKE_SICONOS_IO_PROPERTIES(SP::MatrixIntegrator);
  MAKE_SICONOS_IO_PROPERTIES(SP::PluggedObject);
  MAKE_SICONOS_IO_PROPERTIES(SP::OneStepIntegrator);
  MAKE_SICONOS_IO_PROPERTIES(SP::SiconosMatrix);
  MAKE_SICONOS_IO_PROPERTIES(SP::SimpleMatrix);
  MAKE_SICONOS_IO_PROPERTIES(SP::SiconosVector);
  MAKE_SICONOS_IO_SP_PROPERTIES(MatrixIntegrator);
  MAKE_SICONOS_IO_SP_PROPERTIES(PluggedObject);
  MAKE_SICONOS_IO_SP_PROPERTIES(OneStepIntegrator);
  MAKE_SICONOS_IO_SP_PROPERTIES(SiconosMatrix);
  MAKE_SICONOS_IO_SP_PROPERTIES(SimpleMatrix);
  MAKE_SICONOS_IO_SP_PROPERTIES(SiconosVector);
  MAKE_SICONOS_IO_PROPERTIES(std::string);
  MAKE_SICONOS_IO_PROPERTIES(unsigned int);
  MAKE_SICONOS_IO_PROPERTIES(double);
  MAKE_SICONOS_IO_PROPERTIES(int);
  MAKE_SICONOS_IO_PROPERTIES(bool);
}

namespace boost { namespace serialization { 
    template <class Archive, class T>
    void serialize(Archive& ar, Siconos::VertexProperties<T, _DynamicalSystemsGraph>& p, unsigned int version)
    {
      Siconos::siconos_io(ar, p, version);
    }

    template <class Archive, class T>
    void serialize(Archive& ar, Siconos::VertexSPProperties<T, _DynamicalSystemsGraph>& p, unsigned int version)
    {
      Siconos::siconos_io(ar, p, version);
    }

    template <class Archive, class T>
    void serialize(Archive& ar, Siconos::EdgeProperties<T, _DynamicalSystemsGraph>& p, unsigned int version)
    {
      siconos_io(ar, p, version);
    }

    template <class Archive, class T>
    void serialize(Archive& ar, Siconos::VertexProperties<T, _InteractionsGraph>& p, unsigned int version)
    {
      siconos_io(ar, p, version);
    }

    template <class Archive, class T>
    void serialize(Archive& ar, Siconos::VertexSPProperties<T, _InteractionsGraph>& p, unsigned int version)
    {
      siconos_io(ar, p, version);
    }

    template <class Archive, class T>
    void serialize(Archive& ar, Siconos::EdgeProperties<T, _InteractionsGraph>& p, unsigned int version)
    {
      siconos_io(ar, p, version);
    }
  } // namespace serialization
} // namespace boost




template <class Archive>
void siconos_io_register_Kernel(Archive& ar)
{
  ar.register_type(static_cast<SimpleMatrix*>(NULL));
  ar.register_type(static_cast<SiconosVector*>(NULL));

  siconos_io_register_generated(ar);

  ar.register_type(static_cast<_DynamicalSystemsGraph*>(NULL));
  ar.register_type(static_cast<_InteractionsGraph*>(NULL));
  ar.register_type(static_cast<std::basic_ofstream<char>*>(NULL));

  //  ar.register_type(static_cast<PluginHandle*>(NULL));
  ar.register_type(static_cast<__mpz_struct*>(NULL));
  ar.register_type(static_cast<FrictionContact*>(NULL));
  ar.register_type(static_cast<LsodarOSI*>(NULL));


}
#endif
#endif
