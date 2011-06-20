#include "KernelTest.hpp"
#include "../Register.hpp"

#include <boost/static_assert.hpp>

#include <SiconosKernel.hpp>


#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>


#include <SiconosVisitor.hpp>


CPPUNIT_TEST_SUITE_REGISTRATION(KernelTest);

template<typename Archive>
struct SiconosSerializer : public SiconosVisitor
{
  Archive* ar;
  unsigned int version;

  SiconosSerializer(Archive& ar, unsigned int version) : ar(&ar)  {};

  void visit(SimpleVector& c)
  {
    *ar & boost::serialization::make_nvp("_dense", c._dense);
    if (c._dense)
    {
      *ar & boost::serialization::make_nvp("vect", c.vect.Dense);
    }
    if (!c._dense)
    {
      *ar & boost::serialization::make_nvp("vect", c.vect.Sparse);
    }
  }

  void visit(SimpleMatrix& c)
  {
    *ar & boost::serialization::make_nvp("num", c.num);
    switch (c.num)
    {
    case 1:
    {
      *ar & boost::serialization::make_nvp("mat", c.mat.Dense);
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
      *ar & boost::serialization::make_nvp("mat", c.mat.Sparse);
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
  }
};

void KernelTest::setUp() {};
void KernelTest::tearDown() {};

void KernelTest::t0()
{
  SP::SimpleVector q(new SimpleVector(3));
  SP::SimpleVector q0(new SimpleVector(3));

  (*q)(0) = 1.0;
  (*q)(1) = 2.0;
  (*q)(2) = 2.0;


  std::ofstream ofs("Kernelt0.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    SiconosSerializer<boost::archive::xml_oarchive> Serializer = SiconosSerializer<boost::archive::xml_oarchive>(oa, 0);
    (*q).acceptModifier(Serializer);
  }

  std::ifstream ifs("Kernelt0.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    SiconosSerializer<boost::archive::xml_iarchive> Serializer = SiconosSerializer<boost::archive::xml_iarchive>(ia, 0);
    (*q0).acceptModifier(Serializer);
  }



  CPPUNIT_ASSERT((*q0) == (*q));
}


void KernelTest::t1()
{
  SP::SimpleMatrix m1(new SimpleMatrix(3, 3));
  SP::SimpleMatrix m2(new SimpleMatrix());

  m1->eye();
  (*m1)(1, 0) = 3.0;
  (*m1)(2, 1) = -7;

  std::ofstream ofs("Kernelt1.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    SiconosSerializer<boost::archive::xml_oarchive> Serializer = SiconosSerializer<boost::archive::xml_oarchive>(oa, 0);
    (*m1).acceptModifier(Serializer);
  }

  std::ifstream ifs("Kernelt1.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    SiconosSerializer<boost::archive::xml_iarchive> Serializer = SiconosSerializer<boost::archive::xml_iarchive>(ia, 0);
    (*m2).acceptModifier(Serializer);
  }

  CPPUNIT_ASSERT((*m1) == (*m2));

}

void KernelTest::t2()
{
}

void KernelTest::t3()
{
}

