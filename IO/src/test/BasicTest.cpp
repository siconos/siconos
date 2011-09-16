#include "BasicTest.hpp"
#include "../Register.hpp"

#include <boost/static_assert.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>


CPPUNIT_TEST_SUITE_REGISTRATION(BasicTest);

class point
{
public:
  virtual void dummy() {};
  int x;
  int y;
};

class colored_point : public point
{
public:
  virtual void dummy() {};
  std::string color;
};

struct label
{
  std::string name;
  double value;
};

class line
{
private:
  point* _p1;
  point* _p2;

  template<class Archive>
  friend void siconos_io(Archive&, line&, const unsigned int);

public:
  line() {};
  line(point* a, point*b) : _p1(a), _p2(b) {};

  const point& p1() const
  {
    return *_p1;
  };
  const point& p2() const
  {
    return *_p2;
  };

  void setp1Ptr(point* v)
  {
    _p1 = v;
  };
  void setp2Ptr(point* v)
  {
    _p2 = v;
  };
};



BOOST_CLASS_EXPORT_GUID(point, "point")
BOOST_CLASS_EXPORT_GUID(colored_point, "colored_point")
SICONOS_IO_REGISTER(point, (x)(y));
SICONOS_IO_REGISTER_WITH_BASE(colored_point, point, (color));

SICONOS_IO_REGISTER(label, (name)(value))

SICONOS_IO_REGISTER(line, (_p1)(_p2));




void BasicTest::setUp() {};
void BasicTest::tearDown() {};

void BasicTest::t0()
{
  point p0;
  p0.x = 3.;
  p0.y = 4.;

  point p1;

  std::ostringstream os;
  {
    boost::archive::binary_oarchive oa(os);
    oa << BOOST_SERIALIZATION_NVP(p0);
  }
  std::istringstream is(os.str());
  {
    boost::archive::binary_iarchive ia(is);
    ia >> BOOST_SERIALIZATION_NVP(p1);
  }

  CPPUNIT_ASSERT((p0.x == p1.x));
  CPPUNIT_ASSERT((p0.y == p1.y));

}


void BasicTest::t1()
{
  label l0 = { "l0", 1.0 };
  label l1;

  std::ofstream ofs("t1.xml");
  CPPUNIT_ASSERT(ofs.good());
  {
    boost::archive::xml_oarchive oa(ofs);
    boost::serialization::serialize(oa, l0, 0);
  }
  std::ifstream ifs("t1.xml");
  CPPUNIT_ASSERT(ifs.good());
  {
    boost::archive::xml_iarchive ia(ifs);
    boost::serialization::serialize(ia, l1, 0);
  }

  CPPUNIT_ASSERT((l0.name == l1.name));
  CPPUNIT_ASSERT((l0.value == l1.value));

}

void BasicTest::t2()
{

  colored_point p0, p1;
  p0.color = "red";
  p0.x = 10;
  p0.y = 99;

  std::ofstream ofs("t2.xml");
  CPPUNIT_ASSERT(ofs.good());
  {
    boost::archive::xml_oarchive oa(ofs);
    oa << boost::serialization::make_nvp("colored_point", p0);
  }
  std::ifstream ifs("t2.xml");
  CPPUNIT_ASSERT(ifs.good());
  {
    boost::archive::xml_iarchive ia(ifs);
    ia >> boost::serialization::make_nvp("colored_point", p1);
  }

  CPPUNIT_ASSERT((p0.color == p1.color));
  CPPUNIT_ASSERT((p0.x == p1.x));
  CPPUNIT_ASSERT((p0.y == p1.y));

}



void BasicTest::t3()
{

  point *p;

  point *pn;

  p = new colored_point();
  pn = new colored_point();

  p->x = 1.;
  p->y = 2.;
  static_cast<colored_point *>(p)->color = "red";


  std::ofstream ofs("t3.xml");
  CPPUNIT_ASSERT(ofs.good());
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<colored_point*>(NULL));
    oa << BOOST_SERIALIZATION_NVP(p);
  }

  std::ifstream ifs("t3.xml");
  CPPUNIT_ASSERT(ifs.good());
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<colored_point*>(NULL));
    ia >> BOOST_SERIALIZATION_NVP(pn);
  }

}

void BasicTest::t4()
{

  colored_point p;

  colored_point p1, p2;

  p1.x = 32.33444;
  p1.y = 45.34434;
  p2.x = 36.33443;
  p2.y = 34.34345;

  p1.color = "red";
  p2.color = "black";



  line l0(&p1, &p2);
  line l1(&p, &p);

  std::ofstream ofs("t4.xml");
  CPPUNIT_ASSERT(ofs.good());
  {
    boost::archive::xml_oarchive oa(ofs);
    boost::serialization::serialize(oa, l0, 0);
  }

  std::ifstream ifs("t4.xml");
  CPPUNIT_ASSERT(ifs.good());
  {
    boost::archive::xml_iarchive ia(ifs);
    boost::serialization::serialize(ia, l1, 0);
  }

  CPPUNIT_ASSERT((l0.p1().x == l1.p1().x));
  CPPUNIT_ASSERT((l0.p1().y == l1.p1().y));

  CPPUNIT_ASSERT((l0.p2().x == l1.p2().x));
  CPPUNIT_ASSERT((l0.p2().y == l1.p2().y));

  CPPUNIT_ASSERT((static_cast<const colored_point*>(&l0.p2())->color ==
                  static_cast<const colored_point*>(&l1.p2())->color));

}

