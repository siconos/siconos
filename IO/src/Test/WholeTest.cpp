#include "WholeTest.hpp"
#include "../Register.hpp"

#include <boost/static_assert.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>


CPPUNIT_TEST_SUITE_REGISTRATION(WholeTest);

struct point
{
  int x;
  int y;
};

class line
{
private:
  point* _p1;
  point* _p2;

  template<class Archive>
  friend void save(Archive&, line&, const unsigned int);
  template<class Archive>
  friend void load(Archive&, line&, const unsigned int);

public:
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



struct colored_point : public point
{
  std::string color;
};

struct label
{
  std::string name;
  double value;
};


SICONOS_IO_REGISTER(point, (x)(y));

SICONOS_IO_REGISTER(label, (name)(value));

SICONOS_IO_REGISTER(line, (_p1)(_p2));

SICONOS_IO_REGISTER_WITH_BASE(colored_point, point, (color));



void WholeTest::setUp() {};
void WholeTest::tearDown() {};

void WholeTest::t0()
{
  point p0 = {3, 7};

  point p1;

  std::ostringstream os;
  {
    boost::archive::binary_oarchive oa(os);
    save(oa, p0, 0);
  }
  std::istringstream is(os.str());
  {
    boost::archive::binary_iarchive ia(is);
    load(ia, p1, 0);
  }

  CPPUNIT_ASSERT((p0.x == p1.x));
  CPPUNIT_ASSERT((p0.y == p1.y));

}


void WholeTest::t1()
{
  label l0 = { "l0", 1.0 };
  label l1;

  std::ofstream ofs("t1.xml");
  CPPUNIT_ASSERT(ofs.good());
  {
    boost::archive::xml_oarchive oa(ofs);
    save(oa, l0, 0);
  }
  std::ifstream ifs("t1.xml");
  CPPUNIT_ASSERT(ifs.good());
  {
    boost::archive::xml_iarchive ia(ifs);
    load(ia, l1, 0);
  }

  CPPUNIT_ASSERT((l0.name == l1.name));
  CPPUNIT_ASSERT((l0.value == l1.value));

}

void WholeTest::t2()
{

  colored_point p0, p1;
  p0.color = "red";
  p0.x = 10;
  p0.y = 99;

  std::ofstream ofs("t2.xml");
  CPPUNIT_ASSERT(ofs.good());
  {
    boost::archive::xml_oarchive oa(ofs);
    save(oa, p0, 0);
  }
  std::ifstream ifs("t2.xml");
  CPPUNIT_ASSERT(ifs.good());
  {
    boost::archive::xml_iarchive ia(ifs);
    load(ia, p1, 0);
  }

  CPPUNIT_ASSERT((p0.color == p1.color));
  CPPUNIT_ASSERT((p0.x == p1.x));
  CPPUNIT_ASSERT((p0.y == p1.y));

}



void WholeTest::t3()
{
  point p;

  point p1 = { 2, 3};
  point p2 = { 4, 5};

  line l0(&p1, &p2);

  line l1(&p, &p);

  std::ofstream ofs("t3.xml");
  CPPUNIT_ASSERT(ofs.good());
  {
    boost::archive::xml_oarchive oa(ofs);
    save(oa, l0, 0);
  }

  std::ifstream ifs("t3.xml");
  CPPUNIT_ASSERT(ifs.good());
  {
    boost::archive::xml_iarchive ia(ifs);
    load(ia, l1, 0);
  }

  CPPUNIT_ASSERT((l0.p1().x == l1.p1().x));
}

