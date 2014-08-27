#include "OccTest.hpp"

#include <MechanicsFwd.hpp>
#include "OccContactShape.hpp"
#include <TopoDS_Shape.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>

void OccTest::setUp()
{
}

void OccTest::tearDown()
{
}

void OccTest::t1()
{

  BRepPrimAPI_MakeSphere mksphere(1.0);

  OccContactShape sphere(mksphere.Shape());

  std::string s1 = sphere.exportBRepAsString();

  std::stringstream out;

  BRepTools::Write(mksphere.Shape(), out);

  CPPUNIT_ASSERT(out.str() == s1);

}

void OccTest::t2()
{
}
