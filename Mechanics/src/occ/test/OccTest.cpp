#include "OccTest.hpp"

#include <MechanicsFwd.hpp>
#include "OccContactShape.hpp"
#include <TopoDS_Shape.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <cmath>

CPPUNIT_TEST_SUITE_REGISTRATION(OccTest);

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

  BRepPrimAPI_MakeSphere mksphere(1.0);

  TopExp_Explorer exp;

  exp.Init(mksphere.Shape(), TopAbs_SHELL);

  TopoDS_Shell shell = TopoDS::Shell(exp.Current().Composed(mksphere.Shape().Orientation()));

  exp.Init(shell, TopAbs_FACE);

  TopoDS_Face face = TopoDS::Face(exp.Current().Composed(shell.Orientation()));

  OccContactShape sphere(face);

  sphere.computeUVBounds();

  CPPUNIT_ASSERT(std::abs(sphere.bsup1[0] - 6.28319) < 1e-4);
  CPPUNIT_ASSERT(std::abs(sphere.bsup1[1] - 1.5708) < 1e-4);

  CPPUNIT_ASSERT(std::abs(sphere.binf1[0] - 0.) < 1e-4);
  CPPUNIT_ASSERT(std::abs(sphere.binf1[1] + 1.5708) < 1e-4);

  std::cout << sphere.bsup1[0] << "," << sphere.bsup1[1] << std::endl;
  std::cout << sphere.binf1[0] << "," << sphere.binf1[1] << std::endl;

}
