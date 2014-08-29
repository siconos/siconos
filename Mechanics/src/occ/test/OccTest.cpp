#include "OccTest.hpp"

#include <MechanicsFwd.hpp>
#include "OccContactShape.hpp"
#include <TopoDS_Shape.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_XYZ.hxx>
#include <gp_Quaternion.hxx>

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


#include <SiconosVector.hpp>
#include <SimpleMatrix.hpp>
#include <OccBody.hpp>
void OccTest::t3()
{
  BRepPrimAPI_MakeSphere mksphere(1.0);

  TopExp_Explorer exp;

  exp.Init(mksphere.Shape(), TopAbs_SHELL);

  TopoDS_Shell shell = TopoDS::Shell(exp.Current().Composed(mksphere.Shape().Orientation()));

  exp.Init(shell, TopAbs_FACE);

  TopoDS_Face face = TopoDS::Face(exp.Current().Composed(shell.Orientation()));

  OccContactShape sphere(face);

  SP::SiconosVector position(new SiconosVector(7));
  SP::SiconosVector velocity(new SiconosVector(6));
  SP::SimpleMatrix inertia(new SimpleMatrix(3,3));
  position->zero();
  (*position)(0) = 1.;
  (*position)(1) = 2.;
  (*position)(2) = 3.;

  /* unit quaternion from 4,5,6,7 */
  (*position)(3) = 0.35634832254989918;
  (*position)(4) = 0.44543540318737401;
  (*position)(5) = 0.53452248382484879;
  (*position)(6) = 0.62360956446232352;

  velocity->zero();
  inertia->eye();

  SP::OccBody body(new OccBody(position, velocity, 1, inertia));

  body->addContactShape(createSPtrOccContactShape(sphere));

  body->updateContactShapes();

  gp_XYZ translat = body->contactShape(0).data().Location().Transformation().
    TranslationPart();

  std::cout << translat.X() << "," << translat.Y() << "," << translat.Z()
            << std::endl;

  CPPUNIT_ASSERT(translat.X() == 1.);
  CPPUNIT_ASSERT(translat.Y() == 2.);
  CPPUNIT_ASSERT(translat.Z() == 3.);

  gp_Quaternion rotat = body->contactShape(0).data().Location().Transformation().
    GetRotation();

  CPPUNIT_ASSERT(abs(rotat.X() - 0.44543540318737401) < 1e-9);
  CPPUNIT_ASSERT(abs(rotat.Y() - 0.53452248382484879) < 1e-9);
  CPPUNIT_ASSERT(abs(rotat.Z() - 0.62360956446232352) < 1e-9);
  CPPUNIT_ASSERT(abs(rotat.W() - 0.35634832254989918) < 1e-9);

}


