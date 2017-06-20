#include "OccTest.hpp"

#include "MechanicsFwd.hpp"
#include "OccContactShape.hpp"
#include "OccContactFace.hpp"
#include "OccContactEdge.hpp"
#include "ContactShapeDistance.hpp"
#include "ContactPoint.hpp"
#include "WhichGeometer.hpp"
#include "WhichGeometer.hpp"

#include <TopoDS_Shape.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_XYZ.hxx>
#include <gp_Quaternion.hxx>

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include <iostream>

CPPUNIT_TEST_SUITE_REGISTRATION(OccTest);

void OccTest::setUp()
{
}

void OccTest::tearDown()
{
}

void OccTest::exportBRepAsString()
{

  BRepPrimAPI_MakeSphere mksphere(1.0);

  OccContactShape sphere(mksphere.Shape());

  std::string s1 = sphere.exportBRepToString();

  std::stringstream out;

  BRepTools::Write(mksphere.Shape(), out);

  CPPUNIT_ASSERT(out.str() == s1);

}

void OccTest::computeUVBounds()
{

  TopExp_Explorer exp;
  BRepPrimAPI_MakeSphere mksphere(1.0);

  OccContactShape sphere_shape = OccContactShape(mksphere.Shape());

  OccContactFace sphere_face(sphere_shape, 0);

  sphere_face.computeUVBounds();

  CPPUNIT_ASSERT(std::abs(sphere_face.bsup1[0] - 6.28319) < 1e-4);
  CPPUNIT_ASSERT(std::abs(sphere_face.bsup1[1] - 1.5708) < 1e-4);

  CPPUNIT_ASSERT(std::abs(sphere_face.binf1[0] - 0.) < 1e-4);
  CPPUNIT_ASSERT(std::abs(sphere_face.binf1[1] + 1.5708) < 1e-4);

  std::cout << sphere_face.bsup1[0] << "," << sphere_face.bsup1[1] << std::endl;
  std::cout << sphere_face.binf1[0] << "," << sphere_face.binf1[1] << std::endl;

}


#include <SiconosVector.hpp>
#include <SimpleMatrix.hpp>
#include <OccBody.hpp>
void OccTest::move()
{
  BRepPrimAPI_MakeSphere mksphere(1.0);

  TopExp_Explorer exp;

  exp.Init(mksphere.Shape(), TopAbs_SHELL);

  TopoDS_Shell shell = TopoDS::Shell(exp.Current().Composed(mksphere.Shape().Orientation()));

  exp.Init(shell, TopAbs_FACE);

  TopoDS_Face face = TopoDS::Face(exp.Current().Composed(shell.Orientation()));

  OccContactShape sphere(mksphere.Shape());
  OccContactFace sphere_contact(sphere, 0);

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

  body->addContactShape(createSPtrOccContactShape(sphere_contact));

  gp_XYZ translat = body->contactShape(0).data().Location().Transformation().
    TranslationPart();

  std::cout << translat.X() << "," << translat.Y() << "," << translat.Z()
            << std::endl;

  CPPUNIT_ASSERT(translat.X() == 1.);
  CPPUNIT_ASSERT(translat.Y() == 2.);
  CPPUNIT_ASSERT(translat.Z() == 3.);

  gp_Quaternion rotat = body->contactShape(0).data().Location().Transformation().
    GetRotation();

  CPPUNIT_ASSERT(std::abs(rotat.X() - 0.44543540318737401) < 1e-9);
  CPPUNIT_ASSERT(std::abs(rotat.Y() - 0.53452248382484879) < 1e-9);
  CPPUNIT_ASSERT(std::abs(rotat.Z() - 0.62360956446232352) < 1e-9);
  CPPUNIT_ASSERT(std::abs(rotat.W() - 0.35634832254989918) < 1e-9);

}

void OccTest::distance()
{
  const double pi = boost::math::constants::pi<double>();

  BRepPrimAPI_MakeSphere mksphere1(1, pi);
  BRepPrimAPI_MakeSphere mksphere2(1, pi);

  OccContactShape sphere1(mksphere1.Shape());
  OccContactShape sphere2(mksphere2.Shape());

  OccContactFace sphere1_contact(sphere1, 0);
  OccContactFace sphere2_contact(sphere2, 0);

  SP::SiconosVector position1(new SiconosVector(7));
  SP::SiconosVector position2(new SiconosVector(7));
  SP::SiconosVector velocity(new SiconosVector(6));
  SP::SimpleMatrix inertia(new SimpleMatrix(3,3));
  position1->zero();
  (*position1)(0) = 0.;
  (*position1)(1) = 0.;
  (*position1)(2) = 0.;

  (*position1)(3) = 1;

  position2->zero();
  (*position2)(0) = 3.;
  (*position2)(1) = 0.;
  (*position2)(2) = 0.;

  (*position2)(3) = cos(pi/2.);
  (*position2)(5) = sin(pi/2.);

  velocity->zero();
  inertia->eye();

  SP::OccBody body1(new OccBody(position1, velocity, 1, inertia));
  SP::OccBody body2(new OccBody(position2, velocity, 1, inertia));

  body1->addContactShape(createSPtrOccContactShape(sphere1_contact));
  body2->addContactShape(createSPtrOccContactShape(sphere2_contact));


  std::cout << "umin1:" << body1->contactShape(0).binf1[0] << std::endl;
  std::cout << "umax1:" << body1->contactShape(0).bsup1[0] << std::endl;
  std::cout << "vmin1:" << body1->contactShape(0).binf1[1] << std::endl;
  std::cout << "vmax1:" << body1->contactShape(0).bsup1[1] << std::endl;

  std::cout << "umin2:" << body2->contactShape(0).binf1[0] << std::endl;
  std::cout << "umax2:" << body2->contactShape(0).bsup1[0] << std::endl;
  std::cout << "vmin2:" << body2->contactShape(0).binf1[1] << std::endl;
  std::cout << "vmax2:" << body2->contactShape(0).bsup1[1] << std::endl;

  gp_XYZ translat1 = body1->contactShape(0).data().Location().Transformation().
    TranslationPart();

  std::cout << translat1.X() << "," << translat1.Y() << "," << translat1.Z()
            << std::endl;


  gp_XYZ translat2 = body2->contactShape(0).data().Location().Transformation().
    TranslationPart();

  std::cout << translat2.X() << "," << translat2.Y() << "," << translat2.Z()
            << std::endl;

  SP::Geometer geometer = ask<WhichGeometer<CadmbtbDistanceType> >(body1->contactShape(0));

  body2->contactShape(0).accept(*geometer);

  ContactShapeDistance& dist = geometer->answer;

  std::cout << dist.value << std::endl;

  std::cout << dist.x1 << "," << dist.y1 << "," << dist.z1 << std::endl;

  std::cout << dist.x2 << "," << dist.y2 << "," << dist.z2 << std::endl;

  std::cout << dist.nx << "," << dist.ny << "," << dist.nz << std::endl;

  CPPUNIT_ASSERT(std::abs(dist.value - 1.0) < 1e-9);

}
