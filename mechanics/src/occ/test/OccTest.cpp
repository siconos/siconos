#include "OccTest.hpp"

#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <gp_Quaternion.hxx>
#include <gp_XYZ.hxx>
#include <iostream>

#include "ContactPoint.hpp"
#include "ContactShapeDistance.hpp"
#include "Geometer.hpp"
#include "OccBody.hpp"
#include "OccContactEdge.hpp"
#include "OccContactFace.hpp"
#include "OccContactShape.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(OccTest);

void OccTest::setUp() {}

void OccTest::tearDown() {}

void OccTest::exportBRepAsString() {
  BRepPrimAPI_MakeSphere mksphere(1.0);

  siconos::mechanics::occ::OccContactShape sphere(mksphere.Shape());

  auto s1 = sphere.exportBRepToString();

  std::stringstream out;

  BRepTools::Write(mksphere.Shape(), out);

  CPPUNIT_ASSERT(out.str() == s1);
}

void OccTest::computeUVBounds() {
  TopExp_Explorer exp;
  BRepPrimAPI_MakeSphere mksphere(1.0);

  siconos::mechanics::occ::OccContactShape sphere_shape =
      siconos::mechanics::occ::OccContactShape{mksphere.Shape()};

  siconos::mechanics::occ::OccContactFace sphere_face(sphere_shape, 0);

  sphere_face.computeUVBounds();

  CPPUNIT_ASSERT(std::abs(sphere_face.bsup1[0] - 6.28319) < 1e-4);
  CPPUNIT_ASSERT(std::abs(sphere_face.bsup1[1] - 1.5708) < 1e-4);

  CPPUNIT_ASSERT(std::abs(sphere_face.binf1[0] - 0.) < 1e-4);
  CPPUNIT_ASSERT(std::abs(sphere_face.binf1[1] + 1.5708) < 1e-4);

  std::cout << sphere_face.bsup1[0] << "," << sphere_face.bsup1[1] << std::endl;
  std::cout << sphere_face.binf1[0] << "," << sphere_face.binf1[1] << std::endl;
}

void OccTest::move() {
  BRepPrimAPI_MakeSphere mksphere(1.0);

  TopExp_Explorer exp;

  exp.Init(mksphere.Shape(), TopAbs_SHELL);

  TopoDS_Shell shell = TopoDS::Shell(exp.Current().Composed(mksphere.Shape().Orientation()));

  exp.Init(shell, TopAbs_FACE);

  TopoDS_Face face = TopoDS::Face(exp.Current().Composed(shell.Orientation()));

  siconos::mechanics::occ::OccContactShape sphere(mksphere.Shape());
  siconos::mechanics::occ::OccContactShapeV sphere_contact{
      std::make_shared<siconos::mechanics::occ::OccContactFace>(sphere, 0)};

  auto position = std::make_shared<siconos::algebra::SiconosVector>(7);
  auto velocity = std::make_shared<siconos::algebra::SiconosVector>(6);
  auto inertia = std::make_shared<siconos::algebra::SimpleMatrix>(3, 3);
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

  auto body =
      std::make_shared<siconos::mechanics::occ::OccBody>(position, velocity, 1, inertia);

  body->addContactShape(sphere_contact);

  auto data = [](const auto& obj) { return obj->data(); };

  gp_XYZ translat =
      std::visit(data, body->contactShape(0)).Location().Transformation().TranslationPart();

  std::cout << translat.X() << "," << translat.Y() << "," << translat.Z() << std::endl;

  CPPUNIT_ASSERT(translat.X() == 1.);
  CPPUNIT_ASSERT(translat.Y() == 2.);
  CPPUNIT_ASSERT(translat.Z() == 3.);

  gp_Quaternion rotat =
      std::visit(data, body->contactShape(0)).Location().Transformation().GetRotation();

  CPPUNIT_ASSERT(std::abs(rotat.X() - 0.44543540318737401) < 1e-9);
  CPPUNIT_ASSERT(std::abs(rotat.Y() - 0.53452248382484879) < 1e-9);
  CPPUNIT_ASSERT(std::abs(rotat.Z() - 0.62360956446232352) < 1e-9);
  CPPUNIT_ASSERT(std::abs(rotat.W() - 0.35634832254989918) < 1e-9);
}
#ifdef HAS_FORTRAN
void OccTest::distance() {
  const double pi = boost::math::constants::pi<double>();

  BRepPrimAPI_MakeSphere mksphere1(1, pi);
  BRepPrimAPI_MakeSphere mksphere2(1, pi);

  siconos::mechanics::occ::OccContactShape sphere1{mksphere1.Shape()};
  siconos::mechanics::occ::OccContactShape sphere2{mksphere2.Shape()};

  siconos::mechanics::occ::OccContactShapeV sphere1_contact{
      std::make_shared<siconos::mechanics::occ::OccContactFace>(sphere1, 0)};
  siconos::mechanics::occ::OccContactShapeV sphere2_contact{
      std::make_shared<siconos::mechanics::occ::OccContactFace>(sphere2, 0)};

  auto position1 = std::make_shared<siconos::algebra::SiconosVector>(7);
  auto position2 = std::make_shared<siconos::algebra::SiconosVector>(7);
  auto velocity = std::make_shared<siconos::algebra::SiconosVector>(6);
  auto inertia = std::make_shared<siconos::algebra::SimpleMatrix>(3, 3);
  position1->zero();
  (*position1)(0) = 0.;
  (*position1)(1) = 0.;
  (*position1)(2) = 0.;

  (*position1)(3) = 1;

  position2->zero();
  (*position2)(0) = 3.;
  (*position2)(1) = 0.;
  (*position2)(2) = 0.;

  (*position2)(3) = cos(pi / 2.);
  (*position2)(5) = sin(pi / 2.);

  velocity->zero();
  inertia->eye();

  auto body1 =
      std::make_shared<siconos::mechanics::occ::OccBody>(position1, velocity, 1, inertia);
  auto body2 =
      std::make_shared<siconos::mechanics::occ::OccBody>(position2, velocity, 1, inertia);
  body1->addContactShape(sphere1_contact);
  body2->addContactShape(sphere2_contact);

  auto binf0 = [](const auto& obj) { return obj->binf1[0]; };
  auto bsup0 = [](const auto& obj) { return obj->bsup1[0]; };
  auto binf1 = [](const auto& obj) { return obj->binf1[1]; };
  auto bsup1 = [](const auto& obj) { return obj->bsup1[1]; };

  std::cout << "umin1:" << std::visit(binf0, body1->contactShape(0)) << "\n";
  std::cout << "umax1:" << std::visit(bsup0, body1->contactShape(0)) << "\n";
  std::cout << "vmin1:" << std::visit(binf1, body1->contactShape(0)) << "\n";
  std::cout << "vmax1:" << std::visit(bsup1, body1->contactShape(0)) << "\n";
  std::cout << "umin2:" << std::visit(binf0, body2->contactShape(0)) << "\n";
  std::cout << "umax2:" << std::visit(bsup0, body2->contactShape(0)) << "\n";
  std::cout << "vmin2:" << std::visit(binf1, body2->contactShape(0)) << "\n";
  std::cout << "vmax2:" << std::visit(bsup1, body2->contactShape(0)) << "\n";

  auto data = [](const auto& obj) { return obj->data(); };

  gp_XYZ translat1 =
      std::visit(data, body1->contactShape(0)).Location().Transformation().TranslationPart();

  std::cout << "t1 " << translat1.X() << "," << translat1.Y() << "," << translat1.Z() << "\n";

  gp_XYZ translat2 =
      std::visit(data, body2->contactShape(0)).Location().Transformation().TranslationPart();

  std::cout << "t2 " << translat2.X() << "," << translat2.Y() << "," << translat2.Z() << "\n";

  auto dist = std::visit(
      siconos::mechanics::occ::Geometer<siconos::mechanics::occ::CadmbtbDistanceType>{},
      body1->contactShape(0), body2->contactShape(0));

  std::cout << dist->value << "\n";

  std::cout << dist->x1 << "," << dist->y1 << "," << dist->z1 << "\n";

  std::cout << dist->x2 << "," << dist->y2 << "," << dist->z2 << "\n";

  std::cout << dist->nx << "," << dist->ny << "," << dist->nz << "\n";

  CPPUNIT_ASSERT(std::abs(dist->value - 1.0) < 1e-9);
}
#endif
