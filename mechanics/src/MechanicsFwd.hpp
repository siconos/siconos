
#ifndef MechanicsFwd_hpp
#define MechanicsFwd_hpp
#include <SiconosPointers.hpp>

#define MECHANICS_CLASSES()                     \
  REGISTER(SpaceFilter)                         \
  REGISTER(SiconosBodies)                       \
  REGISTER(ExternalBody)                        \
  REGISTER(Disk)                                \
  REGISTER(Circle)                              \
  REGISTER(SphereNEDS)                          \
  REGISTER(SphereLDS)                           \
  REGISTER(CircleCircleR)                       \
  REGISTER(DiskDiskR)                           \
  REGISTER(DiskPlanR)                           \
  REGISTER(DiskMovingPlanR)                     \
  REGISTER(CircularDS)                          \
  REGISTER(CircularR)                           \
  REGISTER(SphereNEDSPlanR)                     \
  REGISTER(SphereNEDSSphereNEDSR)               \
  REGISTER(SphereLDSPlanR)                      \
  REGISTER(SphereLDSSphereLDSR)                 \
  REGISTER(KneeJointR)                          \
  REGISTER(PivotJointR)                         \
  REGISTER(PrismaticJointR)                     \
  REGISTER(CylindricalJointR)                   \
  REGISTER(NewtonEulerJointR)                   \
  REGISTER(FixedJointR)                         \
  REGISTER(CouplerJointR)                       \
  REGISTER(JointStopR)                          \
  REGISTER(JointFrictionR)                      \
  REGISTER(FMatrix)                             \
  REGISTER(NSLawMatrix)                         \
  REGISTER(OccR)                                \
  REGISTER(OccBody)                             \
  REGISTER(TopoDS_Shape)                        \
  REGISTER(TopoDS_Face)                         \
  REGISTER(TopoDS_Edge)                         \
  REGISTER(TopoDS_Shapes)                       \
  REGISTER(OccContactShape)                     \
  REGISTER(OccContactFace)                      \
  REGISTER(OccContactEdge)                      \
  REGISTER(DistanceCalculatorType)              \
  REGISTER(OccDistanceType)                     \
  REGISTER(CadmbtbDistanceType)                 \
  REGISTER(ContactShapes)                       \
  REGISTER(ContactPoint)                        \
  REGISTER(ContactPoints)                       \
  REGISTER(ContactShapeDistance)                \
  REGISTER(Geometer)                            \
  REGISTER(BulletR)                             \
  REGISTER(Bullet5DR)                           \
  REGISTER(Bullet1DR)                           \
  REGISTER(Bullet2dR)                           \
  REGISTER(Bullet2d3DR)                         \
  REGISTER(RigidBodyDS)                         \
  REGISTER(RigidBody2dDS)                       \
  REGISTER(StaticBody)                          \
  REGISTER(ContactR)                            \
  REGISTER(Contact5DR)                          \
  REGISTER(Contact2dR)                          \
  REGISTER(Contact2d3DR)                        \
  REGISTER(SiconosContactor)                    \
  REGISTER(SiconosContactorSet)                 \
  REGISTER(SiconosContactorBase)                \
  REGISTER(SiconosShape)                        \
  REGISTER(SiconosSphere)                       \
  REGISTER(SiconosBox)                          \
  REGISTER(SiconosCylinder)                     \
  REGISTER(SiconosCone)                         \
  REGISTER(SiconosCapsule)                      \
  REGISTER(SiconosConvexHull)                   \
  REGISTER(SiconosPlane)                        \
  REGISTER(SiconosMesh)                         \
  REGISTER(SiconosHeightMap)                    \
  REGISTER(SiconosDisk)                         \
  REGISTER(SiconosBox2d)                        \
  REGISTER(SiconosConvexHull2d)                  \
  REGISTER(SiconosCollisionQueryResult)          \
  REGISTER(SiconosCollisionManager)             \
  REGISTER(SiconosBulletCollisionManager)       \
  REGISTER(BodyShapeRecord)

#include <SiconosVisitables.hpp>

#undef SICONOS_VISITABLES
#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()                              \
  MECHANICS_CLASSES()

#undef REGISTER
#undef REGISTER_BASE
#define REGISTER(X) DEFINE_SPTR(X);
#define REGISTER_BASE(X, Y) DEFINE_SPTR(X);
MECHANICS_CLASSES();
#undef REGISTER

#endif
