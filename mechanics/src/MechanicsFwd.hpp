
#ifndef MechanicsFwd_hpp
#define MechanicsFwd_hpp
#include <SiconosPointers.hpp>

#define MECHANICS_CLASSES()\
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
  REGISTER(FMatrix)                             \
  REGISTER(NSLawMatrix)                         \
  REGISTER(OccR)                                \
  REGISTER(OccBody)                             \
  REGISTER(OccContactShape)                     \
  REGISTER(OccContactFace)                      \
  REGISTER(OccContactEdge)                      \
  REGISTER(ContactShapes)                       \
  REGISTER(ContactPoint)                        \
  REGISTER(ContactPoints)                       \
  REGISTER(ContactShapeDistance)                \
  REGISTER(Geometer)                            \
  REGISTER(BulletDS)                            \
  REGISTER(BulletR)                             \
  REGISTER(BulletFrom1DLocalFrameR)             \
  REGISTER(BulletSpaceFilter)                   \
  REGISTER(BulletTimeStepping)                  \
  REGISTER(MBTB_FC3DContactRelation)            \
  REGISTER(MBTB_ContactRelation)                \
                                                \
  /* Proposed new Mechanics API */              \
  REGISTER(BodyDS)                              \
  REGISTER(SiconosContactor)                    \
  REGISTER(SiconosContactorSet)                 \
  REGISTER(SiconosContactorBase)                \
  REGISTER(SiconosShape)                        \
  REGISTER(SiconosSphere)                       \
  REGISTER(SiconosBox)                          \
  REGISTER(SiconosConvexHull)                   \
  REGISTER(SiconosPlane)                        \
  REGISTER(SiconosCollisionManager)             \
  REGISTER(SiconosBulletCollisionManager)

#include <SiconosVisitables.hpp>

#undef SICONOS_VISITABLES
#define SICONOS_VISITABLES() \
  KERNEL_CLASSES() \
  MECHANICS_CLASSES()

#undef REGISTER
#undef REGISTER_BASE
#define REGISTER(X) DEFINE_SPTR(X);
#define REGISTER_BASE(X, Y) DEFINE_SPTR(X);
MECHANICS_CLASSES();
#undef REGISTER

#endif
