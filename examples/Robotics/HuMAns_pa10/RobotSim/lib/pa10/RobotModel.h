// Robot parametrisation
#define NDOF      7

#define VFRIC1    10
#define VFRIC2    10
#define VFRIC3     5
#define VFRIC4     5
#define VFRIC5     2
#define VFRIC6     2
#define VFRIC7     2

#define SFRIC1     0
#define SFRIC2     0
#define SFRIC3     0
#define SFRIC4     0
#define SFRIC5     0
#define SFRIC6     0
#define SFRIC7     0


SICONOS_EXPORT void
Ndof(int *ndof);

SICONOS_EXPORT void
NLEffects(double N[NDOF], double q[NDOF], double qdot[NDOF]);

SICONOS_EXPORT void
Inertia(double M[NDOF*NDOF], const double q[NDOF]);


SICONOS_EXPORT void
Friction(double F[NDOF],
         const double q[NDOF],
         const double qdot[NDOF]);

