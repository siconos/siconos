// Robot parametrisation
#define NDOF      2
#define FRIC1     50
#define FRIC2     50


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

