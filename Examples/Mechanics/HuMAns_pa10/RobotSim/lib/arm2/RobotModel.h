// Robot parametrisation
#define NDOF      2
#define FRIC1     50
#define FRIC2     50


extern void
Ndof(int *ndof);

extern void
NLEffects(double N[NDOF], double q[NDOF], double qdot[NDOF]);

extern void
Inertia(double M[NDOF*NDOF], const double q[NDOF]);


extern void
Friction(double F[NDOF],
         const double q[NDOF],
         const double qdot[NDOF]);

