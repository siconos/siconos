// Robot parametrisation
#define NDDL      2
#define L1        1
#define L2        1
#define M1        10
#define M2        5
#define FRIC1     50
#define FRIC2     50

#define GRAV      9.81

extern void
modele_nddl(int *nddl);

extern void
modele_coriolis(const double q[NDDL],
                const double qdot[NDDL],
                double N[NDDL*NDDL]);

extern void
modele_gravite(const double q[NDDL],
               double G[NDDL]);

extern void
modele_inertie(const double q[NDDL],
               double M[NDDL*NDDL]);


extern void
modele_frottements(const double q[NDDL],
                   const double qdot[NDDL],
                   double F[NDDL]);

