#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <math.h>

/* link masses (kg) */
#define m2 8.41
#define m3 3.51
#define m4 4.31
#define m5 3.45
#define m6 1.46
#define m7 0.24

/* center of mass (m) */
#define rz2 0.06325
#define ry3 0.08944
#define rz4 0.04609
#define ry5 0.1647
#define rz6 -0.03
#define rz7 -0.029

/* offset distance (m) */
#define d1 0.115
#define d3  0.450
#define d5  0.50
#define d7  0.08

/*  motor constants (link-referenced in kg-m^2) */
#define Im1 0.75
#define Im2 0.75
#define Im3 0.2125
#define Im4 0.2125
#define Im5 0.00575
#define Im6 0.00575
#define Im7 0.00575

/* The gravitational constant (m/s^2) */
#define g 9.81

/* Lumped Inertias of Mitsubishi PA-10 */
#define I1 0.0162
#define I2 0.125
#define I3 0.177
#define I4 0.0191
#define I5 0.0430
#define I6 -0.0241
#define I7 0.112
#define I8 1.38
#define I9 0.0063
#define I10 -0.00174
#define I11 -0.783
#define I12 -0.071
#define I13 4.16
#define I14 0.63
#define I15 -1.74
#define I16 1.01
#define I17 0.0982
#define I18 0.0895
#define I19 1.7
#define I20 0.565
#define I21 1.55
#define I22 0.372
#define I23 1.3
#define I24 0.496

/* Visquous friction */
#define FRIC1 1.0
#define FRIC2 1.0
#define FRIC3 1.0
#define FRIC4 1.0
#define FRIC5 1.0
#define FRIC6 1.0
#define FRIC7 1.0


#define N_DOF 7

#define GRAV      9.81

SICONOS_EXPORT void
modele_nddl(int *nddl);

SICONOS_EXPORT void
modele_coriolis(const double q[N_DOF],
                const double qdot[N_DOF],
                double N[N_DOF*N_DOF]);

SICONOS_EXPORT void
modele_gravite(const double q[N_DOF],
               double G[N_DOF]);

SICONOS_EXPORT void
modele_inertie(const double q[N_DOF],
               double M[N_DOF*N_DOF]);


SICONOS_EXPORT void
modele_frottements(const double q[N_DOF],
                   const double qdot[N_DOF],
                   double F[N_DOF]);

