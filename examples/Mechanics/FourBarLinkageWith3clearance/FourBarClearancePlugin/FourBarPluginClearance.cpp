

#include <math.h>
#include <iostream>

#include "RuntimeException.hpp"
#define PI 3.14159265
using namespace std;

// parameters according to Table 1
// geometrical characteristics
double l1 = 1.0;
double l2 = 4.0;
double l3 = 2.5;
double l0 = 3.0;
extern double r1;
double r2 = 0.06;
extern double r3;
double r4 = 0.06;
extern double r5;
double r6 = 0.06;
extern double Kp;
extern double lmd;
// inertial properties
double m1 = 1.0;
double m2 = 1.0;
double m3 = 1.0;
double J1 = m1*l1*l1/12.0;
double J2 = m2*l2*l2/12.0;
double J3 = m3*l3*l3/12.0;

double fcnExpression1(const double *q)
{
double result1;
result1 = sqrt(pow((q[3] - 0.5 * l2 * cos(q[1]) - l1 * cos(q[0])),2) + pow((q[4] - 0.5 * l2 * sin(q[1]) - l1 * sin(q[0])),2));
return result1;
}

double fcnExpression2(const double *q)
{
double result2;
result2 = sqrt(pow((l0 + q[5] + 0.5*l3 * cos(q[2]) - 0.5*l2 * cos(q[1])-q[3]),2) + pow((-q[4] + q[6]+0.5*l3 * sin(q[2]) - 0.5*l2 * sin(q[1])),2));
return result2;
}

double fcnExpression3(const double *q)
{
double result3;
result3 = sqrt((pow((q[5] - 0.5 * l3 * cos(q[2])),2) + pow((q[6] - 0.5 * l3 * sin(q[2])),2)));
return result3;
}
// force elements
double gravity = 9.81;

double Ix1 = m1*l1*l1/12.0;
double Ix2 = m2*l2*l2/12.0;
double Ix3 = m3*l3*l3/12.0;
double Jx1 = 0.5*(0.25*m1*l1*l1 + Ix1 + m2*l1*l1);
double Jx2 = 0.5*(0.25*m2*l2*l2 + Ix2);
double Jx3 = 0.5*(0.25*m3*l3*l3 + Ix3);
double P1 = 0.5*m2*l1*l2;


double Fphi(const double *q)
{
double phi, k1, k2, k3;
k1 = -2.0*l1*l3*sin(q[0]);
k2 = 2.0*l3*(l0-l1*cos(q[0]));
k3 = l0*l0 + l1*l1 - l2*l2 + l3*l3 - 2*l0*l1*cos(q[0]);
phi = 2.0*(PI + atan2((-k1- sqrt(k1*k1 + k2*k2 - k3*k3)),(k3-k2)));
//cout << "phi" << phi << endl;
return phi;
}

double Falfa(const double *q)
{
double alfa;
double v = Fphi(q);
alfa = atan2(-l1*sin(q[0]) + l3*sin(v), l0-l1*cos(q[0])+l3*cos(v));
//cout << "alfa " << alfa << endl;
return alfa;
}

double FS1(const double *q)
{
double s1;
double v = Fphi(q);
double v1 = Falfa(q);
s1 = (l1*sin(v - q[0]))/ (l2*sin(v1 - v));
return s1;
}

double FS2(const double *q)
{
double s2;
double v = Fphi(q);
double v1 = Falfa(q);
s2 = (l1*sin(v1 - q[0]))/ (l3*sin(v1 - v));
return s2;
}

double Fc1(const double *q)
{
double c1;
double v = Fphi(q);
double v1 = Falfa(q);
c1 = cos(q[0] - v1);
return c1;
}

double Fdtc1(const double *q)
{
double dtc1;
double v = Fphi(q);
double v1 = Falfa(q);
dtc1 = sin(v1 - q[0]);
return dtc1;
}

double Fdac1(const double *q)
{
double dac1;
double v = Fphi(q);
double v1 = Falfa(q);
dac1 = -sin(v1 - q[0]);
return dac1;
}

double Fdtg(const double *q)
{
double dtg;
double v = Fphi(q);
double v1 = Falfa(q);
dtg = -(0.5*m1*l1 + m2*l1)*gravity*cos(q[0]);
return dtg;
}

double Fdag(const double *q)
{
double dag;
double v = Fphi(q);
double v1 = Falfa(q);
dag = -0.5*m2*l2*gravity*cos(v1);
return dag;
}

double Fdpg(const double *q)
{
double dpg;
double v = Fphi(q);
double v1 = Falfa(q);
dpg = -0.5*m3*l3*gravity*cos(v);
return dpg;
}

double Ft1(const double *q)
{
double t1;
double v = Fphi(q);
double v1 = Falfa(q);
t1 = (-l1*cos(v-q[0]))/(l2*sin(v1-v));
//cout << "t1 " << t1 << endl;
return t1;
}

double Fta1(const double *q)
{
double ta1;
double v = Fphi(q);
double v1 = Falfa(q);
ta1 = (l1*cos(v-q[0])*cos(v1-v))/(l2*sin(v1-v)*sin(v1-v));
return ta1;
}

double Ftp1(const double *q)
{
double tp1;
double v = Fphi(q);
double v1 = Falfa(q);
tp1 = (2.0*l1*cos(v-q[0]))/(l2*cos(2.0*v1-2.0*v));
return tp1;
}

double Ft2(const double *q)
{
double t2;
double v = Fphi(q);
double v1 = Falfa(q);
t2 = (-l1*sin(v-q[0])*cos(v1-v))/(l2*sin(v1-v)*sin(v1-v));
//cout << "t2 " << t2 << endl;
return t2;
}

double Ftt2(const double *q)
{
double tt2;
double v = Fphi(q);
double v1 = Falfa(q);
tt2 = (l1*cos(v-q[0])*cos(v1-v))/(l2*sin(v1-v)*sin(v1-v));
return tt2;
}

double Fta2(const double *q)
{
double ta2;
double v = Fphi(q);
double v1 = Falfa(q);
ta2 = (l1*(-sin(v1-v)*sin(v1-v) + 2.0)*sin(v-q[0]))/(l2*sin(v1-v)*sin(v1-v)*sin(v1-v));
return ta2;
}

double Ftp2(const double *q)
{
double tp2;
double v = Fphi(q);
double v1 = Falfa(q);
tp2 = (l1*(sin(v1-v)*sin(v1-v)*sin(v-q[0]) - sin(v1-v)*cos(v1-v)*cos(v-q[0]) - 2.0*sin(v-q[0])))/(l2*sin(v1-v)*sin(v1-v)*sin(v1-v));
return tp2;
}

double Ft3(const double *q)
{
double t3;
double v = Fphi(q);
double v1 = Falfa(q);
t3 = (-2.0*l1*sin(v1-q[0]))/(l2*cos(2.0*v1-2.0*v) - l2);
//cout << "t3 " << t3 << endl;
return t3;
}

double Ftt3(const double *q)
{
double tt3;
double v = Fphi(q);
double v1 = Falfa(q);
tt3 = (2.0*l1*cos(-v1+q[0]))/(l2*cos(-2.0*v1+2.0*v) - l2);
return tt3;
}

double Fta3(const double *q)
{
double ta3;
double v = Fphi(q);
double v1 = Falfa(q);
ta3 = (l1*(-4.0*sin(-2.0*v1+2.0*v)*sin(-v1+q[0]) - 2.0*cos(-2.0*v1+2.0*v)*cos(-v1+q[0]) + 2.0*cos(-v1+q[0])))/(l2*(-sin(-2.0*v1+2.0*v) - 2.0*cos(-2.0*v1+2.0*v) + 2.0));
return ta3;
}

double Ftp3(const double *q)
{
double tp3;
double v = Fphi(q);
double v1 = Falfa(q);
tp3 = (4.0*l1*l2*sin(-2.0*v1+2.0*v)*sin(-v1 + q[0]))/((l2*cos(-2.0*v1+2.0*v) - l2)*(l2*cos(-2.0*v1+2.0*v) - l2));
return tp3;
}

double Ft4(const double *q)
{
double t4;
double v = Fphi(q);
double v1 = Falfa(q);
t4 = (-l1*cos(v1-q[0]))/(l3*sin(v1-v));
return t4;
}

double Ftt4(const double *q)
{
double tt4;
double v = Fphi(q);
double v1 = Falfa(q);
tt4 = (-l1*sin(v1-q[0]))/(l3*sin(v1-v));
return tt4;
}

double Fta4(const double *q)
{
double ta4;
double v = Fphi(q);
double v1 = Falfa(q);
ta4 = (l1*(-sin(-v1+q[0])*sin(v1-v) + cos(-v1+q[0])*cos(v1-v)))/(l3*sin(v1-v)*sin(v1-v));
return ta4;
}

double Ftp4(const double *q)
{
double tp4;
double v = Fphi(q);
double v1 = Falfa(q);
tp4 = (-l1*cos(-v1+q[0])*cos(v1-v))/(l3*sin(v1-v)*sin(v1-v));
return tp4;
}

double Ft5(const double *q)
{
double t5;
double v = Fphi(q);
double v1 = Falfa(q);
t5 = (2.0*l1*sin(v-q[0]))/(l3*cos(2.0*v1-2.0*v) -l3);
return t5;
}

double Ftt5(const double *q)
{
double tt5;
double v = Fphi(q);
double v1 = Falfa(q);
tt5 = (-2.0*l1*cos(v-q[0]))/(l3*cos(2.0*v1-2.0*v) -l3);
return tt5;
}

double Fta5(const double *q)
{
double ta5;
double v = Fphi(q);
double v1 = Falfa(q);
ta5 = (-4.0*l1*sin(-2.0*v1+2.0*v)*sin(v-q[0]))/(l3*(-sin(-2.0*v1+2.0*v)*sin(-2.0*v1+2.0*v) - 2.0*cos(-2.0*v1+2.0*v) + 2.0));
return ta5;
}

double Ftp5(const double *q)
{
double tp5;
double v = Fphi(q);
double v1 = Falfa(q);
tp5 = (2.0*l1*(2.0*sin(-2.0*v1+2.0*v)*sin(v-q[0]) + cos(-2.0*v1+2.0*v)*cos(v-q[0]) - cos(v-q[0])))/(l3*(-sin(-2.0*v1+2.0*v)*sin(-2.0*v1+2.0*v) - 2.0*cos(-2.0*v1+2.0*v) + 2.0));
return tp5;
}

double Ft6(const double *q)
{
double t6;
double v = Fphi(q);
double v1 = Falfa(q);
t6 = (l1*sin(v1-q[0])*cos(v1-v))/(l3*sin(v1-v)*sin(v1-v));
return t6;
}

double Ftt6(const double *q)
{
double tt6;
double v = Fphi(q);
double v1 = Falfa(q);
tt6 = (-l1*cos(-v1+q[0])*cos(v1-v))/(l3*sin(v1-v)*sin(v1-v));
return tt6;
}

double Fta6(const double *q)
{
double ta6;
double v = Fphi(q);
double v1 = Falfa(q);
ta6 = (l1*(-sin(-v1+q[0])*sin(v1-v)*sin(v1-v) + 2.0*sin(-v1+q[0]) + sin(v1-v)*cos(-v1+q[0])*cos(v1-v)))/(l3*sin(v1-v)*sin(v1-v)*sin(v1-v));
return ta6;
}

double Ftp6(const double *q)
{
double tp6;
double v = Fphi(q);
double v1 = Falfa(q);
tp6 = (l1*(sin(v1-v)*sin(v1-v) -2.0)*sin(-v1+q[0]))/(l3*sin(v1-v)*sin(v1-v)*sin(v1-v));
return tp6;
}

double Fft9(const double *q)
{
double ft9;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double Ta1 = Fta1(q);
double Tp1 = Ftp1(q);
ft9 = -s11 + (Ta1*s11) + (Tp1*s21);
return ft9;
}

double Fft10(const double *q)
{
double ft10;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double Tt2 = Ftt2(q);
double Ta2 = Fta2(q);
double Tp2 = Ftp2(q);
ft10 = Tt2 + (Ta2*s11) + (Tp2*s21);
return ft10;
}

double Fft11(const double *q)
{
double ft11;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double T1 = Ft1(q);
double T2 = Ft2(q);
double T3 = Ft3(q);
ft11 = T1 + (T2*s11) + (T3*s21);
return ft11;
}

double Fft12(const double *q)
{
double ft12;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double Tt3 = Ftt3(q);
double Ta3 = Fta3(q);
double Tp3 = Ftp3(q);
ft12 = Tt3 + (Ta3*s11) + (Tp3*s21);
return ft12;
}

double Fft13(const double *q)
{
double ft13;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double T4 = Ft4(q);
double T5 = Ft5(q);
double T6 = Ft6(q);
ft13 = T4 + (T5*s11) + (T6*s21);
return ft13;
}

double Fft14(const double *q)
{
double ft14;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double Tt4 = Ftt4(q);
double Ta4 = Fta4(q);
double Tp4 = Ftp4(q);
ft14 = Tt4 + (Ta4*s11) + (Tp4*s21);
return ft14;
}

double Fft15(const double *q)
{
double ft15;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double Tt5 = Ftt5(q);
double Ta5 = Fta5(q);
double Tp5 = Ftp5(q);
ft15 = Tt5 + (Ta5*s11) + (Tp5*s21);
return ft15;
}

double Fft16(const double *q)
{
double ft16;
double v = Fphi(q);
double v1 = Falfa(q);
double s11 = FS1(q);
double s21 = FS2(q);
double Tt6 = Ftt6(q);
double Ta6 = Fta6(q);
double Tp6 = Ftp6(q);
ft16 = Tt6 + (Ta6*s11) + (Tp6*s21);
return ft16;
}


double Fft17(const double *q)
{
double ft17;
double v = Fphi(q);
double v1 = Falfa(q);
ft17 = -cos(v1 - q[0]);
return ft17;
}
double Fft18(const double *q)
{
double ft18;
double v = Fphi(q);
double v1 = Falfa(q);
ft18 = cos(v1 - q[0]);
return ft18;
}

double Fft19(const double *q)
{
double ft19;
double v = Fphi(q);
double v1 = Falfa(q);
ft19 = cos(v1 - q[0]);
return ft19;
}
double Fft20(const double *q)
{
double ft20;
double v = Fphi(q);
double v1 = Falfa(q);
ft20 = -cos(v1 - q[0]);
return ft20;
}

// plugins for smooth equations of motion
double MASS1(const double *q)
{
  // columnwise definition of mass matrix
double mass1;
double s11 = FS1(q);
double s21 = FS2(q);
double c11 = Fc1(q);
//double st1 = Fdts1(q);
//double sa1 = Fdas1(q);
//double sp1 = Fdps1(q);
//double st2 = Fdts2(q);
//double sa2 = Fdas2(q);
//double sp2 = Fdps2(q);
//double ct1 = Fdtc1(q);
//double ca1 = Fdac1(q);
//double gt1 = Fdtg(q);
//double ga1 = Fdag(q);
//double gp1 = Fdpg(q);

  mass1 = 2*(J1 + (J2*s11*s11) + (J3*s21*s21) + (P1*c11*s11));
return mass1;
//cout << "mass" << mass[0] << endl;
}

double NonNL1(const double *q,const double *velocity)
{
  // nonlinear inertia terms (negative in h according to LagrangianDS)
double nonnl1;
double s11 = FS1(q);
double s21 = FS2(q);
double c11 = Fc1(q);
double T7 = Fdtc1(q);
double T8 = Fdac1(q);
double T11 = Fft11(q);
double T13 = Fft13(q);
double ct1 = Fdtc1(q);
double ca1 = Fdac1(q);
double gt1 = Fdtg(q);
double ga1 = Fdag(q);
double gp1 = Fdpg(q);

  nonnl1 = (2*J2*s11*T11 + 2*J3*s21*T13 + P1*(c11*T11 + s11*(T7 + s11*T8)))*velocity[0]*velocity[0];
return nonnl1;
 //NNL[0] = (2*J2*s11*(st1+s11*sa1+s21*sp1) + 2*J3*s21*(st2+s11*sa2+s21*sp2) + P1*(c11*(st1+s11*sa1+s21*sp1) + s11*(T7 + s11*T8)))*velocity[0]*velocity[0];
//cout << "NNL" << NNL[0] << endl;


}


// plugins for smooth equations of motion
extern "C" void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  // columnwise definition of mass matrix
  mass[0] = J1 + (0.25 * m1) * l1 * l1;
  mass[1] = 0.0;
  mass[2] = 0.0;
  mass[3] = 0.0;
  mass[4] = 0.0;
  mass[5] = 0.0;
  mass[6] = 0.0;

  mass[7] = 0.0;
  mass[8] = J2;
  mass[9] = 0.0;
  mass[10] = 0.0;
  mass[11] = 0.0;
  mass[12] = 0.0;
  mass[13] = 0.0;

  mass[14] = 0.0;
  mass[15] = 0.0;
  mass[16] = J3 ;
  mass[17] = 0.0;
  mass[18] = 0.0;
  mass[19] = 0.0;
  mass[20] = 0.0;

  mass[21] = 0.0 ;
  mass[22] = 0.0;
  mass[23] = 0.0;
  mass[24] = m2;
  mass[25] = 0.0;
  mass[26] = 0.0;
  mass[27] = 0.0 ;

  mass[28] = 0.0;
  mass[29] = 0.0;
  mass[30] = 0.0;
  mass[31] = 0.0;
  mass[32] = m2;
  mass[33] = 0.0;
  mass[34] = 0.0;

  mass[35] = 0.0;
  mass[36] = 0.0;
  mass[37] = 0.0;
  mass[38] = 0.0;
  mass[39] = 0.0;
  mass[40] = m3;
  mass[41] = 0.0;

  mass[42] = 0.0;
  mass[43] = 0.0;
  mass[44] = 0.0;
  mass[45] = 0.0;
  mass[46] = 0.0;
  mass[47] = 0.0;
  mass[48] = m3;

}

extern "C" void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, unsigned int sizeZ, double* z)
{
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  NNL[0] = 0.0;
  NNL[1] = 0.0;
  NNL[2] = 0.0;
  NNL[3] = 0.0;
  NNL[4] = 0.0;
  NNL[5] = 0.0;
  NNL[6] = 0.0;
}

extern "C" void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of nonlinear inertia terms with respect to q (columnwise)
  jacob[0] = 0.0;
  jacob[1] = 0.0;
  jacob[2] = 0.0;
  jacob[3] = 0.0;
  jacob[4] = 0.0;
  jacob[5] = 0.0;
  jacob[6] = 0.0;

  jacob[7] = 0.0;
  jacob[8] = 0.0;
  jacob[9] = 0.0;
  jacob[10] = 0.0;
  jacob[11] = 0.0;
  jacob[12] = 0.0;
  jacob[13] = 0.0;

  jacob[14] = 0.0;
  jacob[15] = 0.0;
  jacob[16] = 0.0;
  jacob[17] = 0.0;
  jacob[18] = 0.0;
  jacob[19] = 0.0;
  jacob[20] = 0.0;

  jacob[21] = 0.0;
  jacob[22] = 0.0;
  jacob[23] = 0.0;
  jacob[24] = 0.0;
  jacob[25] = 0.0;
  jacob[26] = 0.0;
  jacob[27] = 0.0;

  jacob[28] = 0.0;
  jacob[29] = 0.0;
  jacob[30] = 0.0;
  jacob[31] = 0.0;
  jacob[32] = 0.0;
  jacob[33] = 0.0;
  jacob[34] = 0.0;

  jacob[35] = 0.0;
  jacob[36] = 0.0;
  jacob[37] = 0.0;
  jacob[38] = 0.0;
  jacob[39] = 0.0;
  jacob[40] = 0.0;
  jacob[41] = 0.0;

  jacob[42] = 0.0;
  jacob[43] = 0.0;
  jacob[44] = 0.0;
  jacob[45] = 0.0;
  jacob[46] = 0.0;
  jacob[47] = 0.0;
  jacob[48] = 0.0;

}

extern "C" void jacobianFGyrVelocity(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of nonlinear inertia terms with respect to velocity (columnwise)
  jacob[0] = 0.0;
  jacob[1] = 0.0;
  jacob[2] = 0.0;
  jacob[3] = 0.0;
  jacob[4] = 0.0;
  jacob[5] = 0.0;
  jacob[6] = 0.0;

  jacob[7] = 0.0;
  jacob[8] = 0.0;
  jacob[9] = 0.0;
  jacob[10] = 0.0;
  jacob[11] = 0.0;
  jacob[12] = 0.0;
  jacob[13] = 0.0;

  jacob[14] = 0.0;
  jacob[15] = 0.0;
  jacob[16] = 0.0;
  jacob[17] = 0.0;
  jacob[18] = 0.0;
  jacob[19] = 0.0;
  jacob[20] = 0.0;

  jacob[21] = 0.0;
  jacob[22] = 0.0;
  jacob[23] = 0.0;
  jacob[24] = 0.0;
  jacob[25] = 0.0;
  jacob[26] = 0.0;
  jacob[27] = 0.0;

  jacob[28] = 0.0;
  jacob[29] = 0.0;
  jacob[30] = 0.0;
  jacob[31] = 0.0;
  jacob[32] = 0.0;
  jacob[33] = 0.0;
  jacob[34] = 0.0;

  jacob[35] = 0.0;
  jacob[36] = 0.0;
  jacob[37] = 0.0;
  jacob[38] = 0.0;
  jacob[39] = 0.0;
  jacob[40] = 0.0;
  jacob[41] = 0.0;

  jacob[42] = 0.0;
  jacob[43] = 0.0;
  jacob[44] = 0.0;
  jacob[45] = 0.0;
  jacob[46] = 0.0;
  jacob[47] = 0.0;
  jacob[48] = 0.0;
}

extern "C" void FInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double* z)
{
  //std::cout << "compute Fint" << std::endl;
  // internal forces (negative in h according to LagrangianDS)
double s11 = FS1(q);
double s21 = FS2(q);
double c11 = Fc1(q);
double T7 = Fdtc1(q);
double T8 = Fdac1(q);
double T11 = Fft11(q);
double T13 = Fft13(q);
double ct1 = Fdtc1(q);
double ca1 = Fdac1(q);
double gt1 = Fdtg(q);
double ga1 = Fdag(q);
double gp1 = Fdpg(q);
double mass11 = MASS1(q);
double nonnl11 = NonNL1(q,velocity);

  fInt[0] = (0.5 * m1) * gravity * l1 * cos(q[0])-(2*(Jx1 + (Jx2*s11*s11) + (Jx3*s21*s21) + (P1*c11*s11))*(-6.0*0.75*0.75*PI*PI*sin(0.75*PI*time) - lmd*(velocity[0] - 0.75*PI*6.0*cos(0.75*PI*time)))) -((2*Jx2*s11*T11 + 2*J3*s21*T13 + P1*(c11*T11 + s11*(T7 + s11*T8)))*velocity[0])*(velocity[0]-lmd*(q[0]-6.0*sin(0.75*PI*time))) + Kp*(velocity[0] - 0.75*PI*6.0*cos(0.75*PI*time)) + Kp*lmd*(q[0]-6.0*sin(0.75*PI*time)) - (-gt1 - s11*ga1 - s21*gp1);

  fInt[1] = 0.0;
  fInt[2] = 0.0;
  fInt[3] = 0.0;
  fInt[4] = m2 * gravity;
  fInt[5] = 0.0;
  fInt[6] = m3 * gravity;
  z[0]= mass11;
  z[1]= nonnl11;
  z[2]=(2*(Jx1 + (Jx2*s11*s11) + (Jx3*s21*s21) + (P1*c11*s11))*(-6.0*0.75*0.75*PI*PI*sin(0.75*PI*time) - lmd*(velocity[0] - 0.75*PI*6.0*cos(0.75*PI*time)))) +((2*Jx2*s11*T11 + 2*J3*s21*T13 + P1*(c11*T11 + s11*(T7 + s11*T8)))*velocity[0])*(velocity[0]-lmd*(q[0]-6.0*sin(0.75*PI*time))) - Kp*(velocity[0] - 0.75*PI*6.0*cos(0.75*PI*time)) - Kp*lmd*(q[0]-6.0*sin(0.75*PI*time)) + (-gt1 - s11*ga1 - s21*gp1);
}

extern "C" void jacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of internal forces with respect to q (columnwise)
  jacob[0] = -(0.5 * m1) * gravity * l1 * sin(q[0]);
  jacob[1] = 0.0;
  jacob[2] = 0.0;
  jacob[3] = 0.0;
  jacob[4] = 0.0;
  jacob[5] = 0.0;
  jacob[6] = 0.0;

  jacob[7] = 0.0;
  jacob[8] = 0.0;
  jacob[9] = 0.0;
  jacob[10] = 0.0;
  jacob[11] = 0.0;
  jacob[12] = 0.0;
  jacob[13] = 0.0;

  jacob[14] = 0.0;
  jacob[15] = 0.0;
  jacob[16] = 0.0;
  jacob[17] = 0.0;
  jacob[18] = 0.0;
  jacob[19] = 0.0;
  jacob[20] = 0.0;

  jacob[21] = 0.0;
  jacob[22] = 0.0;
  jacob[23] = 0.0;
  jacob[24] = 0.0;
  jacob[25] = 0.0;
  jacob[26] = 0.0;
  jacob[27] = 0.0;

  jacob[28] = 0.0;
  jacob[29] = 0.0;
  jacob[30] = 0.0;
  jacob[31] = 0.0;
  jacob[32] = 0.0;
  jacob[33] = 0.0;
  jacob[34] = 0.0;

  jacob[35] = 0.0;
  jacob[36] = 0.0;
  jacob[37] = 0.0;
  jacob[38] = 0.0;
  jacob[39] = 0.0;
  jacob[40] = 0.0;
  jacob[41] = 0.0;

  jacob[42] = 0.0;
  jacob[43] = 0.0;
  jacob[44] = 0.0;
  jacob[45] = 0.0;
  jacob[46] = 0.0;
  jacob[47] = 0.0;
  jacob[48] = 0.0;

}

extern "C" void jacobianFIntqDot(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of internal forces with respect to velocity (columnwise)
  jacob[0] = 0.0;
  jacob[1] = 0.0;
  jacob[2] = 0.0;
  jacob[3] = 0.0;
  jacob[4] = 0.0;
  jacob[5] = 0.0;
  jacob[6] = 0.0;

  jacob[7] = 0.0;
  jacob[8] = 0.0;
  jacob[9] = 0.0;
  jacob[10] = 0.0;
  jacob[11] = 0.0;
  jacob[12] = 0.0;
  jacob[13] = 0.0;

  jacob[14] = 0.0;
  jacob[15] = 0.0;
  jacob[16] = 0.0;
  jacob[17] = 0.0;
  jacob[18] = 0.0;
  jacob[19] = 0.0;
  jacob[20] = 0.0;

  jacob[21] = 0.0;
  jacob[22] = 0.0;
  jacob[23] = 0.0;
  jacob[24] = 0.0;
  jacob[25] = 0.0;
  jacob[26] = 0.0;
  jacob[27] = 0.0;

  jacob[28] = 0.0;
  jacob[29] = 0.0;
  jacob[30] = 0.0;
  jacob[31] = 0.0;
  jacob[32] = 0.0;
  jacob[33] = 0.0;
  jacob[34] = 0.0;

  jacob[35] = 0.0;
  jacob[36] = 0.0;
  jacob[37] = 0.0;
  jacob[38] = 0.0;
  jacob[39] = 0.0;
  jacob[40] = 0.0;
  jacob[41] = 0.0;

  jacob[42] = 0.0;
  jacob[43] = 0.0;
  jacob[44] = 0.0;
  jacob[45] = 0.0;
  jacob[46] = 0.0;
  jacob[47] = 0.0;
  jacob[48] = 0.0;
}
///// Clearance ///
extern "C" void g1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  double v1 = fcnExpression1(q);
  g[0] = (r2-r1) - v1;
  g[1] =0.0;

}

extern "C" void W1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  double v1 = fcnExpression1(q);
  W[0] = (-q[3]*l1*sin(q[0]) + 0.5*l1*l2*sin(q[0]-q[1])  + q[4]*l1*cos(q[0]))/v1;
  W[1] = ((q[3]*l1*cos(q[0]) - 0.5*l1*l2*cos(q[0]-q[1]) + q[4]*l1*sin(q[0]) - l1*l1 )/v1) - r1;

  W[2] = (-0.5*q[3]*l2*sin(q[1]) -0.5*l1*l2*sin(q[0]-q[1]) +0.5*q[4]*l2*cos(q[1]))/v1;
  W[3] = ((0.5*q[3]*l2*cos(q[1]) - 0.5*l1*l2*cos(q[0]-q[1])+ 0.5*q[4]*l2*sin(q[1]) - 0.25*l2*l2 )/v1) + r2;

  W[4] = 0.0;
  W[5] = 0.0;

  W[6] = (-q[3] +l1*cos(q[0]) + 0.5*l2*cos(q[1]))/v1;
  W[7] =(q[4] - 0.5*l2*sin(q[1]) - l1*sin(q[0]))/v1;

  W[8] = (-q[4] +l1*sin(q[0]) + 0.5*l2*sin(q[1]))/v1;
  W[9] = -(q[3] - 0.5*l2*cos(q[1]) - l1*cos(q[0]))/v1;

  W[10] = 0.0;
  W[11] = 0.0;

  W[12] = 0.0;
  W[13] = 0.0;
}

extern "C" void g2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  double v2 = fcnExpression2(q);
  g[0] = (r4-r3) - v2;
  g[1] =0.0;
}

extern "C" void W2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  double v2 = fcnExpression2(q);
  W[0] = 0.0;
  W[1] = 0.0;

  W[2] = (-0.5*l0*l2*sin(q[1])-0.5*q[5]*l2*sin(q[1])+0.5*q[3]*l2*sin(q[1])-0.25*l2*l3*sin(q[1]-q[2])-0.5*q[4]*l2*cos(q[1])+0.5*q[6]*l2*cos(q[1]))/v2;
  W[3] = ((0.5*l0*l2*cos(q[1])+0.5*q[5]*l2*cos(q[1])-0.5*q[3]*l2*cos(q[1])+0.25*l2*l3*cos(q[1]-q[2])-0.5*q[4]*l2*sin(q[1])+0.5*q[6]*l2*sin(q[1])-0.25*l2*l2)/v2) - r3;

  W[4] = (0.5*l0*l3*sin(q[2])+0.5*q[5]*l3*sin(q[2])-0.5*q[3]*l3*sin(q[2])+0.25*l2*l3*sin(q[1]-q[2])+0.5*q[4]*l3*cos(q[2])-0.5*q[6]*l3*cos(q[2]))/v2;
  W[5] = ((-0.5*l0*l3*cos(q[2])-0.5*q[5]*l3*cos(q[2])+0.5*q[3]*l3*cos(q[2])+0.25*l2*l3*cos(q[1]-q[2])+0.5*q[4]*l3*sin(q[2])-0.5*q[6]*l3*sin(q[2])-0.25*l3*l3)/v2) +r4;

  W[6] = (-q[3]+l0+q[5]+0.5*l3*cos(q[2])-0.5*l2*cos(q[1]))/v2;
  W[7] = (q[4]-q[6]-0.5*l3*sin(q[2])+0.5*l2*sin(q[1]))/v2;

  W[8] = (-q[4]+q[6]+0.5*l3*sin(q[2])-0.5*l2*sin(q[1]))/v2;
  W[9] = (q[5]+l0-q[3]+0.5*l3*cos(q[2])-0.5*l2*cos(q[1]))/v2;

  W[10] = (-q[5]-l0-0.5*l3*cos(q[2])+0.5*l2*cos(q[1])+q[3])/v2;
  W[11] = (q[6]-q[4]+0.5*l3*sin(q[2])-0.5*l2*sin(q[1]))/v2;

  W[12] = (-q[6]+q[4]-0.5*l3*sin(q[2])+0.5*l2*sin(q[1]))/v2;
  W[13] = (-q[5]-l0-0.5*l3*cos(q[2])+0.5*l2*cos(q[1])+q[3])/v2;
}

extern "C" void g3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  double v3 = fcnExpression3(q);
  g[0] = (r6-r5)-v3;
  g[1] =0.0;
}

extern "C" void W3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  double v3 = fcnExpression3(q);
  W[0] = 0.0;
  W[1] = 0.0;

  W[2] = 0.0;
  W[3] = 0.0;

  W[4] = (-0.5*q[5]*l3*sin(q[2])+0.5*q[6]*l3*cos(q[2]))/v3;
  W[5] = ((0.5*q[5]*l3*cos(q[2])+0.5*q[6]*l3*sin(q[2])-0.25*l3*l3)/v3)-r5;

  W[6] = 0.0;
  W[7] = 0.0;

  W[8] = 0.0;
  W[9] = 0.0;

  W[10] = (-q[5] + 0.5*l3 * cos(q[2]))/v3;
  W[11] = (q[6] -0.5* l3 * sin(q[2]))/v3;

  W[12] = (-q[6] + 0.5*l3 * sin(q[2]))/v3;
  W[13] = (-q[5] + 0.5*l3 * cos(q[2]))/v3;
}
