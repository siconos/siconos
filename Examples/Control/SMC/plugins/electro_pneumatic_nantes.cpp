#ifdef _WIN32
 #define SICONOS_EXPORT extern "C" __declspec(dllexport)
 #else
 #define SICONOS_EXPORT extern "C"
 #endif

#include <stdio.h>
 #include <cmath>

#define restrict __restrict
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#define CHECK(x)
SICONOS_EXPORT void computeh(double t, unsigned n, double* restrict x, unsigned p, double* restrict lambda, double* restrict h, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double l1 = lambda[0];
double l2 = lambda[1];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  double x0=0.1*t; CHECK(x0);
  double x1=v - 0.004*cos(x0); CHECK(x1);
  double x2=sin(x0); CHECK(x2);
  h[0]=alpha*(-0.04*x2 + y) + x1; CHECK(h[0]);
  h[1]=alpha*x1 - 0.00132352941176471*pn + 0.00132352941176471*pp - 14.7058823529412*v + 0.0004*x2; CHECK(h[1]);
}
SICONOS_EXPORT void computeJachx(double t, unsigned n, double* restrict x, unsigned p, double* restrict lambda, double* restrict Jachx, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double l1 = lambda[0];
double l2 = lambda[1];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  Jachx[0]=0; CHECK(Jachx[0]);
  Jachx[2]=0; CHECK(Jachx[2]);
  Jachx[4]=1; CHECK(Jachx[4]);
  Jachx[6]=alpha; CHECK(Jachx[6]);
  Jachx[1]=0.00132352941176471; CHECK(Jachx[1]);
  Jachx[3]=-0.00132352941176471; CHECK(Jachx[3]);
  Jachx[5]=alpha - 14.7058823529412; CHECK(Jachx[5]);
  Jachx[7]=0; CHECK(Jachx[7]);
}
SICONOS_EXPORT void computeJachlambda(double t, unsigned n, double* restrict x, unsigned p, double* restrict lambda, double* restrict Jachlambda, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double l1 = lambda[0];
double l2 = lambda[1];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  Jachlambda[0]=0; CHECK(Jachlambda[0]);
  Jachlambda[2]=0; CHECK(Jachlambda[2]);
  Jachlambda[1]=0; CHECK(Jachlambda[1]);
  Jachlambda[3]=0; CHECK(Jachlambda[3]);
}
SICONOS_EXPORT void computef(double t, unsigned n, double* restrict x, double* restrict f, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  double x0=0.0054*v; CHECK(x0);
  double x1=0.0045*y; CHECK(x1);
  double x2=1.0/(x1 + 0.000374193137473672); CHECK(x2);
  double x3=1.0/(-x1 + 0.000374193137473672); CHECK(x3);
  f[0]=-pp*x0*x2 + 133.641661725*x2*(6.95536270702858e-32*pp*pp*pp*pp*pp - 1.54181893025638e-25*pp*pp*pp*pp + 1.21454472133664e-19*pp*pp*pp - 4.47823017035973e-14*pp*pp + 6.78962532904015e-9*pp + 8.53342339766871e-5); CHECK(f[0]);
  f[1]=pn*x0*x3 + 133.641661725*x3*(6.95536270702858e-32*pn*pn*pn*pn*pn - 1.54181893025638e-25*pn*pn*pn*pn + 1.21454472133664e-19*pn*pn*pn - 4.47823017035973e-14*pn*pn + 6.78962532904015e-9*pn + 8.53342339766871e-5); CHECK(f[1]);
  f[2]=-0.00132352941176471*pn + 0.00132352941176471*pp - 14.7058823529412*v; CHECK(f[2]);
  f[3]=v; CHECK(f[3]);
}
SICONOS_EXPORT void computeJacfx(double t, unsigned n, double* restrict x, double* restrict Jacfx, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  double x0=0.0045*y; CHECK(x0);
  double x1=x0 + 0.000374193137473672; CHECK(x1);
  double x2=1.0/(x1); CHECK(x2);
  double x3=0.0054*x2; CHECK(x3);
  double x4=pp*pp; CHECK(x4);
  double x5=pp*pp*pp; CHECK(x5);
  double x6=pp*pp*pp*pp; CHECK(x6);
  double x7=2.43e-5*v; CHECK(x7);
  double x8=1./(x1*x1); CHECK(x8);
  double x9=-x0 + 0.000374193137473672; CHECK(x9);
  double x10=1.0/(x9); CHECK(x10);
  double x11=0.0054*x10; CHECK(x11);
  double x12=pn*pn; CHECK(x12);
  double x13=pn*pn*pn; CHECK(x13);
  double x14=pn*pn*pn*pn; CHECK(x14);
  double x15=1./(x9*x9); CHECK(x15);
  Jacfx[0]=-v*x3 + 133.641661725*x2*(-8.95646034071946e-14*pp + 3.64363416400992e-19*x4 - 6.1672757210255e-25*x5 + 3.47768135351429e-31*x6 + 6.78962532904015e-9); CHECK(Jacfx[0]);
  Jacfx[4]=0; CHECK(Jacfx[4]);
  Jacfx[8]=-pp*x3; CHECK(Jacfx[8]);
  Jacfx[12]=pp*x7*x8 - 0.6013874777625*x8*(6.95536270702858e-32*pp*pp*pp*pp*pp + 6.78962532904015e-9*pp - 4.47823017035973e-14*x4 + 1.21454472133664e-19*x5 - 1.54181893025638e-25*x6 + 8.53342339766871e-5); CHECK(Jacfx[12]);
  Jacfx[1]=0; CHECK(Jacfx[1]);
  Jacfx[5]=v*x11 + 133.641661725*x10*(-8.95646034071946e-14*pn + 3.64363416400992e-19*x12 - 6.1672757210255e-25*x13 + 3.47768135351429e-31*x14 + 6.78962532904015e-9); CHECK(Jacfx[5]);
  Jacfx[9]=pn*x11; CHECK(Jacfx[9]);
  Jacfx[13]=pn*x15*x7 + 0.6013874777625*x15*(6.95536270702858e-32*pn*pn*pn*pn*pn + 6.78962532904015e-9*pn - 4.47823017035973e-14*x12 + 1.21454472133664e-19*x13 - 1.54181893025638e-25*x14 + 8.53342339766871e-5); CHECK(Jacfx[13]);
  Jacfx[2]=0.00132352941176471; CHECK(Jacfx[2]);
  Jacfx[6]=-0.00132352941176471; CHECK(Jacfx[6]);
  Jacfx[10]=-14.7058823529412; CHECK(Jacfx[10]);
  Jacfx[14]=0; CHECK(Jacfx[14]);
  Jacfx[3]=0; CHECK(Jacfx[3]);
  Jacfx[7]=0; CHECK(Jacfx[7]);
  Jacfx[11]=1; CHECK(Jacfx[11]);
  Jacfx[15]=0; CHECK(Jacfx[15]);
}
SICONOS_EXPORT void computeg(double t, unsigned n, double* restrict x, unsigned p, double* restrict lambda, double* restrict g, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double l1 = lambda[0];
double l2 = lambda[1];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  double x0=beta*l2 + l1; CHECK(x0);
  double x1=0.0045*y; CHECK(x1);
  double x2=133.641661725*G*x0/(x1 + 0.000374193137473672); CHECK(x2);
  double x3=pp*pp; CHECK(x3);
  double x4=pp*pp*pp; CHECK(x4);
  double x5=pp*pp*pp*pp; CHECK(x5);
  double x6=pp*pp*pp*pp*pp; CHECK(x6);
  double x7=G*x0; CHECK(x7);
  int    x8=x7 >= 0.0;
  int    x9=x7 < 0.0;
  double x10=133.641661725*G*x0/(-x1 + 0.000374193137473672); CHECK(x10);
  double x11=pn*pn; CHECK(x11);
  double x12=pn*pn*pn; CHECK(x12);
  double x13=pn*pn*pn*pn; CHECK(x13);
  double x14=pn*pn*pn*pn*pn; CHECK(x14);
  if (x8)
  {
    g[0]=x2*(-6.27567493976828e-5*pp + 4.32959552123787e-10*x3 - 1.36184599900831e-15*x4 + 1.96424682992079e-21*x5 - 1.10088063828155e-27*x6 + 14.462054466034); CHECK(g[0]);
    g[1]=-x10*(0.000105066798974752*pn - 4.59203872247267e-10*x11 + 1.14855459001783e-15*x12 - 1.38377833290208e-21*x13 + 6.42465433109007e-28*x14 - 6.59951125132568); CHECK(g[1]);
  };
  if (x9)
  {
    g[0]=x2*(0.000105066798974752*pp - 4.59203872247267e-10*x3 + 1.14855459001783e-15*x4 - 1.38377833290208e-21*x5 + 6.42465433109007e-28*x6 - 6.59951125132568); CHECK(g[0]);
    g[1]=-x10*(-6.27567493976828e-5*pn + 4.32959552123787e-10*x11 - 1.36184599900831e-15*x12 + 1.96424682992079e-21*x13 - 1.10088063828155e-27*x14 + 14.462054466034); CHECK(g[1]);
  };
  g[2]=0; CHECK(g[2]);
  g[3]=0; CHECK(g[3]);
}
SICONOS_EXPORT void computeJacgx(double t, unsigned n, double* restrict x, unsigned p, double* restrict lambda, double* restrict Jacgx, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double l1 = lambda[0];
double l2 = lambda[1];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  double x0=beta*l2 + l1; CHECK(x0);
  double x1=0.0045*y; CHECK(x1);
  double x2=x1 + 0.000374193137473672; CHECK(x2);
  double x3=133.641661725*G*x0/x2; CHECK(x3);
  double x4=pp*pp; CHECK(x4);
  double x5=pp*pp*pp; CHECK(x5);
  double x6=pp*pp*pp*pp; CHECK(x6);
  double x7=G*x0; CHECK(x7);
  int    x8=x7 >= 0.0;
  int    x9=x7 < 0.0;
  double x10=0.6013874777625*G*x0/x2*x2; CHECK(x10);
  double x11=pp*pp*pp*pp*pp; CHECK(x11);
  double x12=-x1 + 0.000374193137473672; CHECK(x12);
  double x13=133.641661725*G*x0/x12; CHECK(x13);
  double x14=pn*pn; CHECK(x14);
  double x15=pn*pn*pn; CHECK(x15);
  double x16=pn*pn*pn*pn; CHECK(x16);
  double x17=0.6013874777625*G*x0/x12*x12; CHECK(x17);
  double x18=pn*pn*pn*pn*pn; CHECK(x18);
  if (x8)
  {
    Jacgx[0]=x3*(8.65919104247574e-10*pp - 4.08553799702493e-15*x4 + 7.85698731968318e-21*x5 - 5.50440319140777e-27*x6 - 6.27567493976828e-5); CHECK(Jacgx[0]);
    Jacgx[12]=-x10*(-6.27567493976828e-5*pp - 1.10088063828155e-27*x11 + 4.32959552123787e-10*x4 - 1.36184599900831e-15*x5 + 1.96424682992079e-21*x6 + 14.462054466034); CHECK(Jacgx[12]);
    Jacgx[5]=-x13*(-9.18407744494533e-10*pn + 3.44566377005348e-15*x14 - 5.5351133316083e-21*x15 + 3.21232716554503e-27*x16 + 0.000105066798974752); CHECK(Jacgx[5]);
    Jacgx[13]=-x17*(0.000105066798974752*pn - 4.59203872247267e-10*x14 + 1.14855459001783e-15*x15 - 1.38377833290208e-21*x16 + 6.42465433109007e-28*x18 - 6.59951125132568); CHECK(Jacgx[13]);
  };
  if (x9)
  {
    Jacgx[0]=x3*(-9.18407744494533e-10*pp + 3.44566377005348e-15*x4 - 5.5351133316083e-21*x5 + 3.21232716554503e-27*x6 + 0.000105066798974752); CHECK(Jacgx[0]);
    Jacgx[12]=-x10*(0.000105066798974752*pp + 6.42465433109007e-28*x11 - 4.59203872247267e-10*x4 + 1.14855459001783e-15*x5 - 1.38377833290208e-21*x6 - 6.59951125132568); CHECK(Jacgx[12]);
    Jacgx[5]=-x13*(8.65919104247574e-10*pn - 4.08553799702493e-15*x14 + 7.85698731968318e-21*x15 - 5.50440319140777e-27*x16 - 6.27567493976828e-5); CHECK(Jacgx[5]);
    Jacgx[13]=-x17*(-6.27567493976828e-5*pn + 4.32959552123787e-10*x14 - 1.36184599900831e-15*x15 + 1.96424682992079e-21*x16 - 1.10088063828155e-27*x18 + 14.462054466034); CHECK(Jacgx[13]);
  };
  Jacgx[4]=0; CHECK(Jacgx[4]);
  Jacgx[8]=0; CHECK(Jacgx[8]);
  Jacgx[1]=0; CHECK(Jacgx[1]);
  Jacgx[9]=0; CHECK(Jacgx[9]);
  Jacgx[2]=0; CHECK(Jacgx[2]);
  Jacgx[6]=0; CHECK(Jacgx[6]);
  Jacgx[10]=0; CHECK(Jacgx[10]);
  Jacgx[14]=0; CHECK(Jacgx[14]);
  Jacgx[3]=0; CHECK(Jacgx[3]);
  Jacgx[7]=0; CHECK(Jacgx[7]);
  Jacgx[11]=0; CHECK(Jacgx[11]);
  Jacgx[15]=0; CHECK(Jacgx[15]);
}
SICONOS_EXPORT void computeJacglambda(double t, unsigned n, double* restrict x, unsigned p, double* restrict lambda, double* restrict Jacglambda, unsigned sizeZ, double* restrict z)
{
double pp = x[0];
double pn = x[1];
double v = x[2];
double y = x[3];
double l1 = lambda[0];
double l2 = lambda[1];
double G = z[0];
double beta = z[1];
double alpha = z[2];
  double x0=0.0045*y; CHECK(x0);
  double x1=133.641661725*G/(x0 + 0.000374193137473672); CHECK(x1);
  double x2=pp*pp; CHECK(x2);
  double x3=pp*pp*pp; CHECK(x3);
  double x4=pp*pp*pp*pp; CHECK(x4);
  double x5=pp*pp*pp*pp*pp; CHECK(x5);
  double x6=x1*(-6.27567493976828e-5*pp + 4.32959552123787e-10*x2 - 1.36184599900831e-15*x3 + 1.96424682992079e-21*x4 - 1.10088063828155e-27*x5 + 14.462054466034); CHECK(x6);
  double x7=G*(beta*l2 + l1); CHECK(x7);
  int    x8=x7 >= 0.0;
  double x9=x1*(0.000105066798974752*pp - 4.59203872247267e-10*x2 + 1.14855459001783e-15*x3 - 1.38377833290208e-21*x4 + 6.42465433109007e-28*x5 - 6.59951125132568); CHECK(x9);
  int    x10=x7 < 0.0;
  double x11=133.641661725*G/(-x0 + 0.000374193137473672); CHECK(x11);
  double x12=pn*pn; CHECK(x12);
  double x13=pn*pn*pn; CHECK(x13);
  double x14=pn*pn*pn*pn; CHECK(x14);
  double x15=pn*pn*pn*pn*pn; CHECK(x15);
  double x16=x11*(0.000105066798974752*pn - 4.59203872247267e-10*x12 + 1.14855459001783e-15*x13 - 1.38377833290208e-21*x14 + 6.42465433109007e-28*x15 - 6.59951125132568); CHECK(x16);
  double x17=x11*(-6.27567493976828e-5*pn + 4.32959552123787e-10*x12 - 1.36184599900831e-15*x13 + 1.96424682992079e-21*x14 - 1.10088063828155e-27*x15 + 14.462054466034); CHECK(x17);
  if (x8)
  {
    Jacglambda[0]=x6; CHECK(Jacglambda[0]);
    Jacglambda[4]=beta*x6; CHECK(Jacglambda[4]);
    Jacglambda[1]=-x16; CHECK(Jacglambda[1]);
    Jacglambda[5]=-beta*x16; CHECK(Jacglambda[5]);
  };
  if (x10)
  {
    Jacglambda[0]=x9; CHECK(Jacglambda[0]);
    Jacglambda[4]=beta*x9; CHECK(Jacglambda[4]);
    Jacglambda[1]=-x17; CHECK(Jacglambda[1]);
    Jacglambda[5]=-beta*x17; CHECK(Jacglambda[5]);
  };
  Jacglambda[2]=0; CHECK(Jacglambda[2]);
  Jacglambda[6]=0; CHECK(Jacglambda[6]);
  Jacglambda[3]=0; CHECK(Jacglambda[3]);
  Jacglambda[7]=0; CHECK(Jacglambda[7]);
}
// initial state: [501571.0  501571.0  0.0  -2.0705e-27]
