#include <math.h>
#include <assert.h>
#include <op3x3.h>
#include <stdlib.h>

#define Sign(x) ((x>0) - (x<0))
#define Max fmax
#define Heaviside(x) (.5*Sign(x) + .5)
#define Rand(x) ((double) rand()/ (double) RAND_MAX)
#define CHECK(x) x = (isfinite(x) ? x : 0)

// temporary bug fix for overloaded pow. Sympy generates code with long double
// and it is not clear for all compiler which cast should be applied.
// The real fix is to prevent sympy from adding the type specifier
#ifdef __cplusplus
#define pow(x, y) std::pow(static_cast<double>(x), static_cast<double>(y))
#else
#define pow(x,y) pow((double)x, (double)y)
#endif

void frictionContact3D_FischerBurmeisterFABGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0=pow(ut1, 2); CHECK(x0);
  double x1=pow(ut2, 2); CHECK(x1);
  double x2=sqrt(x0 + x1); CHECK(x2);
  double x3=mu*x2; CHECK(x3);
  double x4=un + x3; CHECK(x4);
  double x5=mu*rn; CHECK(x5);
  double x6=pow(mu, 2); CHECK(x6);
  double x7=x0*x6; CHECK(x7);
  double x8=x1*x6; CHECK(x8);
  double x9=pow(rn, 2)*x6 + pow(rt1, 2) + pow(rt2, 2) + pow(x4, 2) + x7 + x8; CHECK(x9);
  double x10=mu*ut1; CHECK(x10);
  double x11=rt1*x5 + x10*x4; CHECK(x11);
  double x12=pow(x11, 2); CHECK(x12);
  double x13=mu*ut2; CHECK(x13);
  double x14=rt2*x5 + x13*x4; CHECK(x14);
  double x15=pow(x14, 2); CHECK(x15);
  double x16=x12 + x15; CHECK(x16);
  double x17=sqrt(x16); CHECK(x17);
  double x18=2*x17; CHECK(x18);
  double x19=-x18 + x9; CHECK(x19);
  double x20=sqrt(Max(0, x19)); CHECK(x20);
  double x21=0.5*x20; CHECK(x21);
  double x22=x18 + x9; CHECK(x22);
  double x23=sqrt(Max(0, x22)); CHECK(x23);
  double x24=0.5*x23; CHECK(x24);
  double x25=2*un + 2*x3; CHECK(x25);
  double x26=x10*x11; CHECK(x26);
  double x27=x13*x14; CHECK(x27);
  double x28=1.0/x17; CHECK(x28);
  double x29=2*x28; CHECK(x29);
  double x30=x29*(x26 + x27); CHECK(x30);
  double x31=0.25*Heaviside(x19)/x20; CHECK(x31);
  double x32=x31*(x25 - x30); CHECK(x32);
  double x33=0.25*Heaviside(x22)/x23; CHECK(x33);
  double x34=x33*(x25 + x30); CHECK(x34);
  double x35=1.0/x2; CHECK(x35);
  double x36=x10*x35; CHECK(x36);
  double x37=2*x6; CHECK(x37);
  double x38=ut1*x37 + x25*x36; CHECK(x38);
  double x39=ut1*ut2*x35*x6; CHECK(x39);
  double x40=x14*x39; CHECK(x40);
  double x41=mu*x4; CHECK(x41);
  double x42=2*x41; CHECK(x42);
  double x43=2*x35*x6; CHECK(x43);
  double x44=(1.0L/2.0L)*x11*(x0*x43 + x42); CHECK(x44);
  double x45=x29*(x40 + x44); CHECK(x45);
  double x46=x31*(x38 - x45); CHECK(x46);
  double x47=x33*(x38 + x45); CHECK(x47);
  double x48=x13*x35; CHECK(x48);
  double x49=ut2*x37 + x25*x48; CHECK(x49);
  double x50=x11*x39; CHECK(x50);
  double x51=(1.0L/2.0L)*x14*(x1*x43 + x42); CHECK(x51);
  double x52=x29*(x50 + x51); CHECK(x52);
  double x53=x31*(x49 - x52); CHECK(x53);
  double x54=x33*(x49 + x52); CHECK(x54);
  double x55=rn*x37; CHECK(x55);
  double x56=mu*rt1; CHECK(x56);
  double x57=x11*x56; CHECK(x57);
  double x58=mu*rt2; CHECK(x58);
  double x59=x14*x58; CHECK(x59);
  double x60=x29*(x57 + x59); CHECK(x60);
  double x61=x31*(x55 - x60); CHECK(x61);
  double x62=x33*(x55 + x60); CHECK(x62);
  double x63=2*rt1; CHECK(x63);
  double x64=2*mu*rn*x28; CHECK(x64);
  double x65=x11*x64; CHECK(x65);
  double x66=x31*(x63 - x65); CHECK(x66);
  double x67=x33*(x63 + x65); CHECK(x67);
  double x68=2*rt2; CHECK(x68);
  double x69=x14*x64; CHECK(x69);
  double x70=x31*(x68 - x69); CHECK(x70);
  double x71=x33*(x68 + x69); CHECK(x71);
  int    x72=x17 > 0;
  double x73=Rand(1); CHECK(x73);
  double x74=Rand(2); CHECK(x74);
  double x75=pow(pow(fabs(x73), 2) + pow(fabs(x74), 2), -1.0L/2.0L); CHECK(x75);
  int    x76=x17 <= 0;
  double x77;
  double x78=-x26 - x27; CHECK(x78);
  double x79=pow(x16, -3.0L/2.0L); CHECK(x79);
  double x80=x11*x79; CHECK(x80);
  double x81;
  double x82=-x40 - x44; CHECK(x82);
  double x83;
  double x84=x28*x39; CHECK(x84);
  double x85=-x50 - x51; CHECK(x85);
  double x86;
  double x87=-x57 - x59; CHECK(x87);
  double x88;
  double x89=x28*x5; CHECK(x89);
  double x90=mu*rn*x79; CHECK(x90);
  double x91;
  double x92;
  double x93;
  double x94;
  double x95=x14*x79; CHECK(x95);
  double x96;
  double x97;
  double x98;
  double x99;
  double x100;
  if (x72)
  {
    x77=x11*x28; CHECK(x77);
    x81=x10*x28 + x78*x80; CHECK(x81);
    x83=x28*(x35*x7 + x41) + x80*x82; CHECK(x83);
    x86=x80*x85 + x84; CHECK(x86);
    x88=x28*x56 + x80*x87; CHECK(x88);
    x91=-x12*x90 + x89; CHECK(x91);
    x92=-x11*x14*x90; CHECK(x92);
    x94=x14*x28; CHECK(x94);
    x96=x13*x28 + x78*x95; CHECK(x96);
    x97=x82*x95 + x84; CHECK(x97);
    x98=x28*(x35*x8 + x41) + x85*x95; CHECK(x98);
    x99=x28*x58 + x87*x95; CHECK(x99);
    x100=-x15*x90 + x89; CHECK(x100);
  };
  if (x76)
  {
    x77=x73*x75; CHECK(x77);
    x81=0; CHECK(x81);
    x83=0; CHECK(x83);
    x86=0; CHECK(x86);
    x88=0; CHECK(x88);
    x91=0; CHECK(x91);
    x92=0; CHECK(x92);
    x94=x74*x75; CHECK(x94);
    x96=0; CHECK(x96);
    x97=0; CHECK(x97);
    x98=0; CHECK(x98);
    x99=0; CHECK(x99);
    x100=0; CHECK(x100);
  };
  if (x72 || x76)
  {
    x93=x21*x92 - x24*x92; CHECK(x93);
    result[1]=rt1 + x10 + x21*x77 - x24*x77; CHECK(result[1]);
    result[4]=x21*x81 - x24*x81 + x32*x77 - x34*x77; CHECK(result[4]);
    result[7]=mu + x21*x83 - x24*x83 + x46*x77 - x47*x77; CHECK(result[7]);
    result[10]=x21*x86 - x24*x86 + x53*x77 - x54*x77; CHECK(result[10]);
    result[13]=x21*x88 - x24*x88 + x61*x77 - x62*x77; CHECK(result[13]);
    result[16]=x21*x91 - x24*x91 + x66*x77 - x67*x77 + 1; CHECK(result[16]);
    result[19]=x70*x77 - x71*x77 + x93; CHECK(result[19]);
    result[2]=rt2 + x13 + x21*x94 - x24*x94; CHECK(result[2]);
    result[5]=x21*x96 - x24*x96 + x32*x94 - x34*x94; CHECK(result[5]);
    result[8]=x21*x97 - x24*x97 + x46*x94 - x47*x94; CHECK(result[8]);
    result[11]=mu + x21*x98 - x24*x98 + x53*x94 - x54*x94; CHECK(result[11]);
    result[14]=x21*x99 - x24*x99 + x61*x94 - x62*x94; CHECK(result[14]);
    result[17]=x66*x94 - x67*x94 + x93; CHECK(result[17]);
    result[20]=x100*x21 - x100*x24 + x70*x94 - x71*x94 + 1; CHECK(result[20]);
  };
  result[0]=-x21 - x24 + x4 + x5; CHECK(result[0]);
  result[3]=-x32 - x34 + 1; CHECK(result[3]);
  result[6]=x36 - x46 - x47; CHECK(result[6]);
  result[9]=x48 - x53 - x54; CHECK(result[9]);
  result[12]=mu - x61 - x62; CHECK(result[12]);
  result[15]=-x66 - x67; CHECK(result[15]);
  result[18]=-x70 - x71; CHECK(result[18]);
}
void frictionContact3D_FischerBurmeisterFGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0=pow(ut1, 2); CHECK(x0);
  double x1=pow(ut2, 2); CHECK(x1);
  double x2=mu*sqrt(x0 + x1) + un; CHECK(x2);
  double x3=mu*rn; CHECK(x3);
  double x4=pow(mu, 2); CHECK(x4);
  double x5=pow(rn, 2)*x4 + pow(rt1, 2) + pow(rt2, 2) + x0*x4 + x1*x4 + pow(x2, 2); CHECK(x5);
  double x6=mu*ut1; CHECK(x6);
  double x7=rt1*x3 + x2*x6; CHECK(x7);
  double x8=mu*ut2; CHECK(x8);
  double x9=rt2*x3 + x2*x8; CHECK(x9);
  double x10=sqrt(pow(x7, 2) + pow(x9, 2)); CHECK(x10);
  double x11=2*x10; CHECK(x11);
  double x12=0.5*sqrt(Max(0, -x11 + x5)); CHECK(x12);
  double x13=0.5*sqrt(Max(0, x11 + x5)); CHECK(x13);
  double x14=1.0/x10; CHECK(x14);
  int    x15=x10 > 0;
  double x16=Rand(1); CHECK(x16);
  double x17=Rand(2); CHECK(x17);
  double x18=pow(pow(fabs(x16), 2) + pow(fabs(x17), 2), -1.0L/2.0L); CHECK(x18);
  int    x19=x10 <= 0;
  double x20;
  double x21;
  if (x15)
  {
    x20=x14*x7; CHECK(x20);
    x21=x14*x9; CHECK(x21);
  };
  if (x19)
  {
    x20=x16*x18; CHECK(x20);
    x21=x17*x18; CHECK(x21);
  };
  if (x15 || x19)
  {
    result[1]=rt1 + x12*x20 - x13*x20 + x6; CHECK(result[1]);
    result[2]=rt2 + x12*x21 - x13*x21 + x8; CHECK(result[2]);
  };
  result[0]=-x12 - x13 + x2 + x3; CHECK(result[0]);
}
void frictionContact3D_FischerBurmeisterABGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0=pow(ut1, 2); CHECK(x0);
  double x1=pow(ut2, 2); CHECK(x1);
  double x2=sqrt(x0 + x1); CHECK(x2);
  double x3=mu*x2; CHECK(x3);
  double x4=2*un + 2*x3; CHECK(x4);
  double x5=mu*rn; CHECK(x5);
  double x6=un + x3; CHECK(x6);
  double x7=mu*x6; CHECK(x7);
  double x8=rt1*x5 + ut1*x7; CHECK(x8);
  double x9=mu*x8; CHECK(x9);
  double x10=ut1*x9; CHECK(x10);
  double x11=rt2*x5 + ut2*x7; CHECK(x11);
  double x12=mu*x11; CHECK(x12);
  double x13=ut2*x12; CHECK(x13);
  double x14=pow(x8, 2); CHECK(x14);
  double x15=pow(x11, 2); CHECK(x15);
  double x16=x14 + x15; CHECK(x16);
  double x17=sqrt(x16); CHECK(x17);
  double x18=1.0/x17; CHECK(x18);
  double x19=2*x18; CHECK(x19);
  double x20=x19*(x10 + x13); CHECK(x20);
  double x21=pow(mu, 2); CHECK(x21);
  double x22=x0*x21; CHECK(x22);
  double x23=x1*x21; CHECK(x23);
  double x24=pow(rn, 2)*x21 + pow(rt1, 2) + pow(rt2, 2) + x22 + x23 + pow(x6, 2); CHECK(x24);
  double x25=2*x17; CHECK(x25);
  double x26=x24 - x25; CHECK(x26);
  double x27=sqrt(Max(0, x26)); CHECK(x27);
  double x28=0.25*Heaviside(x26)/x27; CHECK(x28);
  double x29=x28*(-x20 + x4); CHECK(x29);
  double x30=x24 + x25; CHECK(x30);
  double x31=sqrt(Max(0, x30)); CHECK(x31);
  double x32=0.25*Heaviside(x30)/x31; CHECK(x32);
  double x33=x32*(x20 + x4); CHECK(x33);
  double x34=1.0/x2; CHECK(x34);
  double x35=mu*x34; CHECK(x35);
  double x36=ut1*x35; CHECK(x36);
  double x37=2*x21; CHECK(x37);
  double x38=ut1*x37 + x36*x4; CHECK(x38);
  double x39=ut1*ut2*x21*x34; CHECK(x39);
  double x40=x11*x39; CHECK(x40);
  double x41=2*x7; CHECK(x41);
  double x42=2*x21*x34; CHECK(x42);
  double x43=(1.0L/2.0L)*x8*(x0*x42 + x41); CHECK(x43);
  double x44=x19*(x40 + x43); CHECK(x44);
  double x45=x28*(x38 - x44); CHECK(x45);
  double x46=x32*(x38 + x44); CHECK(x46);
  double x47=ut2*x35; CHECK(x47);
  double x48=ut2*x37 + x4*x47; CHECK(x48);
  double x49=x39*x8; CHECK(x49);
  double x50=(1.0L/2.0L)*x11*(x1*x42 + x41); CHECK(x50);
  double x51=x19*(x49 + x50); CHECK(x51);
  double x52=x28*(x48 - x51); CHECK(x52);
  double x53=x32*(x48 + x51); CHECK(x53);
  double x54=rn*x37; CHECK(x54);
  double x55=rt1*x9; CHECK(x55);
  double x56=rt2*x12; CHECK(x56);
  double x57=x19*(x55 + x56); CHECK(x57);
  double x58=x28*(x54 - x57); CHECK(x58);
  double x59=x32*(x54 + x57); CHECK(x59);
  double x60=2*rt1; CHECK(x60);
  double x61=2*mu*rn*x18; CHECK(x61);
  double x62=x61*x8; CHECK(x62);
  double x63=x28*(x60 - x62); CHECK(x63);
  double x64=x32*(x60 + x62); CHECK(x64);
  double x65=2*rt2; CHECK(x65);
  double x66=x11*x61; CHECK(x66);
  double x67=x28*(x65 - x66); CHECK(x67);
  double x68=x32*(x65 + x66); CHECK(x68);
  double x69=mu*x18; CHECK(x69);
  double x70=-x10 - x13; CHECK(x70);
  double x71=pow(x16, -3.0L/2.0L); CHECK(x71);
  double x72=x71*x8; CHECK(x72);
  int    x73=x17 > 0;
  int    x74=x17 <= 0;
  double x75;
  double x76=Rand(1); CHECK(x76);
  double x77=Rand(2); CHECK(x77);
  double x78=pow(pow(fabs(x76), 2) + pow(fabs(x77), 2), -1.0L/2.0L); CHECK(x78);
  double x79;
  double x80=-x40 - x43; CHECK(x80);
  double x81;
  double x82=x18*x39; CHECK(x82);
  double x83=-x49 - x50; CHECK(x83);
  double x84;
  double x85=-x55 - x56; CHECK(x85);
  double x86;
  double x87=x18*x5; CHECK(x87);
  double x88=mu*rn*x71; CHECK(x88);
  double x89;
  double x90;
  double x91;
  double x92=x11*x71; CHECK(x92);
  double x93;
  double x94;
  double x95;
  double x96;
  double x97;
  double x98;
  if (x73)
  {
    x75=0.5*ut1*x69 + 0.5*x70*x72; CHECK(x75);
    x79=x18*x8; CHECK(x79);
    x81=0.5*x18*(x22*x34 + x7) + 0.5*x72*x80; CHECK(x81);
    x84=0.5*x72*x83 + 0.5*x82; CHECK(x84);
    x86=0.5*rt1*x69 + 0.5*x72*x85; CHECK(x86);
    x89=-0.5*x14*x88 + 0.5*x87; CHECK(x89);
    x90=-0.5*x11*x8*x88; CHECK(x90);
    x93=0.5*ut2*x69 + 0.5*x70*x92; CHECK(x93);
    x94=x11*x18; CHECK(x94);
    x95=0.5*x80*x92 + 0.5*x82; CHECK(x95);
    x96=0.5*x18*(x23*x34 + x7) + 0.5*x83*x92; CHECK(x96);
    x97=0.5*rt2*x69 + 0.5*x85*x92; CHECK(x97);
    x98=-0.5*x15*x88 + 0.5*x87; CHECK(x98);
  };
  if (x74)
  {
    x75=0; CHECK(x75);
    x79=x76*x78; CHECK(x79);
    x81=0; CHECK(x81);
    x84=0; CHECK(x84);
    x86=0; CHECK(x86);
    x89=0; CHECK(x89);
    x90=0; CHECK(x90);
    x93=0; CHECK(x93);
    x94=x77*x78; CHECK(x94);
    x95=0; CHECK(x95);
    x96=0; CHECK(x96);
    x97=0; CHECK(x97);
    x98=0; CHECK(x98);
  };
  if (x73 || x74)
  {
    x91=x27*x90 - x31*x90; CHECK(x91);
    result[1]=x27*x75 + x29*x79 - x31*x75 - x33*x79; CHECK(result[1]);
    result[4]=mu + x27*x81 - x31*x81 + x45*x79 - x46*x79; CHECK(result[4]);
    result[7]=x27*x84 - x31*x84 + x52*x79 - x53*x79; CHECK(result[7]);
    result[10]=x27*x86 - x31*x86 + x58*x79 - x59*x79; CHECK(result[10]);
    result[13]=x27*x89 - x31*x89 + x63*x79 - x64*x79 + 1; CHECK(result[13]);
    result[16]=x67*x79 - x68*x79 + x91; CHECK(result[16]);
    result[2]=x27*x93 + x29*x94 - x31*x93 - x33*x94; CHECK(result[2]);
    result[5]=x27*x95 - x31*x95 + x45*x94 - x46*x94; CHECK(result[5]);
    result[8]=mu + x27*x96 - x31*x96 + x52*x94 - x53*x94; CHECK(result[8]);
    result[11]=x27*x97 - x31*x97 + x58*x94 - x59*x94; CHECK(result[11]);
    result[14]=x63*x94 - x64*x94 + x91; CHECK(result[14]);
    result[17]=x27*x98 - x31*x98 + x67*x94 - x68*x94 + 1; CHECK(result[17]);
  };
  result[0]=-x29 - x33 + 1; CHECK(result[0]);
  result[3]=x36 - x45 - x46; CHECK(result[3]);
  result[6]=x47 - x52 - x53; CHECK(result[6]);
  result[9]=mu - x58 - x59; CHECK(result[9]);
  result[12]=-x63 - x64; CHECK(result[12]);
  result[15]=-x67 - x68; CHECK(result[15]);
}


void frictionContact3D_FischerBurmeisterFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{
  double result[21]; //3 + 2 * 9

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);


  if (f && A && B)
  {

    frictionContact3D_FischerBurmeisterFABGenerated(
      *reaction0, *reaction1, *reaction2,
      *velocity0, *velocity1, *velocity2,
      mu,
      *rho0, *rho1, *rho2,
      result);
    cpy3(result, f);
    cpy3x3(result + 3, A);
    cpy3x3(result + 12, B);
  }

  else
  {
    if (f)
    {
      frictionContact3D_FischerBurmeisterFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      frictionContact3D_FischerBurmeisterABGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3x3(result, A);
      cpy3x3(result + 9, B);
    }
  }
}
