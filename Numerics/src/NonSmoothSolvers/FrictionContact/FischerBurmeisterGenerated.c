#include <math.h>
#include <assert.h>
#include <op3x3.h>

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
  double x0=ut1*ut1;
  double x1=ut2*ut2;
  double x2=sqrt(x0 + x1);
  double x3=mu*x2 + un;
  double x4=mu*rn;
  double x5=mu*mu;
  double x6=x0*x5;
  double x7=x1*x5;
  double x8=pow(rn, 2)*x5 + pow(rt1, 2) + pow(rt2, 2) + pow(x3, 2) + x6 + x7;
  double x9=mu*ut1;
  double x10=rt1*x4 + x3*x9;
  double x11=x10*x10;
  double x12=mu*ut2;
  double x13=rt2*x4 + x12*x3;
  double x14=x13*x13;
  double x15=x11 + x14;
  double x16=sqrt(x15);
  double x17=2*x16;
  double x18=sqrt(-x17 + x8);
  double x19=0.5*x18;
  double x20=sqrt(x17 + x8);
  double x21=0.5*x20;
  double x22=0.5/x20;
  double x23=x10*x9;
  double x24=x12*x13;
  double x25=1./(x16);
  double x26=x25*(x23 + x24);
  double x27=x22*(x26 + x3);
  double x28=0.5/x18;
  double x29=x28*(-x26 + x3);
  double x30=1./(x2);
  double x31=x30*x9;
  double x32=ut1*x5 + x3*x31;
  double x33=ut1*ut2*x30*x5;
  double x34=x13*x33;
  double x35=mu*x3;
  double x36=2*x35;
  double x37=x30*x6;
  double x38=0.5*x10*(x36 + 2*x37);
  double x39=x25*(x34 + x38);
  double x40=x22*(x32 + x39);
  double x41=x28*(x32 - x39);
  double x42=x12*x30;
  double x43=ut2*x5 + x3*x42;
  double x44=x10*x33;
  double x45=x30*x7;
  double x46=0.5*x13*(x36 + 2*x45);
  double x47=x25*(x44 + x46);
  double x48=x22*(x43 + x47);
  double x49=x28*(x43 - x47);
  double x50=rn*x5;
  double x51=mu*rt1;
  double x52=x10*x51;
  double x53=mu*rt2;
  double x54=x13*x53;
  double x55=x25*(x52 + x54);
  double x56=x22*(x50 + x55);
  double x57=x28*(x50 - x55);
  double x58=x25*x4;
  double x59=x10*x58;
  double x60=x22*(rt1 + x59);
  double x61=x28*(rt1 - x59);
  double x62=x13*x58;
  double x63=x22*(rt2 + x62);
  double x64=x28*(rt2 - x62);
  int    x65=x16 > 0;
  int    x66=x16 <= 0;
  double x67;
  double x68=-x23 - x24;
  double x69=sqrt(1./(x15*x15*x15));
  double x70=x10*x69;
  double x71;
  double x72=-x34 - x38;
  double x73;
  double x74=x25*x33;
  double x75=-x44 - x46;
  double x76;
  double x77=-x52 - x54;
  double x78;
  double x79=mu*rn*x69;
  double x80;
  double x81;
  double x82;
  double x83;
  double x84=x13*x69;
  double x85;
  double x86;
  double x87;
  double x88;
  double x89;
  if (x65)
  {
    x67=x10*x25;
    x71=x25*x9 + x68*x70;
    x73=x25*(x35 + x37) + x70*x72;
    x76=x70*x75 + x74;
    x78=x25*x51 + x70*x77;
    x80=-x11*x79 + x58;
    x81=-x10*x13*x79;
    x83=x13*x25;
    x85=x12*x25 + x68*x84;
    x86=x72*x84 + x74;
    x87=x25*(x35 + x45) + x75*x84;
    x88=x25*x53 + x77*x84;
    x89=-x14*x79 + x58;
  };
  if (x66)
  {
    x67=1;
    x71=0;
    x73=0;
    x76=0;
    x78=0;
    x80=0;
    x81=0;
    x83=0;
    x85=0;
    x86=0;
    x87=0;
    x88=0;
    x89=0;
  };
  if (x65 || x66)
  {
    x82=x19*x81 - x21*x81;
    result[1]=rt1 + x19*x67 - x21*x67 + x9;
    result[4]=x19*x71 - x21*x71 - x27*x67 + x29*x67;
    result[7]=mu + x19*x73 - x21*x73 - x40*x67 + x41*x67;
    result[10]=x19*x76 - x21*x76 - x48*x67 + x49*x67;
    result[13]=x19*x78 - x21*x78 - x56*x67 + x57*x67;
    result[16]=x19*x80 - x21*x80 - x60*x67 + x61*x67 + 1;
    result[19]=-x63*x67 + x64*x67 + x82;
    result[2]=rt2 + x12 + x19*x83 - x21*x83;
    result[5]=x19*x85 - x21*x85 - x27*x83 + x29*x83;
    result[8]=x19*x86 - x21*x86 - x40*x83 + x41*x83;
    result[11]=mu + x19*x87 - x21*x87 - x48*x83 + x49*x83;
    result[14]=x19*x88 - x21*x88 - x56*x83 + x57*x83;
    result[17]=-x60*x83 + x61*x83 + x82;
    result[20]=x19*x89 - x21*x89 - x63*x83 + x64*x83 + 1;
  };
  result[0]=-x19 - x21 + x3 + x4;
  result[3]=-x27 - x29 + 1;
  result[6]=x31 - x40 - x41;
  result[9]=x42 - x48 - x49;
  result[12]=mu - x56 - x57;
  result[15]=-x60 - x61;
  result[18]=-x63 - x64;
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
  double x0=ut1*ut1;
  double x1=ut2*ut2;
  double x2=mu*sqrt(x0 + x1) + un;
  double x3=mu*rn;
  double x4=mu*mu;
  double x5=pow(rn, 2)*x4 + pow(rt1, 2) + pow(rt2, 2) + x0*x4 + x1*x4 + pow(x2, 2);
  double x6=mu*ut1;
  double x7=rt1*x3 + x2*x6;
  double x8=mu*ut2;
  double x9=rt2*x3 + x2*x8;
  double x10=sqrt(pow(x7, 2) + pow(x9, 2));
  double x11=2*x10;
  double x12=0.5*sqrt(-x11 + x5);
  double x13=0.5*sqrt(x11 + x5);
  double x14=1./(x10);
  int    x15=x10 > 0;
  int    x16=x10 <= 0;
  double x17;
  double x18;
  if (x15)
  {
    x17=x14*x7;
    x18=x14*x9;
  };
  if (x16)
  {
    x17=1;
    x18=0;
  };
  if (x15 || x16)
  {
    result[1]=rt1 + x12*x17 - x13*x17 + x6;
    result[2]=rt2 + x12*x18 - x13*x18 + x8;
  };
  result[0]=-x12 - x13 + x2 + x3;
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
  double x0=mu*mu;
  double x1=ut1*ut1;
  double x2=x0*x1;
  double x3=ut2*ut2;
  double x4=x0*x3;
  double x5=sqrt(x1 + x3);
  double x6=mu*x5 + un;
  double x7=pow(rn, 2)*x0 + pow(rt1, 2) + pow(rt2, 2) + x2 + x4 + pow(x6, 2);
  double x8=mu*rn;
  double x9=mu*x6;
  double x10=rt1*x8 + ut1*x9;
  double x11=x10*x10;
  double x12=rt2*x8 + ut2*x9;
  double x13=x12*x12;
  double x14=x11 + x13;
  double x15=sqrt(x14);
  double x16=2*x15;
  double x17=sqrt(x16 + x7);
  double x18=0.5/x17;
  double x19=mu*x10;
  double x20=ut1*x19;
  double x21=mu*x12;
  double x22=ut2*x21;
  double x23=1./(x15);
  double x24=x23*(x20 + x22);
  double x25=x18*(x24 + x6);
  double x26=sqrt(-x16 + x7);
  double x27=0.5/x26;
  double x28=x27*(-x24 + x6);
  double x29=1./(x5);
  double x30=mu*x29;
  double x31=ut1*x30;
  double x32=ut1*x0 + x31*x6;
  double x33=ut1*ut2*x0*x29;
  double x34=x12*x33;
  double x35=2*x9;
  double x36=x2*x29;
  double x37=0.5*x10*(x35 + 2*x36);
  double x38=x23*(x34 + x37);
  double x39=x18*(x32 + x38);
  double x40=x27*(x32 - x38);
  double x41=ut2*x30;
  double x42=ut2*x0 + x41*x6;
  double x43=x10*x33;
  double x44=x29*x4;
  double x45=0.5*x12*(x35 + 2*x44);
  double x46=x23*(x43 + x45);
  double x47=x18*(x42 + x46);
  double x48=x27*(x42 - x46);
  double x49=rn*x0;
  double x50=rt1*x19;
  double x51=rt2*x21;
  double x52=x23*(x50 + x51);
  double x53=x18*(x49 + x52);
  double x54=x27*(x49 - x52);
  double x55=x23*x8;
  double x56=x10*x55;
  double x57=x18*(rt1 + x56);
  double x58=x27*(rt1 - x56);
  double x59=x12*x55;
  double x60=x18*(rt2 + x59);
  double x61=x27*(rt2 - x59);
  double x62=0.5*x26;
  double x63=mu*x23;
  double x64=-x20 - x22;
  double x65=sqrt(1./(x14*x14*x14));
  double x66=x10*x65;
  int    x67=x15 > 0;
  int    x68=x15 <= 0;
  double x69;
  double x70=0.5*x17;
  double x71;
  double x72=-x34 - x37;
  double x73;
  double x74=x23*x33;
  double x75=-x43 - x45;
  double x76;
  double x77=-x50 - x51;
  double x78;
  double x79=mu*rn*x65;
  double x80;
  double x81;
  double x82;
  double x83=x12*x65;
  double x84;
  double x85;
  double x86;
  double x87;
  double x88;
  double x89;
  if (x67)
  {
    x69=ut1*x63 + x64*x66;
    x71=x10*x23;
    x73=x23*(x36 + x9) + x66*x72;
    x76=x66*x75 + x74;
    x78=rt1*x63 + x66*x77;
    x80=-x11*x79 + x55;
    x81=-0.5*x10*x12*x79;
    x84=ut2*x63 + x64*x83;
    x85=x12*x23;
    x86=x72*x83 + x74;
    x87=x23*(x44 + x9) + x75*x83;
    x88=rt2*x63 + x77*x83;
    x89=-x13*x79 + x55;
  };
  if (x68)
  {
    x69=0;
    x71=1;
    x73=0;
    x76=0;
    x78=0;
    x80=0;
    x81=0;
    x84=0;
    x85=0;
    x86=0;
    x87=0;
    x88=0;
    x89=0;
  };
  if (x67 || x68)
  {
    x82=-x17*x81 + x26*x81;
    result[1]=-x25*x71 + x28*x71 + x62*x69 - x69*x70;
    result[4]=mu - x39*x71 + x40*x71 + x62*x73 - x70*x73;
    result[7]=-x47*x71 + x48*x71 + x62*x76 - x70*x76;
    result[10]=-x53*x71 + x54*x71 + x62*x78 - x70*x78;
    result[13]=-x57*x71 + x58*x71 + x62*x80 - x70*x80 + 1;
    result[16]=-x60*x71 + x61*x71 + x82;
    result[2]=-x25*x85 + x28*x85 + x62*x84 - x70*x84;
    result[5]=-x39*x85 + x40*x85 + x62*x86 - x70*x86;
    result[8]=mu - x47*x85 + x48*x85 + x62*x87 - x70*x87;
    result[11]=-x53*x85 + x54*x85 + x62*x88 - x70*x88;
    result[14]=-x57*x85 + x58*x85 + x82;
    result[17]=-x60*x85 + x61*x85 + x62*x89 - x70*x89 + 1;
  };
  result[0]=-x25 - x28 + 1;
  result[3]=x31 - x39 - x40;
  result[6]=x41 - x47 - x48;
  result[9]=mu - x53 - x54;
  result[12]=-x57 - x58;
  result[15]=-x60 - x61;
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
