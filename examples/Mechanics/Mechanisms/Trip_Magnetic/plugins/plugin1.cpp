#include <math.h>
#include <stdio.h>
#include <cmath>
#include "MBTB_DATA.hpp"
#include <boost/math/quaternion.hpp>
#include "op3x3.h"
#define EQUIPAGE 4
#define PLATINE 3

//#define VERBOSE_DEBUG_PLUGIN

/*we assume forces are computed before the mometums, thanks siconos.*/


//static double ManetteForce[3];
extern "C" void manetteForcesFromSpring(double t, double* q, double *v,double *f, unsigned int size_z,double *z){
 /*spring K is N.mm/rad*/
  double K=0.382*10;//*10;
  double A0=1.99;
  double angle0=2*acos(z[3]);
  double angle =2*acos(q[3]);
  // point of contact of the spring.
  double P0contactSpringx=0;
  double P0contactSpringy=-5.65;
  double P0contactSpringz=-0.9;
  //double d=5.65;
  double GP=3.19;
  if (q[6]<0)
    angle=-angle;
  if (angle-angle0<-1e7)
    printf("ERROR in externalManetteForces, angle negatif=%lf\n",angle-angle0);
  double C_p=K*(A0-(angle-angle0));
#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("externalManetteForces angle = %e, angle0=%e \n",angle, angle0);
#endif
  double F_p=C_p/fabs(P0contactSpringy);
  //Mf_p e_z= F /\ 

  double GP0x=0-z[0];
  double GP0y=-5.65-z[1];
  double GP0z=-0.9-z[2];
  double Dir0x=1.0;
  double Dir0y=0;//5.65;
  double Dir0z=0;
  ::boost::math::quaternion<double>    quattrf(q[3],q[4],q[5],q[6]);
  ::boost::math::quaternion<double>    quatDir0(0,Dir0x,Dir0y,Dir0z);
  ::boost::math::quaternion<double>    cquattrf(q[3],-q[4],-q[5],-q[6]);
  ::boost::math::quaternion<double>    quatbuf;
  quatbuf=quattrf*quatDir0*cquattrf;
  f[0]=F_p*quatbuf.R_component_2();
  f[1]=F_p*quatbuf.R_component_3();
  f[2]=F_p*quatbuf.R_component_4();
  //ManetteForce[0]=f[0];
  //ManetteForce[1]=f[1];
  //ManetteForce[2]=f[2];
}

extern "C" void externalManetteForces(double t, double* q,double *v,double *f, unsigned int size_z,double *z){
  f[0]=0;
  f[1]=0;
  f[2]=0;
  manetteForcesFromSpring( t,  q,  v, f, size_z, z);
  //printf("%e %e %e\n", f[0], f[1], f[2]);
#ifdef  VERBOSE_DEBUG_PLUGIN  
  printf("externalManetteForces : The spring forces intensity is %e in the direction %e %e %e \n",F_p,quatbuf.R_component_2(),quatbuf.R_component_3(),quatbuf.R_component_4());
#endif
  f[0]=-f[0];
  f[1]=-f[1];
  f[2]=-f[2];
}


/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

extern "C" void externalManetteMomentum(double t, double* q,double *v,double *m, unsigned int size_z,double *z){
  m[0]=0;
  m[1]=0;
  m[2]=0;
  double GP_manette[3];
  double GP0x=0-z[0];
  double GP0y=-5.65-z[1];
  double GP0z=-0.9-z[2];
#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("externalManetteMomentum : GP0= %e, %e, %e\n",GP0x,GP0y,GP0z);
#endif
  ::boost::math::quaternion<double>    quattrf(q[3],q[4],q[5],q[6]);
  ::boost::math::quaternion<double>    quatGP0(0,GP0x,GP0y,GP0z);
  ::boost::math::quaternion<double>    cquattrf(q[3],-q[4],-q[5],-q[6]);
  ::boost::math::quaternion<double>    quatbuf;
  quatbuf=quattrf*quatGP0*cquattrf;
  GP_manette[0]=quatbuf.R_component_2();
  GP_manette[1]=quatbuf.R_component_3();
  GP_manette[2]=quatbuf.R_component_4();
  double f[3];
  manetteForcesFromSpring( t,  q,  v, f, size_z, z);

  cross3(GP_manette,f,m);
  //printf("externalManetteMometum: m[0] = %lf, m[1] = %lf, m[2]=%lf.\n",m[0],m[1],m[2]);
  m[0]=-m[0];
  m[1]=-m[1];
  m[2]=-m[2];

}
/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/




extern "C" void MZPlatineMomemtumDueSpringFct(double *qEquipage, double *qPlatine, double *MZPlatineMomemtumDueSpring, double *  Theta_equipage)
{
  /*PlatineEquipage is K is N.mm/rad*/
  double K2=244.46;

  // double norm2_q =0.0;
  // for (int i =0; i <7; i++){norm2_q +=(qEquipage[i]-sDS[EQUIPAGE]->q()->getValue(i))*(qEquipage[i]-sDS[EQUIPAGE]->q()->getValue(i));
  // }
  // if (sqrt(norm2_q) < 1e-12)
  // {
  //   std::cout << "qEquipage belongs to DS EQUIPAGE. accuracy "<< norm2_q <<std::endl;
  // }
  // norm2_q =0.0;
  // for (int i =0; i <7; i++){norm2_q +=(qPlatine[i]-sDS[PLATINE]->q()->getValue(i))*(qPlatine[i]-sDS[PLATINE]->q()->getValue(i));
  // }
  // if (sqrt(norm2_q) < 1e-12)
  // {
  //   std::cout << "qPlatine belongs to DS PLATINE. accuracy"<< sqrt(norm2_q) <<std::endl;
  // }
  
  *Theta_equipage=2*acos(qEquipage[3]);
  if (qEquipage[6]<0)
    *Theta_equipage=- *Theta_equipage;
/*   printf("Theta_equipage:%e\n",Theta_equipage);*/
  double Theta_platine=2*acos(qPlatine[3]);
  if (qPlatine[6]<0)
    Theta_platine=-Theta_platine;
/*   printf("Theta_platine:%e\n",Theta_platine);*/
  /*diffAngle is zero in the closed position.*/
  double diffAngle = *Theta_equipage-Theta_platine-(-54.54-22.0)*M_PI/180.0;/*-(-54.54-22.0)(-69+2).....FvsD -78.5+2 */
 // printf("diffangle1:%e\n",diffAngle);
 // printf("Equipage:%e\n",*Theta_equipage);
  //printf("Plate:%e\n",Theta_platine);
  /*The Z momentum:*/
  /*if (Theta_equipage < -69.0*M_PI/180.0)
    { 
    MZPlatineMomemtumDueSpring=-K2*(0.5236);
    }
    else
    { 
    MZPlatineMomemtumDueSpring-K2*(diffAngle + 0.5236);
    }*/
  /*The force intensity*/
  *MZPlatineMomemtumDueSpring =K2*(diffAngle + 0.49);/*0.6547, 0.9599*/


}




extern "C" void FPlatineA2(double t, double* q,double *v,double *f, unsigned int size_z,double *z){

  ::boost::math::quaternion<double>    quattrf(q[3],q[4],q[5],q[6]);
  ::boost::math::quaternion<double>    cquattrf(q[3],-q[4],-q[5],-q[6]);
  ::boost::math::quaternion<double>    quatbuf;

  f[0]=0;
  f[1]=0;
  f[2]=0;
  double g=0.0;
  /*spring A2*/
  double K=1.16;
  double L0=18.44;
  /*coordonite of the application on the batie*/
  double PbX=11.8;
  double PbY=9;
  double PbZ=1.95;
  /*P0 coordinate of the application point of the spring on platine.*/
  ::boost::math::quaternion<double>   G0P0(0,
					   0.19-3.319,
					   9-4.365,
					   -7.6+4.7);

  quatbuf=quattrf*G0P0*cquattrf;
  double Px=quatbuf.R_component_2()+q[0];
  double Py=quatbuf.R_component_3()+q[1];
  double Pz=quatbuf.R_component_4()+q[2];
  double lengthSpring = sqrt((PbX-Px)*(PbX-Px)+
			     (PbY-Py)*(PbY-Py)+
			     (PbZ-Pz)*(PbZ-Pz));

/*  printf("plugin, externalPlatineForces, spring A2 lenght:%e\n",lengthSpring);*/
#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("plugin, externalPlatineForces, spring A2 lenght:%e\n",lengthSpring);
#endif
  double nFp=K*(L0-lengthSpring);
  /*direction*/
  ::boost::math::quaternion<double>    quatDir0(0,-1,0,0);
  quatbuf=quattrf*quatDir0*cquattrf;
  f[0]=nFp*quatbuf.R_component_2();
  f[1]=nFp*quatbuf.R_component_3();
  f[2]=nFp*quatbuf.R_component_4();

  
}

extern "C" void externalPlatineForces(double t, double* q,double *v,double *f, unsigned int size_z,double *z){

  FPlatineA2( t, q, v, f, size_z, z);

#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("plugin, externalPlatineForces force due to A2 %e, %e, %e\n",f[0],f[1],f[2]);
#endif

  double MZPlatineMomemtumDueSpring ;
  double Theta_equipage;
  MZPlatineMomemtumDueSpringFct( &(*(sDS[EQUIPAGE]->q()))(0), q, &MZPlatineMomemtumDueSpring, &Theta_equipage);
  
  double F=MZPlatineMomemtumDueSpring/5.48; //updated 26 May 2015 (old value 5.0)// new value 5.48//
/*  printf("externalPlatineForces, F=%e, diffAngle contact spring =%lf\n",F,diffAngle);*/
#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("externalPlatineForces, F=%e, diffAngle contact spring =%lf\n",F,diffAngle);
#endif
  if (MZPlatineMomemtumDueSpring<0)
    printf("externalPlatineForces, M negatif:%lf\n",MZPlatineMomemtumDueSpring);

  ::boost::math::quaternion<double>    quattrf(q[3],q[4],q[5],q[6]);
  ::boost::math::quaternion<double>    cquattrf(q[3],-q[4],-q[5],-q[6]);
  ::boost::math::quaternion<double>    quatbuf;

  ::boost::math::quaternion<double>    quatDir0_2(0,-1.1,0.9,0);
  
  quatbuf=quattrf*quatDir0_2*cquattrf;
  double norm2FP=sqrt(quatbuf.R_component_2()*quatbuf.R_component_2()+quatbuf.R_component_3()*quatbuf.R_component_3()+quatbuf.R_component_4()*quatbuf.R_component_4());
  double FPlatineSpringP[3];
  FPlatineSpringP[0]=quatbuf.R_component_2()* F/norm2FP;
  FPlatineSpringP[1]=quatbuf.R_component_3()* F/norm2FP;
  FPlatineSpringP[2]=quatbuf.R_component_4()* F/norm2FP;
  
  f[0]+=FPlatineSpringP[0];
  f[1]+=FPlatineSpringP[1];
  f[2]+=FPlatineSpringP[2];
/*  printf("externalPlatineForces : The springP force is :%e,%e,%e\n",FPlatineSpringP[0],FPlatineSpringP[1],FPlatineSpringP[2]);

  printf("externalPlatineForces : The ext forces is :%e %e %e\n",f[0],f[1],f[2]);*/

#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("externalPlatineForces : The springP force is :%e,%e,%e\n",FPlatineSpringP[0],FPlatineSpringP[1],FPlatineSpringP[2]);

  printf("externalPlatineForces : The ext forces is :%e %e %e\n",f[0],f[1],f[2]);
#endif
  f[0]=-f[0];
  f[1]=-f[1];
  f[2]=-f[2];
  
  
}

/*/////////////////////////////////////////////////////////////////////////////////////////////////////////*/


extern "C" void externalPlatineMomentum(double t, double* q,double *v,double *m, unsigned int size_z,double *z){
  m[0]=0;
  m[1]=0;
  m[2]=0;

  ::boost::math::quaternion<double>    quattrf(q[3],q[4],q[5],q[6]);
  ::boost::math::quaternion<double>    cquattrf(q[3],-q[4],-q[5],-q[6]);
  ::boost::math::quaternion<double>    quatbuf;

  /*location of the force from A2.*/
  double G0[3];
  G0[0]=3.319088;
  G0[1]=4.365461;
  G0[2]=-4.706395;
  double P0A2[3];
  P0A2[0]=0.19;
  P0A2[1]=9;
  P0A2[2]=-7.6;
  double GP0A2[3] ;
  GP0A2[0]=P0A2[0]-G0[0];
  GP0A2[1]=P0A2[1]-G0[1];
  GP0A2[2]=P0A2[2]-G0[2];
  ::boost::math::quaternion<double>    quatGP0A2(0,GP0A2[0],GP0A2[1],GP0A2[2]);
  quatbuf=quattrf*quatGP0A2*cquattrf;
  double GPA2[3];
  double MA2[3];
  GPA2[0]=quatbuf.R_component_2();
  GPA2[1]=quatbuf.R_component_3();
  GPA2[2]=quatbuf.R_component_4();
  double f[3];
  FPlatineA2( t, q, v, f, size_z, z);
  
  cross3(GPA2,f,MA2);
/*  printf("externalPlatineMomentum : momemtum due to A2:%e,%e,%e\n",MA2[0],MA2[1],MA2[2]);*/

#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("externalPlatineMomentum : momemtum due to A2:%e,%e,%e\n",MA2[0],MA2[1],MA2[2]);
#endif

  
  m[0]=MA2[0];
  m[1]=MA2[1];
  double MZPlatineMomemtumDueSpring ;
  double Theta_equipage;
  double lva1r=0.2513;
  MZPlatineMomemtumDueSpringFct(&(*(sDS[EQUIPAGE]->q()))(0), q, &MZPlatineMomemtumDueSpring, &Theta_equipage);
    
  m[2]=MA2[2]+(lva1r*MZPlatineMomemtumDueSpring);
/*  printf("externalPlatineMomentum :%e,%e,%e\n",m[0],m[1],m[2]);*/
#ifdef  VERBOSE_DEBUG_PLUGIN  
  printf("externalPlatineMomentum :%e,%e,%e\n",m[0],m[1],m[2]);
#endif
  m[0]=-m[0];
  m[1]=-m[1];
  m[2]=-m[2];
}

/*////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
double FEquipageSpringP[3];
extern "C" void externalEquipageForces(double t, double* q,double *v,double *f, unsigned int size_z,double *z){
   f[0]=0;
   f[1]=0;
   f[2]=0;

   ::boost::math::quaternion<double>    quattrf1(q[3],q[4],q[5],q[6]);
   ::boost::math::quaternion<double>    cquattrf1(q[3],-q[4],-q[5],-q[6]);
   ::boost::math::quaternion<double>    quatbuf1;

  double M ;
  double Theta_equipage;
  MZPlatineMomemtumDueSpringFct(q,  &(*(sDS[PLATINE]->q()))(0)   ,  &M, &Theta_equipage);

    
  double d_platine=5.8;
  double d_equipage=6.6;
  double d_gplatine = 5.48;//5(old plugin value validated);/*5.48 new value*/
  double d_gequipage = 6.06;//4.84(old plugin value validated);/*7.99 new value*/
  //double Theta0_g_equi= M_PI+1.35; /* 1.22(old plugin value validated)*//*1.35 new value*/
  //double Theta_g_equip = Theta0_g_equi+Theta_equipage;
  //  double Theta0_g_platine = angleOrthoPG_X;
  //f[0]+=(M/d_gequipage)*cos(Theta_g_equip);
  //f[1]+=(M/d_gequipage)*sin(Theta_g_equip);
  /*printf("externalEquipageForces : The equivalent, spring A3 force is fg= %lf, theta_g_equip = %lf.\n",(M/d_gequipage),Theta_g_equip);*/
/*printf("externalEquipageForces : The equivalent, spring A3 force is fg= %lf, f_0 = %lf,f_1 = %lf,theta_g_equip = %lf.\n",(M/d_gequipage), f[0],f[1],Theta_g_equip);*/
#ifdef  VERBOSE_DEBUG_PLUGIN
  printf("externalEquipageForces : The equivalent, spring A3 force is fg= %lf, theta = %lf.\n",(M/d_gequipage),Theta_g_equip);
#endif

  ::boost::math::quaternion<double>    quatDir0_3(0,0.6,-1.1,-0.01);

  quatbuf1=quattrf1*quatDir0_3*cquattrf1;

  double norm2FE1=sqrt(quatbuf1.R_component_2()*quatbuf1.R_component_2()+quatbuf1.R_component_3()*quatbuf1.R_component_3()+quatbuf1.R_component_4()*quatbuf1.R_component_4());

   FEquipageSpringP[0]=quatbuf1.R_component_2()* (M/d_gequipage)/norm2FE1;
   FEquipageSpringP[1]=quatbuf1.R_component_3()* (M/d_gequipage)/norm2FE1;
   FEquipageSpringP[2]=quatbuf1.R_component_4()* (M/d_gequipage)/norm2FE1;
  
   f[0]+=FEquipageSpringP[0];
   f[1]+=FEquipageSpringP[1];
   f[2]+=FEquipageSpringP[2];

   f[0]=-f[0];
   f[1]=-f[1];
   f[2]=-f[2];
}

/*////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

extern "C" void externalBarreForces(double t, double* q,double *v,double *f, unsigned int size_z,double *z){
 f[0]=0;
 f[1]=0;
 f[2]=0;
}
extern "C" void externalBarreMomentum(double t, double* q,double *v,double *m, unsigned int size_z,double *z){
  if (t <1.0)
  {
    m[0]=0;
    m[1]=0;
    m[2]=15.0;
  }

  else
  {
    m[0]=0;
    m[1]=0;
    m[2]=0.0;
  }

//-3.0;/*FvD bigger clearance-5;23May2015 */
  //printf("externalBarreMomentums: q[6]=%lf \n",q[6]);
  return; //changed on date 08/05/2015
//   if (q[6]<0.349911) //  if (q[6]<0.349911) ...changed on date 08/05/2015
//     f[2]=100;
//   else
//     f[2]=-10;
 
// /*artificial momentum*/
//   /*spring white*/
//   double Theta_equipage=2*acos(sDS[EQUIPAGE]->q()->getValue(3));
//   if (sDS[EQUIPAGE]->q()->getValue(6)<0)
//     Theta_equipage=-Theta_equipage;
  
//   double Theta0_equipage=2*acos(sDS[EQUIPAGE]->q0()->getValue(3));
//   if (sDS[EQUIPAGE]->q0()->getValue(6)<0)
//     Theta0_equipage=-Theta0_equipage;
//   if (Theta_equipage>Theta0_equipage + 4*M_PI/180.0){
//     f[2]-=5;
//    }
//   else{
//     f[2]+=5;
//    }
// //printf("externalBarreMomentums: spring white available :f[2]=%lf. %lf \n",f[2],q[6]);
// #ifdef  VERBOSE_DEBUG_PLUGIN
//   printf("externalBarreMomentums: spring white available :f[2]=%lf. %lf \n",f[2],q[6]);
// #endif
}

extern "C" void externalMomentumx(double t, double* q,double *v,double *m, unsigned int size_z,double *z){

double angle =2*acos(q[3]);
//printf("angle:%e\n",angle);
 if (angle > 0.459)
  {
  	m[0]=0;
  	m[1]=0;
  	m[2]=-0;
  }	
 else
  {
    m[0]=0.0;
    m[1]=0;
    m[2]=0;
  }

 m[0]=-m[0];
 m[1]=-m[1];
 m[2]=-m[2];
 return; 
}

extern "C" void prescribedvelocityB1(double time, unsigned int sizeofprescribedvelocity, double *pv)
{
  /* the plugin implements v(t) = C + A cos(omega *t) */
   if (time <0.7)
    {
      double C = 0.0;
      pv[0] =  C;
    }
   else 
    {
      double C = 0.2;
      pv[0] =  C;
    }
    //printf("prescribed velocity = %e\n", pv[0]);
    
}
