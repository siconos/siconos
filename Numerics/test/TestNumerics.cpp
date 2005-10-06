/* Siconos version 1.0, Copyright INRIA 2005.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "TestNumerics.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(TestNumerics);


TestNumerics::TestNumerics()
{}
TestNumerics::~TestNumerics()
{}

void TestNumerics::setUp()
{}

void TestNumerics::tearDown()
{}

//______________________________________________________________________________

/*
 * tests solverpack
 */
void TestNumerics::testLcp()
{
  int i = system("cd ../src/solverpack/test/; ./test_lcp");
  CPPUNIT_ASSERT_MESSAGE("testLcp ", i != -1);
  printf("TestNumerics >>> testLcp ................................ OK\n");
}

void TestNumerics::testLcp2()
{
  //  //FILE *f1,*f2;
  //  ifstream f1, f2;
  //  char *tmp;
  //  int i,j,nl,nc,nll,it,info,n=31,dimM=n;
  //  double *q,*z,*w,*vec;
  //  //double (*M)[n];
  //  double **M;
  //  double qi,Mij;
  //  char val[14],vall[14];
  //
  //  methode meth;
  //  static methode_lcp meth_lcp  = {"gcp",101, 0.0001,0.6};
  //
  //  meth.lcp = meth_lcp;
  //
  //
  //  f1.open("MM_mmc.dat", ios::in);    // open the streams
  //  f2.open("qq_mmc.dat", ios::in);
  ////  if ((f1=fopen("MM_mmc.dat","r"))==NULL){
  ////  perror("fopen 1");
  ////  exit(1);
  ////  }
  ////
  ////  if ((f2=fopen("qq_mmc.dat","r"))==NULL){
  ////  perror("fopen 2");
  ////  exit(2);
  ////  }
  //
  //  //  f3=fopen("aff_M.dat","w+");
  //
  //  //M=malloc(dimM*dimM*sizeof(double));
  //  *M = new double[dimM];
  //  for(i=0; i<dimM; i++) M[i] = new double [dimM];
  //
  //  //vec=(double*)malloc(dimM*dimM*sizeof(double));
  //  vec = new double [dimM*dimM];
  //
  //  for (i=0;i<dimM;i++)
  //    for (j=0;j<dimM;j++)
  //    M[i][j]=0.;
  //
  //  while (!f1.eof()/*feof(f1)*/)
  //  {
  ////    fscanf(f1,"%d",&nl);
  ////    fscanf(f1,"%d",&nc);
  ////    fscanf(f1,"%s",val);
  //  f1.read(tmp, sizeof(int));
  //  sprintf(tmp, "%d", nl);
  //  f1.read(tmp, sizeof(int));
  //  sprintf(tmp, "%d", nc);
  //  f1.read(tmp, sizeof(int));
  //  sprintf(tmp, "%f", Mij);
  ////    Mij=atof(val);
  //
  //  /////////////       on met la transpos       ////////////////
  //    /*fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
  //      fscanf(f3,"%d %d %.14e\n", &nc, &nl, &Mij);*/
  //    *(*(M+nc-1)+nl-1)=Mij;
  //  //////////////         fin transpos         ////////////////////
  //   }
  //
  //
  //  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  //    for (i=0;i<dimM;i++)
  //      for (j=0;j<dimM;j++)
  //  vec[j*dimM+i]= M[i][j];
  //  //       printf("vec(%d) = %.14e \n",i*dimM+j,M[i][j]);}
  //  ////////////////////////////////////////////////////////////////////////
  //
  //
  //
  ////  if ((f2=fopen("qq_mmc.dat","r"))==NULL){
  ////    perror("fopen 2");
  ////    exit(2);
  ////  }
  //
  //
  //  //  f4=fopen("aff_q.dat","w+");
  //  q=new double[dimM];//malloc(dimM*sizeof(double));
  //  z=new double[dimM];//malloc(dimM*sizeof(double));
  //  w=new double[dimM];//malloc(dimM*sizeof(double));
  //
  //  while (!f2.eof()/*feof(f2)*/){
  ////    fscanf(f2,"%d",&nll);
  ////    fscanf(f2,"%s",vall);
  //  f2.read(tmp, sizeof(int));
  //  sprintf(tmp, "%d", nll);
  //  f2.read(tmp, sizeof(int));
  //  sprintf(tmp, "%f", qi);
  //    //qi=atof(vall);
  //    //fprintf(f4,"%d %.14e\n",nll,qi);
  //    *(q+nll-1)=-qi;
  // }
  //f1.close();   // close the streams
  //  f2.close();
  //
  //  printf("\n we go in the function\n");
  //
  //  info=solve_lcp(vec, q, &n, /*&meth_lcp*/ &meth, z, w);
  //
  //  printf("\n we go out the function and info is %d\n",info);
  //
  //  //fclose(f2);fclose(f1);
  //
  //
  //  //free(M);  free(vec);  free(q);  free(z);  free(w);
  //  for(i=0; i<dimM; i++) delete []M[i];
  //  delete []M;
  //
  //  delete []q;
  //  delete []vec;
  //  delete []z;
  //  delete []w;
}

void TestNumerics::testRp()
{
  int i = system("cd ../src/solverpack/test/; ./test_rp");
  CPPUNIT_ASSERT_MESSAGE("testRp ", i != -1);
  printf("TestNumerics >>> testRp ................................ OK\n");
}

void TestNumerics::testCfd()
{
  int i = system("cd ../src/solverpack/test/; ./test_cfd");
  CPPUNIT_ASSERT_MESSAGE("testCfd ", i != -1);
  printf("TestNumerics >>> testCfd ................................ OK\n");
}

void TestNumerics::testCfp()
{
  int i = system("cd ../src/solverpack/test/; ./test_cfp");
  CPPUNIT_ASSERT_MESSAGE("testCfp ", i != -1);
  printf("TestNumerics >>> testCfp ................................ OK\n");
}

/*
 * tests ODEpack
 */
void TestNumerics::testDLSODE()
{
  int i = system("cd ../src/odepack/test/; ./DLSODE-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODE ", i != -1);
  printf("TestNumerics >>> testDLSODE ................................ OK\n");
}

void TestNumerics::testDLSODES()
{
  int i = system("cd ../src/odepack/test/; ./DLSODES-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODES ", i != -1);
  printf("TestNumerics >>> testDLSODES ................................ OK\n");
}

void TestNumerics::testDLSODA()
{
  int i = system("cd ../src/odepack/test/; ./DLSODA-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODA ", i != -1);
  printf("TestNumerics >>> testDLSODA ................................ OK\n");
}

void TestNumerics::testDLSODAR()
{
  int i = system("cd ../src/odepack/test/; ./DLSODAR-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODAR ", i != -1);
  printf("TestNumerics >>> testDLSODAR ................................ OK\n");
}

void TestNumerics::testDLSODPK()
{
  int i = system("cd ../src/odepack/test/; ./DLSODPK-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODPK ", i != -1);
  printf("TestNumerics >>> testDLSODPK ................................ OK\n");
}

void TestNumerics::testDLSODKR()
{
  int i = system("cd ../src/odepack/test/; ./DLSODKR-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODKR ", i != -1);
  printf("TestNumerics >>> testDLSODKR ................................ OK\n");
}

void TestNumerics::testDLSODI()
{
  int i = system("cd ../src/odepack/test/; ./DLSODI-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODI ", i != -1);
  printf("TestNumerics >>> testDLSODI ................................ OK\n");
}

void TestNumerics::testDLSOIBT()
{
  int i = system("cd ../src/odepack/test/; ./DLSOIBT-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSOIBT ", i != -1);
  printf("TestNumerics >>> testDLSOIBT ................................ OK\n");
}

void TestNumerics::testDLSODIS()
{
  int i = system("cd ../src/odepack/test/; ./DLSODIS-test");
  CPPUNIT_ASSERT_MESSAGE("testDLSODIS ", i != -1);
  printf("TestNumerics >>> testDLSODIS ................................ OK\n");
}
