#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <math.h>


extern "C"   double computeControl1(double time)
{


  double u;
  double alpha = 10000.0;
  double T = 0.0001;
  int oddoreven = 1 ;

  oddoreven = int(time / T);
  if ((oddoreven / 2) == 0) u = alpha;
  else u = -alpha;
  u = 30 * sin(50 * time);

  return u;
}



extern "C"   double computeControl2(double time)
{


  double u;
  double alpha = -50.0;
  int njump;

  if (time <= 1) u = -0.0;
  else if (time >= 3) u = -alpha;
  else // (time > 1 && time <= 3)
  {
    njump = (int)((log(1.0 / (3.0 - time)) / log(2.0)) + 1.0);


    //printf("njump = %i\n",njump);
    if (njump % 2 == 0)
    {
      u = -alpha;
    }
    else
    {
      u = alpha;
    }
  };


  return u;
}
extern "C"   double computeControl(double time)
{


  double u;
  double alpha = 50.0;
  int oddoreven = 1 ;
  int njump;
  double a = 24 / 25.0;

  double timeaccu = 1 / (1.0 - a);


  if (time < 1) u = 0.0;

  else if ((time >= 1) && (time < timeaccu))
  {
    njump = (int)((log(1.0 - (1.0 - a) * time) / log(a)) - 1.0);
    u =  alpha / (pow(2 * njump + 1.0, 1.0 / a));


    if ((njump % 2) == 0) u = -u;

    //    printf("njump = %i\n",njump);
    //printf("time = %e\n",time);
    //u =  -alpha*(1.0+pow(2,njump+1)*(3.0-1.0/(pow(2,njump-1))));
    //printf("u = %e\n",u);
  }
  else // (time >= timeaccu)
  {
    oddoreven = int(time - timeaccu);
    printf("time = %e\n", time);
    printf("oddorven = %i\n", oddoreven);
    if ((oddoreven % 2) == 0) u = alpha / 10;
    else u = -alpha / 10;
    printf("u = %e\n", u);
  }



  return u;
}
extern "C"   double computeControlori2(double time)
{


  double u;
  double alpha = 50.0;
  int oddoreven = 1 ;
  int njump;
  int N = 100;
  if (time < 1) u = -0.0;

  else if ((time >= 1) && (time < 2))
  {
    njump = 0;
    u =  alpha * pow(2.0, njump + 1) * time - alpha * (1.0 + pow(2.0, njump + 1)
         * (3.0 -
            pow(2.0, -njump + 1)
           )
                                                    ) ;
    //printf("njump = %i\n",njump);
    //printf("time = %e\n",time);
    //u =  -alpha*(1.0+pow(2,njump+1)*(3.0-1.0/(pow(2,njump-1))));
    //printf("u = %e\n",u);
  }

  else if ((time >= 2) && (time < 3 - pow(2.0, -N)))
  {
    njump = (int)((log(1.0 / (3.0 - time)) / log(2.0)) + 1.0);
    u =  alpha * pow(2.0, njump + 1) * time - alpha * (1.0 + pow(2.0, njump + 1)
         * (3.0 -
            pow(2.0, -njump + 1)
           )
                                                    ) ;
    printf("njump = %i\n", njump);
    //printf("time = %e\n",time);
    //u =  -alpha*(1.0+pow(2,njump+1)*(3.0-1.0/(pow(2,njump-1))));
    //printf("u = %e\n",u);
  }
  else if ((time >= 3 - pow(2.0, -N)) && (time < 3))
  {
    u = 0.0;
  }
  else // (time >= 3)
  {
    oddoreven = int(time - 3);
    printf("time = %e\n", time);
    printf("oddorven = %i\n", oddoreven);
    if ((oddoreven % 2) == 0) u = alpha;
    else u = -alpha;
    printf("u = %e\n", u);
  }



  return u;
}
extern "C"   double computeControlori(double time)
{


  double u;
  double alpha = 50.0;
  int njump;

  if (time < 1) u = -0.0;
  else if (time >= 3) u = 0.0;
  else if ((time >= 1) && (time < 2))
  {
    njump = 0;
    u =  alpha * pow(2.0, njump + 1) * time - alpha * (1.0 + pow(2.0, njump + 1)
         * (3.0 -
            pow(2.0, -njump + 1)
           )
                                                    ) ;
    //printf("njump = %i\n",njump);
    //printf("time = %e\n",time);
    //u =  -alpha*(1.0+pow(2,njump+1)*(3.0-1.0/(pow(2,njump-1))));
    //printf("u = %e\n",u);
  }

  else // (time >= 2 && time < 3)
  {
    njump = (int)((log(1.0 / (3.0 - time)) / log(2.0)) + 1.0);
    u =  alpha * pow(2.0, njump + 1) * time - alpha * (1.0 + pow(2.0, njump + 1)
         * (3.0 -
            pow(2.0, -njump + 1)
           )
                                                    ) ;
    //printf("njump = %i\n",njump);
    //printf("time = %e\n",time);
    //u =  -alpha*(1.0+pow(2,njump+1)*(3.0-1.0/(pow(2,njump-1))));
    //printf("u = %e\n",u);
  };


  return u;
}




SICONOS_EXPORT void computeU(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{

  double u;
  u = computeControl(time);
  b[0] = 0.1 * u ;
  b[1] = 0.2 * u;
  b[2] = 0.1 * u ;
  b[3] = 0.2 * u;

}


SICONOS_EXPORT void computeE(double time, unsigned int sizeOfB, double* e, unsigned int sizeOfZ, double* z)
{

  // Si b est de taille 2
  e[0] = computeControl(time) ;
  e[1] = computeControl(time) ;
}
