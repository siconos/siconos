double CamState(double time, double rpm, double& CamPosition, double& CamVelocity, double& CamAcceleration)
{
  double t = time;
  //  double rpm = *param;
  double res = -0.0;
  double w, beta, hc, hcp, hcpp, mass, kspring, cfriction, phio;
  double PI = 3.14159265;
  double Rb = 0.40;
  double k1 = 0.45;
  double k2 = 0.40393320723821;
  double rho1 = 0.85;
  double rho3 = 0.15;
  double rho2 = 0.55393320723821;
  double dBetadt, tmp1;
  double beta1 = 0;
  double beta2 = 1.050328174371336;
  double beta3 = 1.22173047639603;


  // char str [80],comando[80];
  //   float f;
  //   FILE * pFile;

  //   pFile = fopen ("parameters.txt","r");
  //   rewind (pFile);
  //   fscanf (pFile, "%f", &f);
  //   fclose (pFile);
  // printf ("%f\n",f);
  //rpm=(double)(f);

  phio = PI / 2;
  //  rpm=358;

  //  printf ("%f\n",rpm);
  //  rpm=150;
  w = 2 * PI * rpm / 60;
  mass = 1.221;
  kspring = 1430.8;
  cfriction = 0;

  beta = w * t;
  beta = fmod(beta + phio, 2 * PI);

  //hc=sin(beta)-0.9;
  //hp=w*cos(beta);
  //hpp=-w*w*sin(beta);

  if (beta > PI)
  {
    beta = 2 * PI - beta;
    w = -w;
  }
  dBetadt = w;
  if (beta <= (PI / 2 - beta1))
  {
    hc = Rb;
    hcp = 0;
    hcpp = 0;
  }
  else if (beta <= (PI / 2 + beta2))
  {
    hc = -k1 * sin(beta) + rho1 * sqrt(1 - pow((cos(beta)) * (k1 / rho1), 2)) ;
    hcp = (-k1 * cos(beta) + (pow(k1, 2) * sin(beta) * cos(beta)) / (rho1 * sqrt(1 - pow((cos(beta)) * (k1 / rho1), 2)))) * dBetadt;
    tmp1 = pow(cos(beta - beta1), 2) - pow(sin(beta - beta1), 2) - (pow(k1 * sin(beta - beta1) * cos(beta - beta1), 2)) / (pow(rho1, 2) * (1 - pow((cos(beta - beta1)) * (k1 / rho1), 2)));
    hcpp = (k1 * sin(beta - beta1) + (pow(k1, 2) / (rho1 * sqrt(1 - pow((cos(beta - beta1)) * (k1 / rho1), 2)))) * tmp1) * (pow(dBetadt, 2));
  }
  else if (beta <= (PI / 2 + beta3))
  {
    hc = -k2 * cos(beta + PI / 9) + rho3 * sqrt(1 - pow((sin(beta + PI / 9)) * (k2 / rho3), 2)) ;
    hcp = (k2 * sin(beta + PI / 9) - (pow(k2, 2) * sin(beta + PI / 9) * cos(beta + PI / 9)) / (rho3 * sqrt(1 - pow((sin(beta + PI / 9)) * (k2 / rho3), 2)))) * dBetadt;
    tmp1 = pow(cos(beta - beta3), 2) - pow(sin(beta - beta3), 2) - (pow(k2 * sin(beta - beta3) * cos(beta - beta3), 2)) / (pow(rho3, 2) * (1 - pow((cos(beta - beta3)) * (k2 / rho3), 2)));
    hcpp = (-k2 * sin(beta - beta3) + (pow(k2, 2) / (rho3 * sqrt(1 - pow((cos(beta - beta3)) * (k2 / rho3), 2)))) * tmp1) * pow(dBetadt, 2);
  }
  else
  {
    hc = rho2;
    hcp = 0;
    hcpp = 0;
  }

  hc = hc - 0.3;

  CamPosition = hc;
  CamVelocity = hcp;
  CamAcceleration = hcpp;
  res = -(mass * hcpp + cfriction * hcp + kspring * hc);
  //res=0;
  return res;
}
