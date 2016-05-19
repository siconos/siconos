/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

// for M_PI
#if defined(_MSC_VER)
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include <assert.h>

double puissance(double a, int n)
{

  double resultat = 1 ;

  for (int i = 0 ; i < n ; i++)
    resultat = resultat * a ;
  return resultat ;
}

double r = 0.005;    // rayons des cylindres du yoyo
double R = 0.03;
const double m = 0.1;   // masse du yoyo
double I = 0.5 * m * R * R; // moment d'inertie du yoyo par rapport à l'axe de rotation
double L = 0.4; // longueur totale du fil
const double g = 9.8; // valeur de champ de gravité
const double epsilon = 0.000; // coefficient de frottement entre le tambour et le fil

double e_eq = (I - m*r*r) / (I + m*r*r); // coefficient de restitution équivalent
double ks = (1 - e_eq) * (3 + e_eq*e_eq) / puissance(1 + e_eq, 3);
double w = ks + 0.005; // coefficient d'accélération
double k0 = (1 - e_eq) / (1 + e_eq);
double c1 = 10;
double c2 = 60 ; //coefficients d'amortissement
double ue = pow(e_eq * w + e_eq - 1 + w, 2) / ((pow(1 + e_eq, 2) * w + pow(1 - e_eq, 2)) * w) ;

double A[] = {1.369, -3.665, 1.748, -0.4906, 0.01883};
double Cy = sqrt(2 * L * (I + m*r*r) / (m*r*r*g));

//double hand(double time){return 0;} double deriveehand(double time) {return 0;}double deriveederiveehand(double time) {return 0;}


double accelerationmain(int n , double A[] , double Cy, double time)
{
  double resultat = 0 ;
  for (int i = 0; i < n ; i++)
    resultat = resultat + A[i] * sin((i + 1) * (M_PI / Cy) * time);
  return resultat;

}

double vitessemain(int n , double A[] , double Cy, double time)
{
  double resultat = 0 ;
  for (int i = 0; i < n ; i++)
    resultat = resultat - A[i] * cos((i + 1) * (M_PI / Cy) * time) / ((i + 1) * M_PI / Cy);
  return resultat;
}

// amplitudes objectifs  avec les instants corrependants
int const G = 5;
double temps[G] = {5, 20, 35, 50, 60} ;
double Som[G] = {L / 2, L / 4, L, L, L / 2};



double thetaset(double Som[], size_t i)
{

  assert(i < G);
  return (1 - ue) * (Som[i]) / r ;

}
