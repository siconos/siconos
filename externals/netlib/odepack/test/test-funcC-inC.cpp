/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <stdio.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "SiconosFortran.h"

void f1(int* neq, double* t, double* y, double* ydot) {
  ydot[0] = ((2.0 * std::log(y[0]) + 8.0) / *t - 5.0) * y[0];
}

void gr1(int* neq, double* t, double* y, int* ng, double* groot) {
  groot[0] = ((2.0 * std::log(y[0]) + 8.0) / *t - 5.0) * y[0];
  groot[1] = log(y[0]) - 2.2491;
}

void jdum(int* neq, double* t, double* y, int* ml, int* mu, double* pd, int* nrowpd) {}

int main(void) {
  // We reproduce in C++ the first example described in DLSODAR-test.f.
  int neq = 1;
  std::vector<double> y = {1., 2, 3};
  double t = 1.;
  double tout = 2.;
  int itol = 1;
  double rtol = 1e-6;
  std::vector<double> atol = {rtol, 0.};
  int istate = 1;
  int iopt = 0;
  int lrw = 44;
  int liw = 21;
  int jt = 2;
  int ng = 2;
  std::vector<int> jroot(ng);
  std::vector<double> rwork(lrw);
  std::vector<int> iwork(liw);
  int nerr = 0;
  double ero = 0.;
  double errt = 0;
  std::cout << "Demonstration program for DLSODAR package \n\n\n\n"
            << "First problem\n\n\n"
            << "Problem is  dy/dt = ((2*log(y)+8)/t - 5)*y,  y(1) = 1\n\n"
            << " Solution is  y(t) = exp(-t**2 + 5*t - 4)\n\n"
            << " Root functions are: \n"
            << " g1 = dy/dt  (root at t = 2.5) \n"
            << " g2 = log(y) - 2.2491  (roots at t = 2.47 and t = 2.53)\n\n"
            << " itol =" << itol << "  rtol = " << rtol << "   atol =" << atol[0]
            << " \n\njt =" << jt << "\n\n\n";

  int iter = 0;
  while (iter < 5) {
    siconos::netlib::lsodar(&f1, &neq, &y.front(), &t, &tout, &itol, &rtol, &atol.front(),
                            &istate, &rwork.front(), &lrw, &iwork.front(), &liw, &jdum, &jt,
                            &gr1, &ng, &jroot.front());
    auto yt = std::exp(-t * t + 5. * t - 4.);
    auto er = y[0] - yt;
    std::cout << " At t =" << t << " y = " << y[0] << " error = " << er << "\n";
    if (istate < 0) {
      nerr += 1;
      auto nfea = iwork[10];
      if (jt == 2) nfea = iwork[11] - neq * iwork[12];

      std::cout << "\n\n\nFinal statistics for this run:\n";
      std::cout << "rwork size " << iwork[16] << "iwork size " << iwork[17]
                << "\n number of steps " << iwork[10] << "\n number of f-s" << iwork[11]
                << "\n (excluding j-s)" << nfea << "\n number of j-s" << iwork[12]
                << "\n number of g-s" << iwork[9] << "\n error overrun" << ero << "\n";

      break;
    }
    er = std::abs(er) / (rtol * abs(y[0]) + atol[0]);
    ero = std::max(ero, er);
    if (er > 1000) {
      nerr += 1;
      std::cout << "Warning: error exceeds 1000 * tolerance.\n";
    }
    if (istate != 3) {
      tout += 1.;
      iter += 1;
      continue;
    }
    // If a root was found, write results and check root location.
    // Then reset istate to 2 and return to DLSODAR call.
    std::cout << "\n Root found at t =" << t << " jroot=" << jroot[0] << "\t" << jroot[1]
              << "\n";
    if (jroot[0] == 1) errt = t - 2.5;
    if (jroot[1] == 1 && t <= 2.5) errt = t - 2.47;
    if (jroot[1] == 1 && t > 2.5) errt = t - 2.53;
    std::cout << "Error in t location of root is " << errt << "\n\n";
    if (std::abs(errt) > 1.e-3) {
      std::cout << " Warning: root error exceeds 1e-3\n";
      nerr += 1;
    }
    istate = 2;
  }
  if (istate < 0) nerr += 1;
  auto nfea = iwork[10];
  if (jt == 2) nfea = iwork[11] - neq * iwork[12];
  std::cout << "\n\n\nFinal statistics for this run:\n";
  std::cout << "rwork size " << iwork[16] << "  iwork size " << iwork[17]
            << "\nnumber of steps " << iwork[10] << "\nnumber of f-s " << iwork[11]
            << "\n(excluding j-s) " << nfea << "\nnumber of j-s " << iwork[12]
            << "\nnumber of g-s " << iwork[9] << "\nerror overrun " << ero << "\n";
}
