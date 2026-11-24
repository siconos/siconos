// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/************************************************************************
                        acetime.h

**************************************************************************/
#include <stdio.h>

#include <iostream>

#include "sys/time.h"
using namespace std;

#ifndef ACETIME_H
#define ACETIME_H
// Class aceTime
//
//

/** \class aceTime
    \brief simple timer utilities
 */
class aceTime {
 public:
  aceTime();
  void start();
  void stop();
  void setName(char const* Name);
  void print(ostream& os = cout);
  virtual ~aceTime();

 protected:
 private:
  long mCall;
  timeval mStart;
  long mCumul;
  char mName[128];
  bool mIsRunning;
};
#endif  // ACETIME_H
