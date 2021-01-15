// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_Timer
#define WAC_Timer
#include <iostream>
#include <chrono>
#include <fstream>
#include "TObject.h"

using namespace std;
using namespace std::chrono;

class Timer
{
public:

  Timer();
  virtual ~Timer();
  void start();
  void stop();
  void print(ostream & os);

  high_resolution_clock::time_point startTime;
  high_resolution_clock::time_point stopTime;
  duration<double> totalDuration;
  int    days;
  int    hours;
  int    minutes;
  double seconds;

  ClassDef(Timer,0)
};

#endif /* WAC_Timer */
