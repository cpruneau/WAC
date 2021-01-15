// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Timer
 \ingroup WAC

 Class defining Timer
 */

#include "Timer.hpp"

ClassImp(Timer);

Timer::Timer()
:
startTime(),
stopTime(),
totalDuration(),
hours(0),
minutes(0),
seconds(0)
{
  start();
}

Timer::~Timer()
{
  // no ops.
}

void Timer::start()
{
  startTime = high_resolution_clock::now();
}
void Timer::stop()
{
  stopTime = high_resolution_clock::now();
  totalDuration = chrono::duration_cast<chrono::seconds>(stopTime - startTime);
  days     = (int)(totalDuration.count()/(24*3600));
  hours    = (int)(totalDuration.count()/3600);
  minutes  = (int)((totalDuration.count() - 3600 * hours)/60);
  seconds  = totalDuration.count() - 60 * minutes - 3600 * hours;
}

void Timer::print(ostream & os)
{
  os << "Time elapsed: " << days << " days, "<< hours << " hours, " << minutes << " minutes, " << seconds << " seconds." <<endl;
}


