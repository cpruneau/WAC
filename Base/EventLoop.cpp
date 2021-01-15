// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */

#include "EventLoop.hpp"
ClassImp(EventLoop);

EventLoop::EventLoop(const TString & name)
:
TaskCollection(name,nullptr,100),
timer(),
nEventRequested(1000000),
nEventReported(10000),
nEventPartialSave(10000),
nEventProcessed(0),
partialSave(0),
subsampleAnalysis(0)
{
  timer.start();
}

EventLoop::~EventLoop()
{
  // no ops
}

void EventLoop::run()
{
  run(nEventRequested,nEventReported);
}

void EventLoop::run(long nEvent, long nReport)
{
  if (reportInfo("EventLoop",getTaskName(),"run(...)")) cout << "Running for nEvent: " << nEvent << endl;
  postTaskOk();
  initialize();
  if (!isTaskOk())
    {
    if (reportWarning("EventLoop",getTaskName(),"run(...)")) cout << "Initialization failed. Abort." << endl;
    return;
    }
  nEventProcessed = 0;
  if (reportInfo("EventLoop",getTaskName(),"run(...)")) cout << "Starting..." << endl;
  for (long iEvent=1; iEvent<=nEvent; iEvent++)
    {
    execute(); if (!isTaskOk()) continue;
    nEventProcessed++;
    if (nEventProcessed%nReport==0 )
      {
      if (reportInfo("EventLoop",getTaskName(),"run(...)")) cout << "Completed event # " << iEvent << endl;
      }
//    cout << "nEventProcessed" << nEventProcessed<< endl;
//    cout << "nEventPartialSave" << nEventPartialSave<< endl;
//    cout << nEventProcessed%nEventPartialSave<<endl;
    if ( (subsampleAnalysis||partialSave) && nEventProcessed%nEventPartialSave==0)
      {
      savePartialResults();
      if (subsampleAnalysis) reset();
      if (!isTaskOk()) continue;
      }
    }
  if (isTaskOk() && !((subsampleAnalysis||partialSave) && nEventProcessed%nEventPartialSave==0)) finalize(); //if (subsampleAnalysis||partialSave) && nEventProcessed%nEventPartialSave==0  is true then this already happens above in savePartialResults()
  timer.stop();
  if (reportInfo("EventLoop",getTaskName(),"run(...)"))
    {
    cout << endl;
    cout << "  Completed with status : " << getTaskStatusName() << endl;
    cout << "        Completed Events: " << nEventProcessed << endl;
    cout << "        Requested Events: " << nEvent << endl;
    cout << "            "; timer.print(cout);
    cout << endl << endl;
    }
}

