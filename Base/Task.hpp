// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_Task
#define WAC_Task
#include <iostream>
#include "TClass.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TList.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom.h"
#include "TString.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TFile.h"
#include "TaskConfiguration.hpp"
#include "MessageLogger.hpp"
#include "Event.hpp"

using namespace std;

class Task : public MessageLogger
{
public:



  ////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////
  TString             taskName;
  TaskConfiguration * taskConfiguration;
  TRandom           * taskRandomGenerator;
  Event             * event;

  long nEventProcessed;
  long nEventAccepted;
  int  subSampleIndex;

  ////////////////////////////////////////////
  // Member functions
  ////////////////////////////////////////////

  Task(const TString & name,
       TaskConfiguration * configuration,
       Event * event,
       LogLevel selectedLevel);
  virtual ~Task();
  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual void reset();
  virtual void clear();
  virtual void savePartialResults();
  
  virtual void createHistograms();
  virtual void loadHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void resetHistograms();
  virtual void clearHistograms();
  virtual void calculateDerivedHistograms();
  virtual void scaleHistograms();
  virtual void scaleHistograms(double factor);
  virtual void saveHistogramsAsText();
  virtual void saveHistograms();
  virtual void saveHistograms(TFile * outputFile);
  virtual void addHistogramsToExtList(TList *list, bool all=false);

  virtual double readParameter(TFile * inputFile, const TString & parameterName);
  TFile * openRootFile(const TString & inputPath, const TString & fileName, const TString & ioOption);

  Event * getEvent();

  TaskConfiguration * getTaskConfiguration();
  void setTaskConfiguration(TaskConfiguration * config);
  virtual void printTaskConfiguration(ostream & output);
  inline TString getTaskName() const
  {
    return taskName;
  }

  inline void setTaskName(const TString & name)
  {
    taskName = name;
  }

  TRandom * getRandomGenerator()
  {
  return taskRandomGenerator;
  }

  inline long incrementEventProcessed()
  {
  return ++nEventProcessed;
  }

  inline long incrementEventAccepted()
  {
  return ++nEventAccepted;
  }

  inline long getNEventProcessed() const
  {
  return nEventProcessed;
  }

  inline long getNEventAccepted() const
  {
  return nEventAccepted;
  }

  inline long getSubSampleIndex() const
  {
  return subSampleIndex;
  }

  void saveNEventProcessed(TFile * outputFile);
  void saveNEventAccepted(TFile * outputFile);
  long loadNEventProcessed(TFile * inputFile);
  long loadNEventAccepted(TFile * inputFile);


  void setRandomGenerator(TRandom * randomGenerator);


  enum TaskStatus   { Unknown, TaskOk, TaskEof, TaskEod, TaskWarning, TaskError, TaskFatal};

private:

  static TaskStatus taskStatus;

public:

  inline static TaskStatus getTaskStatus()
  {
    return taskStatus;
  }

  inline static void setTaskStatus(TaskStatus newStatus)
  {
    taskStatus = newStatus;
  }

  inline static void postTaskOk()
  {
    taskStatus = TaskOk;
  }

  inline static void postTaskEof()
  {
    taskStatus = TaskEof;
  }


  inline static void postTaskEod()
  {
    taskStatus = TaskEod;
  }

  inline static void postTaskWarning()
  {
    taskStatus = TaskWarning;
  }

  inline static void postTaskError()
  {
    taskStatus = TaskError;
  }

  inline static void postTaskFatal()
  {
    taskStatus = TaskFatal;
  }


  inline static bool isTaskOk()
  {
    return (taskStatus == TaskOk);
  }

  static TString getTaskStatusName();

  ClassDef(Task,0)
};

#endif /* WAC_Task */
