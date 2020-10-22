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

#include "Task.hpp"

ClassImp(Task);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Task::Task(const TString & name,
           TaskConfiguration * configuration,
           Event * selectedEvent)
:
taskName             ( name ),
taskConfiguration    ( configuration ),
taskRandomGenerator  ( gRandom),
event                ( selectedEvent ),
eventsProcessed      ( 0 )
{
  if (reportDebug()) cout << "Task::Task() No ops." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Task::~Task()
{
  if (reportDebug()) cout << "Task::~Task() No ops." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Task::initialize()
{
  if (reportDebug()) cout << "Task::initialize() Started." << endl;
  if (isTaskOk() && taskConfiguration->loadHistograms)   loadHistograms();
  if (isTaskOk() && taskConfiguration->createHistograms) createHistograms();
  if (reportDebug()) cout << "Task::initialize() Completed." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Execute task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Task::execute()
{
  if (reportDebug()) cout << "Task::execute() No ops." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Finalize task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Task::finalize()
{
  if (reportDebug()) cout << "Task::finalize() Started task " << getName() << endl;
  if (isTaskOk() && taskConfiguration->scaleHistograms)            scaleHistograms();
  if (isTaskOk() && taskConfiguration->calculateDerivedHistograms) calculateDerivedHistograms();
  if (isTaskOk() && taskConfiguration->saveHistograms)             saveHistograms();
  if (reportDebug()) cout << "Task::finalize() Completed task " << getName()  << endl;
}

//  if (analysisConfiguration->scaleHistograms)
//    {
//    cout << "Task::finalize() Scale  histograms by 1/nEventsAccepted with nEventsAccepted = " << nEventsAccepted << endl;
//    double scalingFactor = 1.0/nEventsAccepted;
//    if (analysisConfiguration->externalScaling != 1.0)
//      {
//      if (reportInfo()) cout << "Task::finalize() Additionally scale  histograms by external factor = " << analysisConfiguration->externalScaling << endl;
//      scalingFactor *= analysisConfiguration->externalScaling;
//      }
//    scale(scalingFactor);
//    /* now be consistent with the stored information */
//    nEventsAccepted = 1;
//    }
//  if (analysisConfiguration->calculateDerived) calculateDerivedHistograms();
//  if (analysisConfiguration->saveHistosToRootFile) saveHistograms(analysisConfiguration->outputPath + analysisConfiguration->rootOutputFileName, analysisConfiguration->forceRewrite);
//  if (analysisConfiguration->saveHistosToTextFile) saveHistogramsAsText(); //analysisConfiguration->outputPath + "histograms.txt");
//}
//else
//{
//  if (reportWarning()) cout  << "Task::finalize() Cannot proceed because nEventsAccepted = " << nEventsAccepted << endl;
//    }
//    if (reportDebug())  cout << "Task::finalize() Completed" << endl;
//    }

void Task::reset()
{
  if (reportDebug()) cout << "Task::reset() Started." << endl;
  if (isTaskOk() && taskConfiguration->resetHistograms)   resetHistograms();
  if (reportDebug()) cout << "Task::reset() Completed." << endl;
}

////////////////////////////////////////////
// Reset this task
////////////////////////////////////////////
void Task::clear()
{
  if (reportDebug()) cout << "Task::reset() Started." << endl;
  if (isTaskOk() && taskConfiguration->clearHistograms)   clearHistograms();
  if (reportDebug()) cout << "Task::reset() Completed." << endl;
}

TaskConfiguration * Task::getTaskConfiguration()
{
  return taskConfiguration;
}

void Task::setTaskConfiguration(TaskConfiguration * config)
{
  taskConfiguration = config;
}

////////////////////////////////////////////
// Print task configuration
////////////////////////////////////////////
void Task::printConfiguration(ostream & output)
{
  if (reportInfo(output))
    {
    output << "Task Name : " << taskName <<  endl;
    taskConfiguration->printTaskConfiguration(output);
    }
}

TString Task::getName() const
{
  return taskName;
}

void Task::setName(const TString & name)
{
  taskName = name;
}


void Task::setRandomGenerator(TRandom * randomGenerator)
{
  if (!randomGenerator)
    taskRandomGenerator = randomGenerator;
  else
    {
    if (reportError()) cout << "Task:: setRandomGenerator(...) Null pointer. Random generator will not be set.";
    }
}

void Task::createHistograms()
{
  if (reportDebug()) cout << "Task::createHistograms() No ops" << endl;
}

void Task::loadHistograms()
{
  if (reportDebug()) cout << "Task::loadHistograms() Started." << endl;
  TFile * inputFile;
  TString inputFileName = taskConfiguration->getInputRootFileName();
  if (reportInfo()) cout << "Task::loadHistograms() Opening root input file " << inputFileName << endl;
  inputFile = new TFile(inputFileName,"OLD");
  if (!inputFile)
     {
     if (reportError()) cout << "Task::loadHistograms() Could not open root input file   " << inputFileName << endl;
     postTaskError();
     return;
    }
  loadHistograms(inputFile);
  if (reportDebug()) cout << "Task::loadHistograms() Completed." << endl;
}

void Task::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportError("Task",getName(),"loadHistograms(TFile * inputFile)")) cout << "Null pointer supplied." << endl;
    }
}

void Task::resetHistograms()
{
 if (reportDebug()) cout << "Task::resetHistograms() No ops." << endl;
}

void Task::clearHistograms()
{
 if (reportDebug()) cout << "Task::loadHistograms() No ops." << endl;
}

void Task::scaleHistograms()
{
  if (reportDebug()) cout << "Task::scaleHistograms() Started." << endl;
  if (reportInfo())  cout << "Processed number of events: " <<  eventsProcessed << endl;
  double scalingFactor = 1.0/double(eventsProcessed);
  scaleHistograms(scalingFactor);
if (reportDebug()) cout << "Task::scaleHistograms() Completed." << endl;
}

void Task::scaleHistograms(double factor)
{
  if (reportDebug()) cout << "Task::scaleHistograms(double factor) factor:" << factor << endl;
}

void Task::saveHistogramsAsText()
{
  if (reportDebug()) cout << "Task::saveHistogramsAsText() no ops:"  << endl;
}

void Task::saveHistograms()
{
  if (reportDebug()) cout << "Task::saveHistograms() Started for " << getName() << endl;
  TFile * outputFile;
  TString outputFileName = taskConfiguration->outputPath;
  outputFileName += taskConfiguration->rootOuputFileName;
  outputFileName += getName();
  outputFileName += ".root";

  if (taskConfiguration->forceHistogramsRewrite)
    {
    if (reportInfo()) cout << "Task::saveHistograms()  Opening root output (RECREATE) file  " << outputFileName << endl;
    outputFile = new TFile(outputFileName,"RECREATE"); // obliterate past work...
    if (!outputFile)
      {
      if (reportError()) cout << "Task::saveHistograms(outputFileName)  Could not open (RECREATE) file  " << outputFileName << endl;
      return;
      }
    }
  else
    {
      if (reportInfo()) cout << "Task::saveHistograms()  Opening root output (NEW) file  " << outputFileName << endl;
        outputFile = new TFile(outputFileName,"NEW"); // protect past work...
    if (!outputFile)
      {
      if (reportError()) cout << "Task::saveHistograms(outputFileName)  Could not open (NEW) file  " << outputFileName << endl;
      return;
      }
    }
  saveHistograms(outputFile);
  outputFile->Close();
 if (reportDebug()) cout << "Task::saveHistograms() Completed." << endl;
}

void Task::saveHistograms(TFile * outputFile)
{
  if (reportDebug()) cout << "Task::saveHistograms(TFile * outputFile) No ops  for " << getName() << endl;
  if (!outputFile && reportError())
    {
    cout << "Task::saveHistograms(TFile * outputFile)  Given file pointer is null." << endl;
    }
}

//////////////////////////////////////////////////////////////
// add histograms to an external list
//////////////////////////////////////////////////////////////
void Task::addHistogramsToExtList(TList *list, bool all)
{
  if (reportDebug()) cout << "Task::addHistogramsToExtList(TList *list) No ops" << endl;
  if (!list && reportError())
    {
    cout << "Task::addHistogramsToExtList(TList *list)  Given file pointer is null." << endl;
    }
  all = false; // silence warning.
}

void Task::calculateDerivedHistograms()
{
  if (reportDebug()) cout << "Task::calculateDerivedHistograms() No ops" << endl;
}

Event * Task::getEvent()
{
  return event;
}

//enum TaskStatus   { Unknown, TaskOk, TaskEof, TaskEod, TaskWarning, TaskError, TaskFatal};


Task::TaskStatus Task::taskStatus = Task::TaskOk;

Task::TaskStatus Task::getTaskStatus()
{
  return taskStatus;
}

void Task::setTaskStatus(Task::TaskStatus newStatus)
{
  taskStatus = newStatus;
}

void Task::postTaskOk()
{
  taskStatus = TaskOk;
}

void Task::postTaskEof()
{
  taskStatus = TaskEof;
}


void Task::postTaskEod()
{
  taskStatus = TaskEod;
}

void Task::postTaskWarning()
{
  taskStatus = TaskWarning;
}

void Task::postTaskError()
{
  taskStatus = TaskError;
}

void Task::postTaskFatal()
{
  taskStatus = TaskFatal;
}


bool Task::isTaskOk()
{
  return (taskStatus == TaskOk);
}


TString Task::getTaskStatusName()
{
  TString statusName;
  switch (taskStatus)
    {
      case Unknown:     statusName = "Unknown"; break;
      case TaskOk:      statusName = "TaskOk";  break;
      case TaskEof:     statusName = "TaskEof"; break;
      case TaskEod:     statusName = "TaskEod"; break;
      case TaskWarning: statusName = "TaskWarning"; break;
      case TaskError:   statusName = "TaskError";   break;
      case TaskFatal:   statusName = "TaskFatal";   break;
    }
  return statusName;
}

//enum TaskStatus   { Unknown, TaskOk, TaskEof, TaskEod, TaskWarning, TaskError, TaskFatal};
