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
#include "TParameter.h"
#include "Task.hpp"

ClassImp(Task);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Task::Task(const TString & name,
           TaskConfiguration * configuration,
           Event * selectedEvent,
           LogLevel selectedLevel)
:
MessageLogger(selectedLevel),
taskName             ( name ),
taskConfiguration    ( configuration ),
taskRandomGenerator  ( gRandom),
event                ( selectedEvent ),
nEventProcessed      ( 0 ),
nEventAccepted      ( 0 ),
subSampleIndex       ( 0 )
{
// no ops
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Task::~Task()
{
  // no ops
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Task::initialize()
{
  if (isTaskOk() && taskConfiguration->loadHistograms)   loadHistograms();
  if (isTaskOk() && taskConfiguration->createHistograms) createHistograms();
  nEventProcessed = 0;
  nEventAccepted = 0;
  if (reportEnd("Task",getTaskName(),"initialize()"))
    ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Execute task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Task::execute()
{
  if (reportNoOps("Task",getTaskName(),"execute()"))
    ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Finalize task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Task::finalize()
{
  if (reportStart("Task",getTaskName(),"finalize()"))
    ;
  if (nEventAccepted>0 && nEventProcessed>0)
    {
    if (reportInfo("Task",getTaskName(),"finalize()"))
      {
      cout << " nEventProcessed:" << nEventProcessed << endl;
      cout << " nEventAccepted:" << nEventAccepted << endl;
      }
    if (isTaskOk() && taskConfiguration->scaleHistograms)            scaleHistograms();
    if (isTaskOk() && taskConfiguration->calculateDerivedHistograms) calculateDerivedHistograms();
    if (isTaskOk() && taskConfiguration->saveHistograms)             saveHistograms();
    }
  else
    {
    if (reportWarning("Task",getTaskName(),"finalize()"))
      {
      cout << " nEventProcessed:" << nEventProcessed << endl;
      cout << " nEventAccepted:" << nEventAccepted << endl;
      }
    }
  if (reportEnd("Task",getTaskName(),"finalize()"))
    ;
}

void Task::savePartialResults()
{
  if (reportStart("Task",getTaskName(),"savePartialResults()"))
    ;
  if (nEventAccepted>0 && nEventProcessed>0)
    {
    if (reportInfo("Task",getTaskName(),"savePartialResults()"))
      {
      cout << " nEventProcessed:" << nEventProcessed << endl;
      cout << " nEventAccepted:" << nEventAccepted << endl;
      }
    if (isTaskOk() && taskConfiguration->scaleHistograms)            scaleHistograms();
    if (isTaskOk() && taskConfiguration->calculateDerivedHistograms) calculateDerivedHistograms();
    if (isTaskOk() && taskConfiguration->saveHistograms)             saveHistograms();
    }
  else
    {
    if (reportWarning("Task",getTaskName(),"savePartialResults()"))
      {
      cout << " nEventProcessed:" << nEventProcessed << endl;
      cout << " nEventAccepted:" << nEventAccepted << endl;
      }
    }
  if (reportEnd("Task",getTaskName(),"savePartialResults()"))
    ;
}


void Task::reset()
{
  if (reportStart("Task",getTaskName(),"reset()"))
    ;
  nEventProcessed = 0;
  nEventAccepted  = 0;
  if (isTaskOk())   resetHistograms();
  if (reportEnd("Task",getTaskName(),"reset()"))
    ;
}

////////////////////////////////////////////
// Reset this task
////////////////////////////////////////////
void Task::clear()
{
  if (reportStart("Task",getTaskName(),"clear()"))
    ;
  reset();
  if (reportEnd("Task",getTaskName(),"clear()"))
    ;
}

TaskConfiguration * Task::getTaskConfiguration()
{
  return taskConfiguration;
}

void Task::setTaskConfiguration(TaskConfiguration * config)
{
  if (!config)
    {
    taskConfiguration = config;
    if (reportInfo("Task",getTaskName(),"setTaskConfiguration(...)"))
      {
      cout << "Task configuration was changed." << endl;
      }
    }
  else
    {
    if (reportError("Task",getTaskName(),"setTaskConfiguration(...)"))
      {
      cout << "Attempting to set task configuration to a null pointer."  << endl;
      }
    postTaskError();
    }
}

////////////////////////////////////////////
// Print task configuration
////////////////////////////////////////////
void Task::printTaskConfiguration(ostream & output)
{
  output << "Task Name : " << taskName <<  endl;
  taskConfiguration->printTaskConfiguration(output);
}

void Task::setRandomGenerator(TRandom * randomGenerator)
{
  if (!randomGenerator)
    {
    taskRandomGenerator = randomGenerator;
    if (reportInfo("Task",getTaskName(),"setRandomGenerator(TRandom * randomGenerator)"))
      {
      cout << "Default random generator was changed." << endl;
      }
    }
  else
    {
    postTaskWarning();
    if (reportError("Task",getTaskName(),"setTaskConfiguration(...)"))
      {
      cout << "Null pointer. Random generator will not be set." << endl;
      }
    }
}

void Task::createHistograms()
{
  if (reportNoOps("Task",getTaskName(),"createHistograms()"))
    ;
}

void Task::loadHistograms()
{
  if (reportStart("Task",getTaskName(),"loadHistograms()"))
    ;
  TFile * inputFile;
  TString inputFileName = taskConfiguration->getInputRootFileName();
  if (reportInfo("Task",getTaskName(),"loadHistograms()")) cout << "Opening input file: " << inputFileName << endl;
  inputFile = new TFile(inputFileName,"OLD");
  if (!inputFile)
     {
     if (reportError("Task",getTaskName(),"loadHistograms()")) cout << "Could not open input file:" << inputFileName << endl;
     postTaskError();
     return;
    }
  loadHistograms(inputFile);
  if (reportEnd("Task",getTaskName(),"loadHistograms()"))
    ;
}

void Task::loadHistograms(TFile * inputFile)
{
  inputFile = nullptr;
  if (reportNoOps("Task",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
}

void Task::resetHistograms()
{
  if (reportNoOps("Task",getTaskName(),"resetHistograms()"))
    ;
}

void Task::clearHistograms()
{
  if (reportNoOps("Task",getTaskName(),"clearHistograms()"))
    ;
}

void Task::scaleHistograms()
{
  if (reportStart("Task",getTaskName(),"scaleHistograms()"))
    ;
  if (reportInfo("Task",getTaskName(),"scaleHistograms()"))
    {
    cout << "Accepted number of events: " <<  nEventAccepted << endl;
    }
  if (nEventAccepted>0)
    {
    double scalingFactor = 1.0/double(nEventAccepted);
    scaleHistograms(scalingFactor);
    }
  else
    {
    if (reportWarning("Task",getTaskName(),"scaleHistograms()"))
      {
      cout << "Accepted number of events is null: " <<  nEventAccepted << endl;
      }
    }

  if (reportStart("Task",getTaskName(),"scaleHistograms()"))
    ;
}

void Task::scaleHistograms(double factor)
{
  factor=0;
  if (reportNoOps("Task",getTaskName(),"scaleHistograms(double factor)"))
    ;
}

void Task::saveHistogramsAsText()
{
  if (reportNoOps("Task",getTaskName(),"saveHistogramsAsText()"))
    ;
}

void Task::saveHistograms()
{
  if (reportStart("Task",getTaskName(),"saveHistograms()"))
    ;
  TFile * outputFile;
  TString outputFileName = taskConfiguration->outputPath;
  outputFileName += taskConfiguration->rootOuputFileName;
  outputFileName += getTaskName();
  if (taskConfiguration->subsampleAnalysis)
    {
    outputFileName += "_";
    outputFileName += subSampleIndex++;
    }
  outputFileName += ".root";

  if (taskConfiguration->forceHistogramsRewrite)
    {
    if (reportInfo("Task",getTaskName(),"saveHistograms()")) cout << "Opening output (RECREATE) file  " << outputFileName << endl;
    outputFile = new TFile(outputFileName,"RECREATE"); // obliterate past work...
    if (!outputFile)
      {
      if (reportError("Task",getTaskName(),"saveHistograms()")) cout << "Could not open (RECREATE) file  " << outputFileName << endl;
      postTaskWarning();
      return;
      }
    }
  else
    {
    if (reportInfo("Task",getTaskName(),"saveHistograms()")) cout << "Opening  output (NEW) file  " << outputFileName << endl;
    outputFile = new TFile(outputFileName,"NEW"); // protect past work...
    if (!outputFile)
      {
      if (reportError("Task",getTaskName(),"saveHistograms()")) cout << "Could not open (NEW) file  " << outputFileName << endl;
      postTaskWarning();
      return;
      }
    }
  saveNEventProcessed(outputFile);
  saveNEventAccepted(outputFile);
  saveHistograms(outputFile);
  outputFile->Close();
  if (reportEnd("Task",getTaskName(),"saveHistograms()"))
    ;
}

void Task::saveHistograms(TFile * outputFile)
{
  outputFile = nullptr;
  if (reportNoOps("Task",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
}

void Task::saveNEventProcessed(TFile * outputFile)
{
  outputFile->cd();
  TParameter<Long64_t>("EventProcessed",nEventProcessed,'+').Write();
}

void Task::saveNEventAccepted(TFile * outputFile)
{
  outputFile->cd();
  TParameter<Long64_t>("EventAccepted",nEventAccepted,'+').Write();
}

long Task::loadNEventProcessed(TFile * inputFile)
{
  TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get("nEventProcessed");
  nEventProcessed = par->GetVal();
  delete par;
  return nEventProcessed;
}
long Task::loadNEventAccepted(TFile * inputFile)
{
  TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get("nEventAccepted");
  nEventAccepted = par->GetVal();
  delete par;
  return nEventAccepted;
}

//////////////////////////////////////////////////////////////
// add histograms to an external list
//////////////////////////////////////////////////////////////
void Task::addHistogramsToExtList(TList *list, bool all)
{
  if (reportInfo("Task",getTaskName(),"addHistogramsToExtList(...)"))
    ;
  if (!list && reportError())
    {
    postTaskWarning();
    if (reportError("Task",getTaskName(),"addHistogramsToExtList(...)")) cout << "Given file pointer is null." << endl;
    }
  all = false; // silence warning.
}

void Task::calculateDerivedHistograms()
{
  if (reportNoOps("Task",getTaskName(),"calculateDerivedHistograms()"))
    ;
}

Event * Task::getEvent()
{
  return event;
}

//enum TaskStatus   { Unknown, TaskOk, TaskEof, TaskEod, TaskWarning, TaskError, TaskFatal};


Task::TaskStatus Task::taskStatus = Task::TaskOk;




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


double Task::readParameter(TFile * inputFile, const TString & parameterName)
{
  if (reportDebug ("Task",getTaskName(),"readParameter(...)")) cout << "Reading parameter: " << parameterName << endl;
  TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get(parameterName);
  if (!par)
    {
    if (reportError("Task",getTaskName(),"execute()")) cout << "Parameter not found." << endl;
    postTaskError();
    return -1;
    }
  double value = par->GetVal();
  delete par;
  if (reportDebug ("Task",getTaskName(),"readParameter(...)")) cout << "Parameter value : " << value << endl;
  return value;
}

TFile *  Task::openRootFile(const TString & inputPath, const TString & fileName, const TString & ioOption)
{
  TString inputFileName = inputPath;
  inputFileName += "/";
  inputFileName += fileName;
  inputFileName += ".root";
  if (reportDebug ("Task",getTaskName(),"openRootFile(...)"))
    cout << "Opening file: " << inputFileName << " with option: " << ioOption << endl;
  TFile * inputFile = new TFile(inputFileName,ioOption);
  if (!inputFile)
    {
    if (reportError("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "file: " << inputFileName << " not found." << endl;
    postTaskError();
    return nullptr;
    }
  else
    {
    if (reportDebug ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "file opened successfully." << endl;
    }
  return inputFile;
}


//enum TaskStatus   { Unknown, TaskOk, TaskEof, TaskEod, TaskWarning, TaskError, TaskFatal};
