// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class SubSampleStatCalculator
 \ingroup WAC

 Class defining SubSampleStatCalculator
 */

#include "SubSampleStatCalculator.hpp"
#include "TParameter.h"

ClassImp(SubSampleStatCalculator);


SubSampleStatCalculator::SubSampleStatCalculator(TaskConfiguration * configuration,
                                                 int _nStore,
                                                 Histograms ** _store,
                                                 const TString & _inputPath,
                                                 int _nFiles,
                                                 TString ** _fileNames,
                                                 const TString & _outputPath,
                                                 const TString & _outputFileName,
                                                 MessageLogger::LogLevel debugLevel)
: Task("SubSampleStatCalculator",configuration,nullptr,debugLevel),
store(_store),
reader(_reader),
inputPath(_inputPath),
nFiles(_nFiles),
fileNames(_fileNames),
outputPath(_outputPath),
outputFileName(_outputFileName),
weight(0),
sumWeights(0)
{
  if (reportDebug("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "No ops." << endl;
}

SubSampleStatCalculator::~SubSampleStatCalculator()
{
  if (reportDebug("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "No ops." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Execute task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SubSampleStatCalculator::execute()
{
  if (reportInfo("SubSampleStatCalculator",getTaskName(),"execute()"))
    {
    cout << "Starting subsample analysis of nFiles:" << nFiles << endl;
    cout << "    Input path: " << inputPath << endl;
    cout << "    Sit down, relax, this may take a while..." << endl;
    }

  // store->createHistograms();

  TString inputFileName;
  TFile * inputFile;
  double sumWeights = 0.0;
  double weight = 0.0;
  for (int iFile=0; iFile<nFiles; iFile++)
  {
  inputFileName = inputPath;
  inputFileName += "/";
  inputFileName += *fileNames[iFile];
  inputFileName += ".root";
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Opening file: " << inputFileName << endl;
  inputFile = new TFile(inputFileName,"READ");
  if (!inputFile)
    {
    if (reportError("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "file: " << inputFileName << " not found." << endl;
    return;
    }
  else
    {
    if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "file successfully opened." << endl;
    }
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Reading EventAccepted." << endl;
  TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get("EventAccepted");
  if (!par)
    {
    if (reportError("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "EventAccepted not found." << endl;
    return;
    }
  weight = par->GetVal();
  delete par;
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "EventAccepted : " << weight << endl;
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Loading histograms" << endl;

  HistogramCollection * collection = new HistogramCollection("temp",400);
  collection->loadHistosInList(inputFile, store);
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Load completed."  << endl;
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Compute square differences..."  << endl;
  store->squareDifferenceCollection(*collection, sumWeights, weight, (iFile==(nFiles-1)) ? nFiles : 0);
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Computation completed. Delete collection."  << endl;
  delete collection;
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Computation completed. Delete input file object."  << endl;
  delete inputFile;
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Next clear."  << endl;
  //reader->clear();
  sumWeights += weight;
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "sumWeights : " <<  sumWeights << endl;
  }

  TString fileName = outputPath;
  fileName += "/";
  fileName += outputFileName;
  fileName += ".root";
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Opening output file: " << fileName << endl;
  TFile * outputFile = new TFile(fileName,"RECREATE");
  if (!outputFile)
    {
    if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Could not open output file: " << fileName << endl;
    return;
    }
  store->saveHistograms(outputFile);
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Subsample analysis completed/saved to  file: " << fileName << endl;
}

