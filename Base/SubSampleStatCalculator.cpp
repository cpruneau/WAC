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
#include "HistogramCollection.hpp"
#include "SubSampleStatCalculator.hpp"
#include "TParameter.h"

ClassImp(SubSampleStatCalculator);


SubSampleStatCalculator::SubSampleStatCalculator(const TString & _inputPath,
                                                 int _nFiles,
                                                 TString ** _fileNames,
                                                 const TString & _outputPath,
                                                 const TString & _outputFileName,
                                                 MessageLogger::LogLevel debugLevel)
: Task("SubSampleStatCalculator",nullptr,nullptr,debugLevel),
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////// "EventAccepted"


void SubSampleStatCalculator::execute()
{
  if (reportInfo("SubSampleStatCalculator",getTaskName(),"execute()"))
    {
    cout << "Starting subsample analysis of nFiles:" << nFiles << endl;
    cout << "    Input path: " << inputPath << endl;
    }
  TString inputFileName;
  TFile * inputFile;
  TFile * firstInputFile;
  sumWeights = 0.0;
  weight = 0.0;

  HistogramCollection * collectionAvg = new HistogramCollection("Sum",400);
  HistogramCollection * collection;
  TString parameterName = "EventAccepted";

  for (int iFile=0; iFile<nFiles; iFile++)
  {
  inputFileName = *fileNames[iFile];
  inputFile = openRootFile(inputPath, inputFileName, "READ");
  if (!isTaskOk()) return;
  weight = readParameter(inputFile,parameterName);
  if (!isTaskOk()) return;
  if (reportDebug ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Loading histograms" << endl;
  if (iFile==0)
    {
    firstInputFile = inputFile;
    collectionAvg->loadCollection(inputFile);
    if (reportDebug ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "First Load completed."  << endl;
    sumWeights = weight;
  }
  else
    {
    collection = new HistogramCollection(inputFileName,400);;
    collection->loadCollection(inputFile);
    collectionAvg->squareDifferenceCollection(*collection, sumWeights, weight, (iFile==(nFiles-1)) ? nFiles : -iFile);
    sumWeights += weight;
    delete collection;
    delete inputFile;
    if (reportDebug ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "sumWeights : " <<  sumWeights << endl;
    }
  }
  TFile * outputFile = openRootFile(outputPath, outputFileName, "RECREATE");
  if (!isTaskOk()) return;
  collectionAvg->saveHistograms(outputFile);
  delete collectionAvg;
  delete outputFile;
  delete firstInputFile;
  if (reportInfo ("SubSampleStatCalculator",getTaskName(),"execute()")) cout << "Subsample analysis completed." << endl;
}

