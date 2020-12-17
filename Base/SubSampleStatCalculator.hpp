// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_SubSampleStatCalculator
#define WAC_SubSampleStatCalculator
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
#include "TString.h"
#include "Task.hpp"
using namespace std;

class SubSampleStatCalculator : public Task
{
public:

  TString inputPath;
  int nFiles;
  TString ** fileNames;
  TString outputPath;
  TString outputFileName;
  double weight;
  double sumWeights;
  
  SubSampleStatCalculator(const TString & _inputPath,
                          int _nFiles,
                          TString ** _fileNames,
                          const TString & _outputPath,
                          const TString & _outputFileName,
                          MessageLogger::LogLevel debugLevel);
  virtual ~SubSampleStatCalculator();
  virtual void execute();

  ClassDef(SubSampleStatCalculator,0)
};

#endif /* WAC_SubSampleStatCalculator */
