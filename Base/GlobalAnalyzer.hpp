// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_GlobalAnalyzer
#define WAC_GlobalAnalyzer
#include <TParameter.h>
#include "TFile.h"
#include "TList.h"
#include "Task.hpp"
#include "Event.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "GlobalHistos.hpp"
#include "GlobalAnalyzerConfiguration.hpp"

class GlobalAnalyzer : public Task
{
public:

  //GlobalAnalyzer();

  GlobalAnalyzer(const TString &  name,
                 GlobalAnalyzerConfiguration * configuration,
                 Event * event,
                 EventFilter * eventFilter,
                 int nParticleFilters,
                 ParticleFilter ** particleFilters,
                 LogLevel requiredLevel=Info);
  virtual ~GlobalAnalyzer();
  virtual void execute();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void saveHistograms(TFile * outputFile);
  virtual void scaleHistograms(double factor);
  virtual void resetHistograms();

  EventFilter     *  eventFilter;
  int nParticleFilters;
  ParticleFilter  ** particleFilters;
  GlobalHistos    *  globalHistos;
  TString         ** filterNames;
  double          *  n;
  double          *  e;
  double          *  q;
  double          *  b;

  ClassDef(GlobalAnalyzer,0)
};


#endif /* WAC_GlobalAnalyzer */
