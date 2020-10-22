// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_ParticleAnalyzer
#define WAC_ParticleAnalyzer
#include <TParameter.h>
#include "TFile.h"
#include "TList.h"
#include "Task.hpp"
#include "Event.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "ParticleHistos.hpp"

class ParticleAnalyzer : public Task
{
public:

  //ParticleAnalyzer();

  ParticleAnalyzer(const TString &  name,
                   TaskConfiguration * configuration,
                   Event * event,
                   EventFilter * eventFilter,
                   int nParticleFilters,
                   ParticleFilter ** particleFilters);
  virtual ~ParticleAnalyzer();
  virtual void execute();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void saveHistograms(TFile * outputFile);
  virtual void scaleHistograms(double factor);

  int nParticleFilters;
  EventFilter     *  eventFilter;
  ParticleFilter  ** particleFilters;
  ParticleHistos  ** particleHistos;
  TString         ** partNames;
  double          *  nAccepted;

  ClassDef(ParticleAnalyzer,0)
};


#endif /* WAC_ParticleAnalyzer */
