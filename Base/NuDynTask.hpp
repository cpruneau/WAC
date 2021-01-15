// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_NuDynTask
#define WAC_NuDynTask
#include <TParameter.h>
#include "TFile.h"
#include "TList.h"
#include "Task.hpp"
#include "Event.hpp"
#include "EventFilter.hpp"
#include "NuDynHistos.hpp"
#include "NuDynDerivedHistos.hpp"
#include "NuDynConfiguration.hpp"
#include "ParticleFilter.hpp"
#include "NuDynConfiguration.hpp"

class NuDynTask : public Task
{
public:

  //////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////
  NuDynTask(const TString &  name,
            NuDynConfiguration * configuration,
            Event * event,
            EventFilter * eventFilter,
            ParticleFilter * particleFilter1,
            ParticleFilter * particleFilter2,
            ParticleFilter * particleFilter3,
            ParticleFilter * particleFilter4,
            LogLevel selectedLevel );

  virtual ~NuDynTask();
  virtual void execute();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  //virtual void loadBaseHistograms(TFile * inputFile);
  virtual void saveHistograms(TFile * outputFile);
  //  virtual void addHistogramsToExtList(TList *list, bool all=false);
  virtual void scaleHistograms(double factor);
  virtual void calculateDerivedHistograms();
  virtual void resetHistograms();

  virtual void createIdentical();

  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////
  NuDynHistos        * nuDynHistos;
  NuDynDerivedHistos * nuDynDerivedHistos;
  EventFilter        * eventFilter;
  int nParticleFilters;
  ParticleFilter     ** particleFilters;
  TString            ** partNames;
  int identical[16];

  ClassDef(NuDynTask,0)
};


#endif /* WAC_NuDynTask */
