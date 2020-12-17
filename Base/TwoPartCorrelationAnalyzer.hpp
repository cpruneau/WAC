// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_TwoPartCorrelationAnalyzer
#define WAC_TwoPartCorrelationAnalyzer
#include <TParameter.h>
#include "TFile.h"
#include "TList.h"
#include "Task.hpp"
#include "Event.hpp"
#include "EventFilter.hpp"
#include "EventHistos.hpp"
#include "ParticleFilter.hpp"
#include "ParticlePairFilter.hpp"
#include "ParticleHistos.hpp"
#include "ParticlePairHistos.hpp"
#include "ParticlePairDerivedHistos.hpp"
#include "ParticlePairCombinedHistos.hpp"
#include "ParticlePairAnalyzerConfiguration.hpp"

class TwoPartCorrelationAnalyzer : public Task
{
public:

  //////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////
  TwoPartCorrelationAnalyzer(const TString &  name,
                             ParticlePairAnalyzerConfiguration * configuration,
                             Event * event,
                             EventFilter * eventFilter,
                             ParticleFilter * particleFilter1,
                             ParticleFilter * particleFilter2,
                             LogLevel selectedLevel=Info);
  virtual ~TwoPartCorrelationAnalyzer();
  virtual void execute();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void loadBaseHistograms(TFile * inputFile);
  virtual void saveHistograms(TFile * outputFile);
  virtual void addHistogramsToExtList(TList *list, bool all=false);
  virtual void scaleHistograms(double factor);
  virtual void calculateDerivedHistograms();
  virtual void resetHistograms();

  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////

  EventFilter    * eventFilter;
  ParticleFilter * particleFilter1;
  ParticleFilter * particleFilter2;
  ParticleHistos             * particle1_Histos;
  ParticleHistos             * particle2_Histos;
  ParticlePairHistos         * pair11_Histos;
  ParticlePairHistos         * pair22_Histos;
  ParticlePairHistos         * pair12_Histos;
  ParticlePairDerivedHistos  * pair11_DerivedHistos;
  ParticlePairDerivedHistos  * pair22_DerivedHistos;
  ParticlePairDerivedHistos  * pair12_DerivedHistos;
  ParticlePairCombinedHistos * pair12_CIHistos;
  ParticlePairCombinedHistos * pair12_CDHistos;
  TString partName1;
  TString partName2;

   ClassDef(TwoPartCorrelationAnalyzer,0)
};


#endif /* WAC_TwoPartCorrelationAnalyzer */
