// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/
#ifndef WAC_EventHistos
#define WAC_EventHistos
#include "Histograms.hpp"
#include "EventAnalyzerConfiguration.hpp"

class EventHistos : public Histograms
{
public:
  
  
  EventHistos(const TString & name,
              EventAnalyzerConfiguration * configuration,
              LogLevel  debugLevel);
  EventHistos(TFile * inputFile,
              const TString & name,
              EventAnalyzerConfiguration * configuration,
              LogLevel  debugLevel);
  virtual ~EventHistos();
  virtual void initialize();
  void loadHistograms(TFile * inputFile);
  void fill(double mult, double weight);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data Members -- Computed
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 * h_nPartTot;

    ClassDef(EventHistos,0)
};

#endif /* WAC_EventHistos */
