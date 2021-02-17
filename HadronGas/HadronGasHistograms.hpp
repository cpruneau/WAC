// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_HadronGasHistograms
#define WAC_HadronGasHistograms
#include "Histograms.hpp"
#include "HadronGasConfiguration.hpp"
#include "HadronGas.hpp"

class HadronGasHistograms : public Histograms
{
public:

  HadronGasHistograms(const TString & collectionName,
                      HadronGasConfiguration * analysisConfiguration,
                      LogLevel  debugLevel);
  virtual ~HadronGasHistograms();
  void createHistograms();
  void loadHistograms(TFile * inputFile);
  void fill(HadronGas & hadronGas);

  TH1 * sysEnergyVsT;
  TH1 * sysEntropyVsT;
  TH1 * sysPressureVsT;
  TH1 * sysNetBVsT;
  TH1 * sysNetQVsT;
  TH2 * sysNetBVsTMuB;
  TH2 * sysNetQVsTMuB;

  ClassDef(HadronGasHistograms,0)

};

#endif /* WAC_HadronGasHistograms  */



