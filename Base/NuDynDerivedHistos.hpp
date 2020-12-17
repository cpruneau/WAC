// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_NuDynDerivedHistos
#define WAC_NuDynDerivedHistos
#include "Histograms.hpp"
#include "NuDynHistos.hpp"


class NuDynDerivedHistos : public Histograms
{
public:

  NuDynDerivedHistos(const TString & collectionName,
                     NuDynConfiguration * analysisConfiguration,
                      LogLevel  debugLevel);
  virtual ~NuDynDerivedHistos() {} 
  void createHistograms();
  void loadHistograms(TFile * inputFile);
  void calculateDerivedHistograms(NuDynHistos * source);

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  // n1 number of particles satisfying filter 1
  // n2 number of particles satisfying filter 2
  // h_Fj  = cumulants of order "j"
  // h_Rj  = normalized cumulants of order "j"

  TH1 ** h_F2;
  TH1 ** h_F3;
  TH1 ** h_F4;
  TH1 ** h_R2;
  TH1 ** h_R3;
  TH1 ** h_R4;
  TH1 ** h_nudyn;

  // vs Mult
  TH1 ** h_F2_vsMult;
  TH1 ** h_F3_vsMult;
  TH1 ** h_F4_vsMult;
  TH1 ** h_R2_vsMult;
  TH1 ** h_R3_vsMult;
  TH1 ** h_R4_vsMult;
  TH1 ** h_nudyn_vsMult;

  ClassDef(NuDynDerivedHistos,0)
};

#endif /* WAC_NuDynDerivedHistos  */



