// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_NuDynHistos
#define WAC_NuDynHistos
#include "Histograms.hpp"

class NuDynHistos : public Histograms
{
public:

  NuDynHistos(const TString & collectionName,
                      AnalysisConfiguration * analysisConfiguration,
                      LogLevel  debugLevel);
  NuDynHistos(TFile * inputFile,
                      const TString & collectionName,
                      AnalysisConfiguration * analysisConfiguration,
                      LogLevel  debugLevel);
  virtual ~NuDynHistos();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void fill(double mult, double cent, double * n,  double weight);


  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  // n1 number of particles satisfying filter 1
  // n2 number of particles satisfying filter 2
  // n3 number of particles satisfying filter 3
  // n4 number of particles satisfying filter 4
  // h_f1_1 = <n1>
  // h_f1_2 = <n2>
  // h_f2_11 = <n1(n1-1)>
  // h_f2_22 = <n2(n2-1)>
  // h_f2_12 = <n2(n2-1)>
  // h_f3_111 = <n1(n1-1)(n1-2)>
  // h_f3_112 = <n1(n1-1)n2>
  // h_f3_113 = <n1(n1-1)n3>
  // h_f3_122 = <n1n2(n2-1)>
  // h_f3_123 = <n1n2n3>
  // etc

  // h_r2_11 = <n1(n1-1)>/<n1><n1>
  // h_r2_22 = <n2(n2-1)>/<n2><n2>
  // h_r2_12 = <n1(n2-1)>/<n1><n2>
  // h_nudyn_12 = h_r2_11 + h_r2_22 -2*h_r2_12
  // etc

  // Min bias all included
  TH1 * h_events;
  TProfile ** h_f1;
  TProfile ** h_f2;
  TProfile ** h_f3;
  TProfile ** h_f4;


  // vs Mult measured in fiducial
   TH1 * h_events_vsMult;
  TProfile ** h_f1_vsMult;
  TProfile ** h_f2_vsMult;
  TProfile ** h_f3_vsMult;
  TProfile ** h_f4_vsMult;

  // vs Centrality
  TH1 * h_events_vsCent;
  TProfile ** h_f1_vsCent;
  TProfile ** h_f2_vsCent;
  TProfile ** h_f3_vsCent;
  TProfile ** h_f4_vsCent;

  ClassDef(NuDynHistos,0)
};

#endif /* WAC_NuDynHistos  */



