// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_HadronGasParticlesHistograms
#define WAC_HadronGasParticlesHistograms
#include "Histograms.hpp"
#include "HadronGasConfiguration.hpp"
#include "HadronGas.hpp"

// ===================================================
// Particle abundance at a specific temperature or beam
// energy or etc...
class HadronGasParticleHistograms : public Histograms
{
public:

  HadronGasParticleHistograms(const TString & collectionName,
                               HadronGasConfiguration * hadronGasConfiguration,
                               LogLevel  debugLevel);
  virtual ~HadronGasParticleHistograms();
  void createHistograms();
  void loadHistograms(TFile * inputFile);
  void fill(HadronGas & hadronGas);
  void calculateDerivedHistograms();

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  TH1 * h_allYields;
  TH1 * h_allYieldsToPion;
  TH1 * h_allToAntiRatio;
  TH1 * h_stableThermalYields;
  TH1 * h_stableThermalYieldsToPion;
  TH1 * h_stableThermalToAntiRatio;
  TH1 * h_stableDecayYields;
  TH1 * h_stableDecayYieldsToPion;
  TH1 * h_stableDecayToAntiRatio;
  TH1 * h_stableDecayToThermalRatio;

  TH2 * h_stableDecayPairsYields;
  TH2 * h_stableDecayCorrelatedPairsYields;
  TH2 * h_stableDecayCorrelatedPairsNorm;


    ClassDef(HadronGasParticleHistograms,0)

};




#endif /* WAC_HadronGasHistograms  */



