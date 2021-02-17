//
//  HadronGasHistograms.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "HadronGasParticleHistograms.hpp"
ClassImp(HadronGasParticleHistograms);

HadronGasParticleHistograms::HadronGasParticleHistograms(const TString & name,
                                         HadronGasConfiguration * configuration,
                                         LogLevel  debugLevel)
:
Histograms(name,configuration,100,debugLevel)
{
  // no ops
}

HadronGasParticleHistograms::~HadronGasParticleHistograms()
{
  //deleteHistograms();
}

// for now use the same boundaries for eta and y histogram
void HadronGasParticleHistograms::createHistograms()
{
  HadronGasConfiguration & ac = *(HadronGasConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();
  h_allYields                  = createHistogram(bn+TString("allYields"),          ac.nTypes,        0.0, double(ac.nTypes+1),  "Species","Yields", false,false,false,false);
  h_allYieldsToPion            = createHistogram(bn+TString("allYieldsToPi"),      ac.nTypes,        0.0, double(ac.nTypes+1),  "Species","YieldsToPi", false,false,false,false);
  h_allToAntiRatio             = createHistogram(bn+TString("allYieldsToAnti"),    ac.nTypes,        0.0, double(ac.nTypes+1),  "Species","YieldsToAnti", false,false,false,false);
  h_stableThermalYields        = createHistogram(bn+TString("stableThermalYields"),       ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","ThermalYields", false,false,false,false);
  h_stableThermalYieldsToPion  = createHistogram(bn+TString("stableThermalYieldsToPi"),   ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","ThermalYieldsToPi", false,false,false,false);
  h_stableThermalToAntiRatio   = createHistogram(bn+TString("stableThermalYieldsToAnti"), ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","ThermalYieldsToAnti", false,false,false,false);
  h_stableDecayYields          = createHistogram(bn+TString("stableDecayYields"),         ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","DecayYields", false,false,false,false);
  h_stableDecayYieldsToPion    = createHistogram(bn+TString("stableDecayYieldsToPi"),     ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","DecayYieldsToPi", false,false,false,false);
  h_stableDecayToAntiRatio     = createHistogram(bn+TString("stableDecayYieldsToAnti"),   ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","DecayYieldsToAnti", false,false,false,false);
  h_stableDecayToThermalRatio  = createHistogram(bn+TString("stableDecayToThermalRatio"), ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","wDecayToThermal", false,false,false,false);

  h_stableDecayPairsYields           = createHistogram(bn+TString("stableDecayPairYields"),
                                                       ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                                       ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                                       "Species","Species","Pairs", false,false,false,false);
  h_stableDecayCorrelatedPairsYields = createHistogram(bn+TString("stableDecayCorrelatedPairYields"),
                                                       ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                                       ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                                       "Species","Species","Correlated Pairs", false,false,false,false);
  h_stableDecayCorrelatedPairsNorm  = createHistogram(bn+TString("stableDecayCorrelatedPairNorm"),
                                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                                      "Species","Species","Correlated Pairs", false,false,false,false);
}

void HadronGasParticleHistograms::calculateDerivedHistograms()
{
  h_stableDecayToThermalRatio->Divide(h_stableDecayYields,h_stableThermalYields);

  // normalize per trigger...
  int n = h_stableDecayCorrelatedPairsYields->GetNbinsX();
  for (int iX = 1; iX<=n;iX++)
  {
  double sum = 0.0;
  for (int iY = 1; iY<=n;iY++)
    {
    double v = h_stableDecayCorrelatedPairsYields->GetBinContent(iX,iY);
    sum += v;
    }
  if (sum>0)
    {
    for (int iY = 1; iY<n;iY++)
      {
      double v = h_stableDecayCorrelatedPairsYields->GetBinContent(iX,iY);
      h_stableDecayCorrelatedPairsNorm->SetBinContent(iX,iY,v/sum);
      h_stableDecayCorrelatedPairsNorm->SetBinError(iX,iY,0.0);
      }
    }
  }
}


//________________________________________________________________________
void HadronGasParticleHistograms::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportFatal()) cout << "Attempting to load HadronGasHistograms from an invalid file pointer" << endl;
    return;
    }
  HadronGasConfiguration & ac = *(HadronGasConfiguration*) getConfiguration();
  TString bn  = getHistoBaseName();
  h_allYields = loadH1(inputFile,bn+TString("allYields"));
  if (!h_allYields)
    {
    if (reportError()) cout << "Could not load histogram: " << bn+TString("allYields") << endl;
    return;
    }
  h_allYieldsToPion    = loadH1(inputFile,bn+TString("allYieldsToPi"));
  h_allToAntiRatio     = loadH1(inputFile,bn+TString("allYieldsToAnti"));
  h_stableThermalYields       = loadH1(inputFile,bn+TString("stableThermalYields"));
  h_stableThermalYieldsToPion = loadH1(inputFile,bn+TString("stableThermalYieldsToPi"));
  h_stableThermalToAntiRatio = loadH1(inputFile,bn+TString("stableThermalYieldsToAnti"));
  h_stableDecayYields        = loadH1(inputFile,bn+TString("stableDecayYields"));
  h_stableDecayYieldsToPion  = loadH1(inputFile,bn+TString("stableDecayYieldsToPi"));
  h_stableDecayToAntiRatio   = loadH1(inputFile,bn+TString("stableDecayYieldsToAnti"));
  h_stableDecayPairsYields             = loadH2(inputFile,bn+TString("stableDecayPairYields"));
  h_stableDecayCorrelatedPairsYields   = loadH2(inputFile,bn+TString("stableDecayCorrelatedPairYields"));
}


void HadronGasParticleHistograms::fill(HadronGas & hadronGas)
{
  int nTypes       = hadronGas.getNHadrons();
  int nStableTypes = hadronGas.getNStableHadrons();
  double thermalYield1, thermalYield2;
  double decayYield1,   decayYield2;
  for (int iType=0; iType<nTypes; iType++)
  {
  thermalYield1 = hadronGas.getParticle(iType)->getParticleYield();
  h_allYields->SetBinContent(iType+1, thermalYield1); h_allYields->SetBinError(iType+1, 0.0);
  }

  for (int iType1=0; iType1<nStableTypes; iType1++)
  {
  Hadron * part1 = hadronGas.getStableParticle(iType1);
  thermalYield1 = part1->getParticleYield();
  decayYield1   = part1->getStableParticleYield();
  h_stableThermalYields->SetBinContent(iType1+1, thermalYield1); h_stableThermalYields->SetBinError(iType1, 0.0);
  h_stableDecayYields  ->SetBinContent(iType1+1, decayYield1);   h_stableDecayYields  ->SetBinError(iType1, 0.0);

  for (int iType2=0; iType2<nStableTypes; iType2++)
    {
    Hadron * part2 = hadronGas.getStableParticle(iType2);
    thermalYield2 = part2->getParticleYield();
    decayYield2   = part2->getStableParticleYield();
    h_stableDecayPairsYields->SetBinContent(iType1+1, iType2+1, decayYield1*decayYield2);
    h_stableDecayPairsYields->SetBinError(  iType1+1, iType2+1, 0.0);
    h_stableDecayCorrelatedPairsYields->SetBinContent(iType1+1, iType2+1, hadronGas.stableParticlePairCorrelatedYield[iType1][iType2]);
    h_stableDecayCorrelatedPairsYields->SetBinError(  iType1+1, iType2+1, 0.0);
    }
  }
}

