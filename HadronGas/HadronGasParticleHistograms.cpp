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
  h_allYieldsVsMass            = createProfile(bn+TString("allYieldsVsMass"),      ac.nMass,         0.0, ac.maxMass,           "Mass (GeV)","Yields", false,false,false);
  h_allYieldsToPion            = createHistogram(bn+TString("allYieldsToPi"),      ac.nTypes,        0.0, double(ac.nTypes+1),  "Species","YieldsToPi", false,false,false,false);
  h_allToAntiRatio             = createHistogram(bn+TString("allYieldsToAnti"),    ac.nTypes,        0.0, double(ac.nTypes+1),  "Species","YieldsToAnti", false,false,false,false);
  h_stableThermalYields        = createHistogram(bn+TString("stableThermalYields"),       ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","ThermalYields", false,false,false,false);
  h_stableThermalVsMass        = createProfile(bn+TString(  "stableThermalYieldsVsMass"),
                                               ac.nMass,0.0, ac.maxMass,"Mass (GeV)","Yields", false,false,false);
  h_stableThermalYieldsToPion  = createHistogram(bn+TString("stableThermalYieldsToPi"),   ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","ThermalYieldsToPi", false,false,false,false);
  h_stableThermalToAntiRatio   = createHistogram(bn+TString("stableThermalYieldsToAnti"), ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","ThermalYieldsToAnti", false,false,false,false);
  h_stableDecayYields          = createHistogram(bn+TString("stableDecayYields"),         ac.nStableTypes,  0.0, double(ac.nStableTypes+1), "Species","DecayYields", false,false,false,false);
  h_stableDecayVsMass          = createProfile(bn+TString(  "stableDecayYieldsVsMass"),   ac.nMass,         0.0, ac.maxMass,                "Mass (GeV)","Yields", false,false,false);
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

  h_stableChargeBF  = createHistogram(bn+TString("stableDecayChargeBF"),
                                      18,  0.0, 18.0,
                                      "Pairs","BF", false,false,false,false);

  h_rho2            = createHistogram(bn+TString("rho2"),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      "Species","Species","#rho_{2}", false,false,false,false);
  h_rho1rho1        = createHistogram(bn+TString("rho1rho1"),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      "Species","Species","#rho_{1}#rho_{1}", false,false,false,false);
  h_rho1thrho1th    = createHistogram(bn+TString("rho1thrho1th"),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      "Species","Species","#rho_{1}^{Th}#rho_{1}^{Th}", false,false,false,false);
  h_C2              = createHistogram(bn+TString("C2"),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      "Species","Species","C_{2}", false,false,false,false);
  h_R2              = createHistogram(bn+TString("R2"),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      ac.nStableTypes,  0.0, double(ac.nStableTypes),
                                      "Species","Species","R_{2}", false,false,false,false);

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

// iX:1    Label:#pi^{0}
// iX:2    Label:#pi^{+}
// iX:3    Label:#pi^{-}
// iX:4    Label:K^{+}
// iX:5    Label:K^{-}
// iX:6    Label:K^{0}
// iX:7    Label:#bar K^{0}
// iX:8    Label:p
// iX:9    Label:#bar p
// iX:10    Label:n
// iX:11    Label:#bar n
// iX:12    Label:#Lambda
// iX:13    Label:#bar#Lambda
// iX:14    Label:#Sigma^{+}
// iX:15    Label:#bar#Sigma^{+}
// iX:16    Label:#Sigma^{-}
// iX:17    Label:#bar#Sigma^{-}
// iX:18    Label:#Xi^{0}
// iX:19    Label:#bar#Xi^{0}
// iX:20    Label:#Xi^{-}
// iX:21    Label:#bar#Xi^{-}
// iX:22    Label:#Omega
// iX:23    Label:#bar#Omega

  double zero = 0.0;
  for (int k=1;k<=18;k++)
  {
  h_stableChargeBF->SetBinError(k,zero);
  }
  // pi- | pi+
  h_stableChargeBF->SetBinContent(1,h_stableDecayCorrelatedPairsNorm->GetBinContent(2,3) - h_stableDecayCorrelatedPairsNorm->GetBinContent(2,2));
  // K- | pi+
  h_stableChargeBF->SetBinContent(2,h_stableDecayCorrelatedPairsNorm->GetBinContent(2,5) - h_stableDecayCorrelatedPairsNorm->GetBinContent(2,4));
  // p- | pi+
  h_stableChargeBF->SetBinContent(3,h_stableDecayCorrelatedPairsNorm->GetBinContent(2,9) - h_stableDecayCorrelatedPairsNorm->GetBinContent(2,8));
  // pi+ | pi-
  h_stableChargeBF->SetBinContent(4,h_stableDecayCorrelatedPairsNorm->GetBinContent(3,2) - h_stableDecayCorrelatedPairsNorm->GetBinContent(3,3));
  // K+ | pi-
  h_stableChargeBF->SetBinContent(5,h_stableDecayCorrelatedPairsNorm->GetBinContent(3,4) - h_stableDecayCorrelatedPairsNorm->GetBinContent(3,5));
  // p+ | pi-
  h_stableChargeBF->SetBinContent(6,h_stableDecayCorrelatedPairsNorm->GetBinContent(3,8) - h_stableDecayCorrelatedPairsNorm->GetBinContent(3,9));
  // pi- | K+
  h_stableChargeBF->SetBinContent(7,h_stableDecayCorrelatedPairsNorm->GetBinContent(4,3) - h_stableDecayCorrelatedPairsNorm->GetBinContent(4,2));
  // K- | K+
  h_stableChargeBF->SetBinContent(8,h_stableDecayCorrelatedPairsNorm->GetBinContent(4,4) - h_stableDecayCorrelatedPairsNorm->GetBinContent(4,4));
  // p- | K+
  h_stableChargeBF->SetBinContent(9,h_stableDecayCorrelatedPairsNorm->GetBinContent(4,8) - h_stableDecayCorrelatedPairsNorm->GetBinContent(4,8));
  // pi+ | K-
  h_stableChargeBF->SetBinContent(10,h_stableDecayCorrelatedPairsNorm->GetBinContent(5,2) - h_stableDecayCorrelatedPairsNorm->GetBinContent(5,3));
  // K+ | K-
  h_stableChargeBF->SetBinContent(11,h_stableDecayCorrelatedPairsNorm->GetBinContent(5,4) - h_stableDecayCorrelatedPairsNorm->GetBinContent(5,5));
  // p+ | K-
  h_stableChargeBF->SetBinContent(12,h_stableDecayCorrelatedPairsNorm->GetBinContent(5,8) - h_stableDecayCorrelatedPairsNorm->GetBinContent(5,9));

  // pi- | p+
  h_stableChargeBF->SetBinContent(13,h_stableDecayCorrelatedPairsNorm->GetBinContent(8,3) - h_stableDecayCorrelatedPairsNorm->GetBinContent(8,2));
  // K- | p+
  h_stableChargeBF->SetBinContent(14,h_stableDecayCorrelatedPairsNorm->GetBinContent(8,5) - h_stableDecayCorrelatedPairsNorm->GetBinContent(8,4));
  // p- | p+
  h_stableChargeBF->SetBinContent(15,h_stableDecayCorrelatedPairsNorm->GetBinContent(8,9) - h_stableDecayCorrelatedPairsNorm->GetBinContent(8,8));
  // pi+ | p-
  h_stableChargeBF->SetBinContent(16,h_stableDecayCorrelatedPairsNorm->GetBinContent(9,2) - h_stableDecayCorrelatedPairsNorm->GetBinContent(9,3));
  // K+ | p-
  h_stableChargeBF->SetBinContent(17,h_stableDecayCorrelatedPairsNorm->GetBinContent(9,4) - h_stableDecayCorrelatedPairsNorm->GetBinContent(9,5));
  // p+ | p-
  h_stableChargeBF->SetBinContent(18,h_stableDecayCorrelatedPairsNorm->GetBinContent(9,8) - h_stableDecayCorrelatedPairsNorm->GetBinContent(9,9));


  for (int i1=1; i1<= n; i1++)
  {
  double v1th = h_stableThermalYields->GetBinContent(i1);
  double v1   = h_stableDecayYields->GetBinContent(i1);

  for (int i2=1; i2<= n; i2++)
    {
    double v2th = h_stableThermalYields->GetBinContent(i2);
    double v2   = h_stableDecayYields->GetBinContent(i2);

    h_rho1rho1->SetBinContent(i1,i2,v1*v2);         h_rho1rho1->SetBinError(i1,i2,zero);
    h_rho1thrho1th->SetBinContent(i1,i2,v1th*v2th); h_rho1thrho1th->SetBinError(i1,i2,zero);

    double c2 = h_stableDecayCorrelatedPairsYields->GetBinContent(i1,i2);
    h_rho2->SetBinContent(i1,i2,c2+v1*v2);        h_rho2->SetBinError(i1,i2,zero);
    h_C2  ->SetBinContent(i1,i2,c2);  h_C2->SetBinError(i1,i2,zero);
    }
  }
  h_R2->Divide(h_C2,h_rho1rho1);

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
  h_allYieldsVsMass           = loadProfile(inputFile,bn+TString("allYieldsVsMass"));
  h_allToAntiRatio            = loadH1(inputFile,bn+TString("allYieldsToAnti"));
  h_stableThermalYields       = loadH1(inputFile,bn+TString("stableThermalYields"));
  h_stableThermalVsMass       = loadProfile(inputFile,bn+TString("stableThermalYieldsVsMass"));
  h_stableThermalYieldsToPion = loadH1(inputFile,bn+TString("stableThermalYieldsToPi"));
  h_stableThermalToAntiRatio  = loadH1(inputFile,bn+TString("stableThermalYieldsToAnti"));
  h_stableDecayYields         = loadH1(inputFile,bn+TString("stableDecayYields"));
  h_stableDecayVsMass         = loadProfile(inputFile,bn+TString(  "stableDecayYieldsVsMass"));
  h_stableDecayYieldsToPion   = loadH1(inputFile,bn+TString("stableDecayYieldsToPi"));
  h_stableDecayToAntiRatio    = loadH1(inputFile,bn+TString("stableDecayYieldsToAnti"));
  h_stableDecayPairsYields             = loadH2(inputFile,bn+TString("stableDecayPairYields"));
  h_stableDecayCorrelatedPairsYields   = loadH2(inputFile,bn+TString("stableDecayCorrelatedPairYields"));
  h_stableChargeBF                     = loadH2(inputFile,bn+TString("stableDecayChargeBF"));

  h_rho2            = loadH2(inputFile,bn+TString("rho2"));
  h_rho1rho1        = loadH2(inputFile,bn+TString("rho1rho1"));
  h_rho1thrho1th    = loadH2(inputFile,bn+TString("rho1thrho1th"));
  h_C2              = loadH2(inputFile,bn+TString("C2"));
  h_R2              = loadH2(inputFile,bn+TString("R2"));

}


void HadronGasParticleHistograms::fill(HadronGas & hadronGas)
{
  int nTypes       = hadronGas.getNHadrons();
  int nStableTypes = hadronGas.getNStableHadrons();
  double thermalYield1, thermalYield2;
  double decayYield1,   decayYield2;
  double mass;
  double g;

  for (int iType=0; iType<nTypes; iType++)
  {
  thermalYield1 = hadronGas.getParticle(iType)->getParticleYield();
  h_allYields->SetBinContent(iType+1, thermalYield1); h_allYields->SetBinError(iType+1, 0.0);
  mass = hadronGas.getParticle(iType)->getMass();
  g    = hadronGas.getParticle(iType)->getSpinfactor();
  h_allYieldsVsMass->Fill(mass,thermalYield1/g);
  }

  for (int iType1=0; iType1<nStableTypes; iType1++)
  {
  Hadron * part1 = hadronGas.getStableParticle(iType1);
  thermalYield1  = part1->getParticleYield();
  decayYield1    = part1->getStableParticleYield();
  h_stableThermalYields->SetBinContent(iType1+1, thermalYield1); h_stableThermalYields->SetBinError(iType1, 0.0);
  h_stableDecayYields  ->SetBinContent(iType1+1, decayYield1);   h_stableDecayYields  ->SetBinError(iType1, 0.0);
  mass = hadronGas.getStableParticle(iType1)->getMass();
  g    = hadronGas.getStableParticle(iType1)->getSpinfactor();

  h_stableThermalVsMass->Fill(mass,thermalYield1/g);
  h_stableDecayVsMass->Fill(mass,decayYield1);

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

