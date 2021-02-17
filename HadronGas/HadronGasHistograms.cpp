//
//  HadronGasHistograms.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "HadronGasHistograms.hpp"

ClassImp(HadronGasHistograms);

HadronGasHistograms::HadronGasHistograms(const TString & name,
                                         HadronGasConfiguration * configuration,
                                         LogLevel  debugLevel)
:
Histograms(name,configuration,100,debugLevel),
sysEnergyVsT(nullptr),
sysEntropyVsT(nullptr),
sysPressureVsT(nullptr),
sysNetBVsT(nullptr),
sysNetQVsT(nullptr),
sysNetBVsTMuB(nullptr),
sysNetQVsTMuB(nullptr)
{
  // no ops
}

HadronGasHistograms::~HadronGasHistograms()
{
  //deleteHistograms();
}

void HadronGasHistograms::createHistograms()
{
  HadronGasConfiguration & ac = *(HadronGasConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();
  int n = 1000*(ac.maxT-ac.minT);
  sysEnergyVsT   = createHistogram(bn+TString("sysEnergyVsT"),  n, ac.minT, ac.maxT, "T (GeV)","e (GeV/fm^{3})", false,false,false,false);
  sysEntropyVsT  = createHistogram(bn+TString("sysEntropyVsT"), n, ac.minT, ac.maxT, "T (GeV)","s (fm^{-3})",    false,false,false,false);
  sysPressureVsT = createHistogram(bn+TString("sysPressureVsT"),n, ac.minT, ac.maxT, "T (GeV)","p",              false,false,false,false);
  sysNetBVsT     = createHistogram(bn+TString("sysNetBVsT"),    n, ac.minT, ac.maxT, "T (GeV)","b (fm^{-3})",    false,false,false,false);
  sysNetQVsT     = createHistogram(bn+TString("sysNetQVsT"),    n, ac.minT, ac.maxT, "T (GeV)","q (fm^{-3})",    false,false,false,false);
  sysNetBVsTMuB  = createHistogram(bn+TString("sysNetBVsTMuB"), n, ac.minT, ac.maxT, ac.nMu, ac.minMu, ac.maxMu, "T (GeV)","#mu (GeV)","b (fm^{-3})", false,false,false,false);
  sysNetQVsTMuB  = createHistogram(bn+TString("sysNetQVsTMuB"), n, ac.minT, ac.maxT, ac.nMu, ac.minMu, ac.maxMu, "T (GeV)","#mu (GeV)","q (fm^{-3})", false,false,false,false);
}


//________________________________________________________________________
void HadronGasHistograms::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportFatal()) cout << "Attempting to load HadronGasHistograms from an invalid file pointer" << endl;
    return;
    }
  TString bn  = getHistoBaseName();
  sysEnergyVsT = loadH1(inputFile,bn+TString("sysEnergyVsT"));
  if (!sysEnergyVsT)
    {
    if (reportError()) cout << "Could not load histogram: " << bn+TString("sysEnergyVsT") << endl;
    return;
    }
  sysEntropyVsT  = loadH1(inputFile,bn+TString("sysEntropyVsT"));
  sysPressureVsT = loadH1(inputFile,bn+TString("sysPressureVsT"));
  sysNetBVsT     = loadH1(inputFile,bn+TString("sysNetBVsT"));
  sysNetQVsT     = loadH1(inputFile,bn+TString("sysNetQVsT"));
  sysNetBVsTMuB  = loadH2(inputFile,bn+TString("sysNetBVsTMuB"));
  sysNetQVsTMuB  = loadH2(inputFile,bn+TString("sysNetQVsTMuB"));
}

void HadronGasHistograms::fill(HadronGas & hadronGas)
{
  double temperature = hadronGas.getTemperature();
  double muB = hadronGas.getMuB();
  int iT = sysEnergyVsT->GetXaxis()->FindBin(temperature);
  //sysEnergyVsT   ->Fill(temperature, hadronGas.getEnergyDensity() ); sysEnergyVsT->SetBinError(iT,0.0);
//  sysEntropyVsT  ->Fill(temperature, hadronGas.getEntropyDensity()); sysEntropyVsT->SetBinError(iT,0.0);
//  sysPressureVsT ->Fill(temperature, hadronGas.getPressure()); sysPressureVsT->SetBinError(iT,0.0);
//  sysNetBVsT     ->Fill(temperature, hadronGas.getNetBaryonDensity()); sysNetBVsT->SetBinError(iT,0.0);
//  sysNetQVsT     ->Fill(temperature, hadronGas.getNetStrangenessDensity()); sysNetQVsT->SetBinError(iT,0.0);
//  sysNetBVsTMuB  ->Fill(temperature, muB, hadronGas.getNetBaryonDensity()); sysNetBVsTMuB->SetBinError(iT,0.0);
//  sysNetQVsTMuB  ->Fill(temperature, muB, hadronGas.getNetChargeDensity()); sysNetQVsTMuB->SetBinError(iT,0.0);

  sysEnergyVsT   ->SetBinContent(iT, hadronGas.getEnergyDensity() ); sysEnergyVsT->SetBinError(iT,0.0);
  sysEntropyVsT  ->SetBinContent(iT, hadronGas.getEntropyDensity()); sysEntropyVsT->SetBinError(iT,0.0);
  sysPressureVsT ->SetBinContent(iT, hadronGas.getPressure()); sysPressureVsT->SetBinError(iT,0.0);
  sysNetBVsT     ->SetBinContent(iT, hadronGas.getNetBaryonDensity()); sysNetBVsT->SetBinError(iT,0.0);
  sysNetQVsT     ->SetBinContent(iT, hadronGas.getNetStrangenessDensity()); sysNetQVsT->SetBinError(iT,0.0);
  sysNetBVsTMuB  ->SetBinContent(iT, muB, hadronGas.getNetBaryonDensity()); sysNetBVsTMuB->SetBinError(iT,0.0);
  sysNetQVsTMuB  ->SetBinContent(iT, muB, hadronGas.getNetChargeDensity()); sysNetQVsTMuB->SetBinError(iT,0.0);


}

