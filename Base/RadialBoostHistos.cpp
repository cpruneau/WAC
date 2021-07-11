//
//  RadialBoostHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "RadialBoostHistos.hpp"
ClassImp(RadialBoostHistos);

RadialBoostHistos::RadialBoostHistos(const TString & collectionName,
                                     RadialBoostConfiguration * radialBoostConfiguration,
                                     LogLevel  debugLevel)
:
Histograms(collectionName,radialBoostConfiguration,10,debugLevel)
{
  // no ops
}

RadialBoostHistos::~RadialBoostHistos()
{
  //deleteHistograms();
}

// for now use the same boundaries for eta and y histogram
void RadialBoostHistos::createHistograms()
{
  RadialBoostConfiguration & ac = *(RadialBoostConfiguration*) getConfiguration();
  TString bn = getHistoBaseName();
  h_phi  = createHistogram(bn+TString("ns_phi"),   ac.nBins_phi,  ac.min_phi,  ac.max_phi,  "#phi",      "N_{s}", 1,1,1,1);
  h_r    = createHistogram(bn+TString("ns_r"),     ac.nBins_r,    ac.min_r,    ac.max_r,    "r",         "N_{s}", 1,1,1,1);
  h_beta = createHistogram(bn+TString("ns_beta"),  ac.nBins_beta, ac.min_beta, ac.max_beta, "#beta_{s}", "N_{s}", 1,1,1,1);
  beta_mult = createProfile(bn+TString("beta_vs_mult"),  ac.nBins_mult, ac.min_mult, ac.max_mult, "mult", "#beta_{s}", 1,1,1);
  h_betaxs = createHistogram(bn+TString("n_betaxs"),       ac.nBins_beta,  ac.min_beta,  ac.max_beta,  "#beta","1/#beta dN/d#beta", 1,1,1,1);
}

//________________________________________________________________________
void RadialBoostHistos::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportFatal()) cout << "-Fatal- Attempting to load RadialBoostHistos from an invalid file pointer" << endl;
    return;
    }
  //RadialBoostConfiguration & ac = *(RadialBoostConfiguration*) getConfiguration();
  TString bn = getHistoBaseName();
  h_phi  = loadH1(inputFile,bn+TString("ns_phi")  );
  h_r    = loadH1(inputFile,bn+TString("ns_r")    );
  h_beta = loadH1(inputFile,bn+TString("ns_beta") );
  beta_mult = loadProfile(inputFile,bn+TString("beta_vs_mult") );
  h_betaxs = loadH1(inputFile,bn+TString("n_betaxs") );

  return;
}

void RadialBoostHistos::fill(double r, double phi, double beta, double weight)
{
  h_phi  ->Fill(phi,  weight);
  h_r    ->Fill(r,    weight);
  h_beta ->Fill(beta, weight);
}

void RadialBoostHistos::fill(double r, double phi, double beta, double mult, double weight)
{
  h_phi  ->Fill(phi,  weight);
  h_r    ->Fill(r,    weight);
  h_beta ->Fill(beta, weight);
  if(beta != 0)  h_betaxs ->Fill(beta, weight/beta);
  beta_mult ->Fill(mult,beta);
}





