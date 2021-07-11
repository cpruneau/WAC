// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_RadialBoostHistos
#define WAC_RadialBoostHistos
#include "Histograms.hpp"
#include "RadialBoostConfiguration.hpp"


class RadialBoostHistos : public Histograms
{
public:

  RadialBoostHistos(const TString & collectionName,
                    RadialBoostConfiguration * radialBoostConfiguration,
                    LogLevel  debugLevel);
  virtual ~RadialBoostHistos();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void fill(double r, double phi, double beta, double weight);
  virtual void fill(double r, double phi, double beta, double mult, double weight);

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  TH1 * h_phi;
  TH1 * h_r;
  TH1 * h_beta;
  TH1 * h_betaxs;
  TProfile * beta_mult;

  

  ClassDef(RadialBoostHistos,0)

};

#endif /* WAC_RadialBoostHistos  */



