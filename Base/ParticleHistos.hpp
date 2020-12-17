// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_ParticleHistos
#define WAC_ParticleHistos
#include "Histograms.hpp"
#include "Particle.hpp"
#include "ParticleAnalyzerConfiguration.hpp"
#include "TLorentzVector.h"


class ParticleHistos : public Histograms
{
public:

  ParticleHistos(const TString & collectionName,
                 ParticleAnalyzerConfiguration * analysisConfiguration,
                  LogLevel  debugLevel);
  virtual ~ParticleHistos();
  void createHistograms();
  void loadHistograms(TFile * inputFile);
  void fill(Particle & particle, double weight);
  void fill(TLorentzVector & p, double weight);
  void completeFill();
  void fillMultiplicity(double nAccepted, double totalEnergy, double weight);
  void calculateAverages();

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  TH1 * h_n1;
  TH1 * h_n1_eTotal;
  TH1 * h_n1_pt;
  TH1 * h_n1_ptXS;  // 1/pt dN/dptdy
  TH1 * h_n1_eta;
  TH1 * h_n1_phi;
  TH2 * h_n1_ptEta;
  TH2 * h_n1_phiEta;
  TH3 * h_n1_ptPhiEta;
  TH1 * h_spt_phi;
  TH1 * h_spt_eta;
  TH2 * h_spt_phiEta;
  TH1 * h_pt_phi;
  TH1 * h_pt_eta;
  TH2 * h_pt_phiEta;

  TH1 * h_n1_y;
  TH2 * h_n1_ptY;
  TH2 * h_n1_phiY;
  TH3 * h_n1_ptPhiY;
  TH1 * h_spt_y;
  TH2 * h_spt_phiY;
  TH1 * h_pt_y;
  TH2 * h_pt_phiY;


    ClassDef(ParticleHistos,0)

};

#endif /* WAC_ParticleHistos  */



