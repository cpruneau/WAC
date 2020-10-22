// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining ParticlePairHistos
 */

#ifndef WAC_ParticlePairHistos
#define WAC_ParticlePairHistos

#include "Histograms.hpp"

class ParticlePairHistos : public Histograms
{
public:

  ParticlePairHistos(const TString & name,
                     AnalysisConfiguration * configuration,
                     LogLevel  debugLevel);
  ParticlePairHistos(TFile * inputFile,
                     const TString & name,
                     AnalysisConfiguration * configuration,
                     LogLevel  debugLevel);
  virtual ~ParticlePairHistos();
  void initialize();
  void fill(Particle & particle1, Particle & particle2, double weight);
  void completeFill();
  void loadHistograms(TFile * inputFile);

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  TH2 * h_n2_ptPt;

  TH2 * h_n2_phiEtaPhiEta;
  TH2 * h_npt_phiEtaPhiEta;
  TH2 * h_ptn_phiEtaPhiEta;
  TH2 * h_ptpt_phiEtaPhiEta;
  TH2 * h_n2_ptPhiEtaPtPhiEta;

  TH2 * h_n2_etaEta;
  TH2 * h_npt_etaEta;
  TH2 * h_ptn_etaEta;
  TH2 * h_ptpt_etaEta;

  TH2 * h_n2_phiPhi;
  TH2 * h_ptn_phiPhi;
  TH2 * h_npt_phiPhi;
  TH2 * h_ptpt_phiPhi;
  
  TH2 * h_n2_yY;
  TH2 * h_n2_phiYPhiY;
  TH2 * h_npt_phiYPhiY;
  TH2 * h_ptn_phiYPhiY;
  TH2 * h_ptpt_phiYPhiY;
  TH2 * h_n2_ptPhiYPtPhiY;
  
  TH2 * h_npt_yY;
  TH2 * h_ptn_yY;
  TH2 * h_ptpt_yY;
  
  TH3 * h_n2_Q3D;   // number of pairs vs Q3D
  TH1 * h_mInv_Lab; // lab fram calculation
  TH1 * h_mInv_PF;  // pair frame
  TH1 * h_beta;

     ClassDef(ParticlePairHistos,0)
};


#endif /* WAC_ParticlePairHistos  */
