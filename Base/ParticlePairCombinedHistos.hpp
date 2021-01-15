//
//  ParticlePairCombinedHistos.hpp
//  MyMC
//
//  Created by Claude Pruneau on 12/18/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#ifndef WAC_ParticlePairCombinedHistos
#define WAC_ParticlePairCombinedHistos
#include "ParticlePairDerivedHistos.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Two-Particle Correlation Combinations (CI, CD)
//////////////////////////////////////////////////////////////////////////////////////////////////////
class ParticlePairCombinedHistos : public Histograms
{
public:

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  ParticlePairCombinedHistos(const TString & name,
                             ParticlePairAnalyzerConfiguration * configuration,
                             LogLevel  debugLevel);
  virtual ~ParticlePairCombinedHistos();
  virtual void loadHistograms(TFile * inputFile);
  virtual void createHistograms();
  virtual void calculate(ParticlePairDerivedHistos * pp, ParticlePairDerivedHistos * mm, ParticlePairDerivedHistos * pm, double app, double amm, double apm);
  //void calculateRR_Q3D(ClusterHistos * ch);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data Members - Computed
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH2* h_R2_etaEta;
  TH2* h_R2_phiPhi;
  TH2* h_R2_ptPt;
  TH2* h_R2_DetaDphi_shft;
  TH2* h_R2_yY;
  TH2* h_R2_DyDphi_shft;

  TH2* h_DptDpt_etaEta;
  TH2* h_DptDpt_phiPhi;
  TH2* h_DptDpt_DetaDphi_shft;
  TH2* h_DptDpt_yY;
  TH2* h_DptDpt_DyDphi_shft;

  TH2* h_P2_etaEta;
  TH2* h_P2_phiPhi;
  TH2* h_P2_DetaDphi_shft;
  TH2* h_P2_yY;
  TH2* h_P2_DyDphi_shft;

  TH2* h_G2_etaEta;
  TH2* h_G2_phiPhi;
  TH2* h_G2_DetaDphi_shft;
  TH2* h_G2_yY;
  TH2* h_G2_DyDphi_shft;

  TH2* h_BF_etaEta;
  TH2* h_BF_phiPhi;
  TH2* h_BF_DetaDphi_shft;
  TH2* h_BF_yY;
  TH2* h_BF_DyDphi_shft;

  TH3* h_R2_Q3D;
  TH2* h_R2_Q3D_xy;
  TH2* h_R2_Q3D_xz;
  TH2* h_R2_Q3D_yz;
  TH1 * h_R2_Q3D_x;
  TH1 * h_R2_Q3D_y;
  TH1 * h_R2_Q3D_z;
  
  TH3* h_RR2_Q3D;
  TH2* h_RR2_Q3D_xy;
  TH2* h_RR2_Q3D_xz;
  TH2* h_RR2_Q3D_yz;
  TH1 * h_RR2_Q3D_x;
  TH1 * h_RR2_Q3D_y;
  TH1 * h_RR2_Q3D_z;
  
    ClassDef(ParticlePairCombinedHistos,1)
};



#endif /* WAC_ParticlePairCombinedHistos */
