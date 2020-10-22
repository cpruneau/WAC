//
// ParticlePairDerivedHistos.hpp
// MyMC
//
// created by Claude Pruneau on 9/23/16.
// Copyright © 2016 Claude Pruneau. All rights reserved.
//

#ifndef WAC_ParticlePairDerivedHistos
#define WAC_ParticlePairDerivedHistos
#include "Histograms.hpp"
#include "ParticleHistos.hpp"
#include "ParticlePairHistos.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Two-Particle Correlation Functions
//////////////////////////////////////////////////////////////////////////////////////////////////////
class ParticlePairDerivedHistos : public  Histograms
{
public:

  ParticlePairDerivedHistos(const TString & name,
                            AnalysisConfiguration * configuration,
                            LogLevel  debugLevel);
  ParticlePairDerivedHistos(TFile * inputFile,
                            const TString & name,
                            AnalysisConfiguration * configuration,
                            LogLevel  debugLevel);
  virtual ~ParticlePairDerivedHistos();
  void loadHistograms(TFile * inputFile);
  void createHistograms();
  double calculateWeightCorrection(TFile * fRaw, TFile * fCal);
  void calculateDerivedHistograms(ParticleHistos * part1Histos,
                                  ParticleHistos * part2Histos,
                                  ParticlePairHistos * pairHistos,
                                  double bincorrection);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data Members - Computed
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  TH2* h_n1n1_phiEtaPhiEta;
  TH2* h_n1n1_etaEta;
  TH2* h_n1n1_phiPhi;
  TH2* h_n1n1_ptPt;

  TH2* h_pt1pt1_phiEtaPhiEta;
  TH2* h_pt1pt1_etaEta;
  TH2* h_pt1pt1_phiPhi;


  TH2* h_R2_phiEtaPhiEta;
  TH2* h_R2_etaEta;
  TH2* h_R2_phiPhi;
    TH2* h_P2_phiPhi;
  TH2* h_R2_ptPt;

  TH2* h_DptDpt_phiEtaPhiEta;
  TH2* h_DptDpt_etaEta;
  TH2* h_DptDpt_phiPhi;

  TH2* h_G2_phiEtaPhiEta;
  TH2* h_G2_etaEta;
  TH2* h_G2_phiPhi;

  /* the balance functions components not R2 based */
  TH2* h_bf12_phiEtaPhiEta;
  TH2* h_bf12_etaEta;
  TH2* h_bf12_phiPhi;
  TH2* h_bf12_phiYPhiY;
  TH2* h_bf12_yY;
  TH2* h_bf12_DetaDphi;
  TH2* h_bf12_DyDphi;
  TH2* h_bf12_DetaDphi_shft;
  TH2* h_bf12_DyDphi_shft;
  TH2* h_bf21_phiEtaPhiEta;
  TH2* h_bf21_etaEta;
  TH2* h_bf21_phiPhi;
  TH2* h_bf21_phiYPhiY;
  TH2* h_bf21_yY;
  TH2* h_bf21_DetaDphi;
  TH2* h_bf21_DyDphi;
  TH2* h_bf21_DetaDphi_shft;
  TH2* h_bf21_DyDphi_shft;

  TH2* h_n1n1_DetaDphi;
  TH2* h_pt1pt1_DetaDphi;
  TH2* h_n2_DetaDphi;
  TH2* h_ptn_DetaDphi;
  TH2* h_npt_DetaDphi;
  TH2* h_ptpt_DetaDphi;

  TH2* h_R2_DetaDphi;
  TH2* h_DptDpt_DetaDphi;
  TH2* h_P2_DetaDphi;
  TH2* h_P2_etaEta;
  TH2* h_G2_DetaDphi;

  TH2* h_n2_DetaDphi_shft;
  TH2* h_ptpt_DetaDphi_shft;
  TH2* h_R2_DetaDphi_shft;
  TH2* h_DptDpt_DetaDphi_shft;
  TH2* h_P2_DetaDphi_shft;
  TH2* h_G2_DetaDphi_shft;

  // vs Y ////////
  TH2* h_n1n1_phiYPhiY;
  TH2* h_n1n1_yY;

  TH2* h_pt1pt1_phiYPhiY;
  TH2* h_pt1pt1_yY;

  TH2* h_R2_phiYPhiY;
  TH2* h_R2_yY;

  TH2* h_DptDpt_phiYPhiY;
  TH2* h_DptDpt_yY;

  TH2* h_G2_phiYPhiY;
  TH2* h_G2_yY;
  TH2* h_P2_yY;

  TH2* h_n1n1_DyDphi;
  TH2* h_pt1pt1_DyDphi;
  TH2* h_n2_DyDphi;
  TH2* h_ptn_DyDphi;
  TH2* h_npt_DyDphi;
  TH2* h_ptpt_DyDphi;

  TH2* h_R2_DyDphi;
  TH2* h_DptDpt_DyDphi;
  TH2* h_P2_DyDphi;
  TH2* h_G2_DyDphi;

  TH2* h_n2_DyDphi_shft;
  TH2* h_R2_DyDphi_shft;
  TH2* h_DptDpt_DyDphi_shft;
  TH2* h_P2_DyDphi_shft;
  TH2* h_G2_DyDphi_shft;


  TH2 * h_n2_Q3D_xy;
  TH2 * h_n2_Q3D_xz;
  TH2 * h_n2_Q3D_yz;
  TH1 * h_n2_Q3D_x;
  TH1 * h_n2_Q3D_y;
  TH1 * h_n2_Q3D_z;

  TH3* h_n1n1_Q3D;
  TH2* h_n1n1_Q3D_xy;
  TH2* h_n1n1_Q3D_xz;
  TH2* h_n1n1_Q3D_yz;

  TH3* h_R2_Q3D;
  TH2* h_R2_Q3D_xy;
  TH2* h_R2_Q3D_xz;
  TH2* h_R2_Q3D_yz;
  TH1 * h_R2_Q3D_x;
  TH1 * h_R2_Q3D_y;
  TH1 * h_R2_Q3D_z;
  
  bool ijNormalization;

    ClassDef(ParticlePairDerivedHistos,1)

};

#endif /* WAC_ParticlePairDerivedHistos */



/*
 double v1, v2, binCorrection;
 v1 = h_intR2Nw_M->GetBinContent(1);
 v2 = h_intR2_M->GetBinContent(1);
 cout << "==========================================" << endl;
 cout << "v1=" << v1 << endl;
 cout << "v2=" << v2 << endl;

 if (v2>0)
 binCorrection = v1/v2;
 else
 binCorrection = 1;
 double invEvt = 1/nEv;
 cout << "binCorrection=" << binCorrection << endl;

 //calculateAveragePt(corr->h_s1pt_vsM, corr->h_n1_vsM, h_pt_vsM);
 //calculateAveragePt(corr->h_s1pt_vsM, corr->h_n1_vsM, h_pt_vsM);
 //cout << "calculateAveragePt(_s1pt_vsEtaPhi...)" << endl;
 */

/*
 //merging correction
 double v, ev;
 h_reverseCorrection = 0;
 cout << "  h_reverseCorrection: " << h_reverseCorrection << endl;
 if (_reverseCorrection)
 {
 for (int k=1;k<10;n++)
 {
 for (int iDeta=15;iDeta<=25;iDeta++)
 {
 v = h_R2_DetaDphi->GetBinContent(iDeta,73-k); ev = h_R2_DetaDphi->GetBinError(iDeta,73-k);
 h_R2_DetaDphi->SetBinContent(iDeta,k+1,v); h_R2_DetaDphi->SetBinError(iDeta,k+1,ev);
 v = h_DptDpt_DetaDphi->GetBinContent(iDeta,73-k); ev = h_DptDpt_DetaDphi->GetBinError(16,73-k);
 h_DptDpt_DetaDphi->SetBinContent(iDeta,k+1,v); h_DptDpt_DetaDphi->SetBinError(iDeta,k+1,ev);
 v = h_G2_DetaDphi->GetBinContent(iDeta,73-k); ev = h_G2_DetaDphi->GetBinError(iDeta,73-k);
 h_G2_DetaDphi->SetBinContent(iDeta,k+1,v); h_G2_DetaDphi->SetBinError(iDeta,k+1,ev);
 }
 }
 }
 else
 {
 for (int k=1;k<10;n++)
 {
 for (int iDeta=15;iDeta<=25;iDeta++)
 {
 v = h_R2_DetaDphi->GetBinContent(iDeta,k+1); ev = h_R2_DetaDphi->GetBinError(iDeta,k+1);
 h_R2_DetaDphi->SetBinContent(iDeta,73-k,v); h_R2_DetaDphi->SetBinError(iDeta,73-k,ev);
 v = h_DptDpt_DetaDphi->GetBinContent(iDeta,k+1); ev = h_DptDpt_DetaDphi->GetBinError(16,k+1);
 h_DptDpt_DetaDphi->SetBinContent(iDeta,73-k,v); h_DptDpt_DetaDphi->SetBinError(iDeta,73-k,ev);
 v = h_G2_DetaDphi->GetBinContent(iDeta,k+1); ev = h_G2_DetaDphi->GetBinError(iDeta,k+1);
 h_G2_DetaDphi->SetBinContent(iDeta,73-k,v); h_G2_DetaDphi->SetBinError(iDeta,73-k,ev);
 }
 }
 }
 
 double vv0, vv1;
 for (int iDeta=15;iDeta<=25;iDeta++)
 {
 vv0 = h_R2_DetaDphi->GetBinContent(iDeta,1);
 vv1 = h_R2_DetaDphi->GetBinContent(iDeta,2);
 cout << "iDeta:" << iDeta << " vv0:" << vv0 << " vv1:" << vv1 << endl;
 if (vv1>vv0 || iDeta==20) h_R2_DetaDphi->SetBinContent(iDeta,1, vv1);
 }
 */
