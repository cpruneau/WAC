// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */

#include "ParticlePairHistos.hpp"

ClassImp(ParticlePairHistos);

ParticlePairHistos::ParticlePairHistos(const TString & name,
                                       ParticlePairAnalyzerConfiguration * configuration,
                                       LogLevel  debugLevel)
:
Histograms(name,configuration,150,debugLevel)
{
  // no ops
}

ParticlePairHistos::~ParticlePairHistos()
{
  // no ops
}

void ParticlePairHistos::createHistograms()
{
  if (reportDebug()) cout << "ParticlePairHistos::initialize() Started." << endl;
  ParticlePairAnalyzerConfiguration & ac = *(ParticlePairAnalyzerConfiguration*) getConfiguration();
  TString bn = getHistoBaseName();
  ac.range_pt       = ac.max_pt  - ac.min_pt;
  ac.range_phi      = ac.max_phi - ac.min_phi;
  ac.range_eta      = ac.max_eta - ac.min_eta;
  ac.range_y        = ac.max_y   - ac.min_y;
  ac.nBins_phiEta   = ac.nBins_eta * ac.nBins_phi;
  ac.nBins_phiEtaPt = ac.nBins_eta * ac.nBins_phi * ac.nBins_pt;
  ac.nBins_phiY     = ac.nBins_y   * ac.nBins_phi;
  ac.nBins_phiYPt   = ac.nBins_y   * ac.nBins_phi * ac.nBins_pt;

  h_n2_ptPt           = createHistogram(bn+TString("n2_ptPt"),          ac.nBins_pt,     ac.min_pt, ac.max_pt, ac.nBins_pt, ac.min_pt, ac.max_pt,      "p_{T,1}",     "p_{T,2}", "N_{2}", 1,1,1,0);
  h_n2_etaEta         = createHistogram(bn+TString("n2_etaEta"),        ac.nBins_eta,    ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}",    "#eta_{2}", "N_{2}", 1,1,1,0);
  h_n2_phiPhi         = createHistogram(bn+TString("n2_phiPhi"),        ac.nBins_phi,    ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "N_{2}", 1,1,1,0);
#ifdef OPTIMIZEADDBINCONTENT
  /* big histograms are forced to be created without sumw2 structure for it will not be used */
  bool defsumw2 = TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2(false);
#endif // OPTIMIZEADDBINCONTENT
  h_n2_phiEtaPhiEta   = createHistogram(bn+TString("n2_phiEtaPhiEta"),  ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}","#eta_{2}x#varphi_{2}", "N_{2}", 1,1,1,0);
  h_npt_phiEtaPhiEta  = createHistogram(bn+TString("npt_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}","#eta_{2}x#varphi_{2}", "Nxp_{T}", 1,1,1,0);
  h_ptn_phiEtaPhiEta  = createHistogram(bn+TString("ptn_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}","#eta_{2}x#varphi_{2}", "p_{T}xN", 1,1,1,0);
  h_ptpt_phiEtaPhiEta = createHistogram(bn+TString("ptpt_phiEtaPhiEta"),ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}","#eta_{2}x#varphi_{2}", "p_{T}xp_{T}", 1,1,1,0);
#ifdef OPTIMIZEADDBINCONTENT
  /* big histograms are forced to be created without sumw2 structure for it will not be used */
  h_n2_phiEtaPhiEta->SetBit(TH1::kIsNotW);
  h_n2_phiEtaPhiEta->Sumw2(false);
  h_npt_phiEtaPhiEta->SetBit(TH1::kIsNotW);
  h_npt_phiEtaPhiEta->Sumw2(false);
  h_ptn_phiEtaPhiEta->SetBit(TH1::kIsNotW);
  h_ptn_phiEtaPhiEta->Sumw2(false);
  h_ptpt_phiEtaPhiEta->SetBit(TH1::kIsNotW);
  h_ptpt_phiEtaPhiEta->Sumw2(false);
  TH1::SetDefaultSumw2(defsumw2);
#endif // OPTIMIZEADDBINCONTENT
  h_npt_etaEta        = createHistogram(bn+TString("npt_etaEta"),       ac.nBins_eta,    ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}",    "#eta_{2}", "n x p_{T}", 1,1,0,0);
  h_ptn_etaEta        = createHistogram(bn+TString("ptn_etaEta"),       ac.nBins_eta,    ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}",    "#eta_{2}", "p_{T} x n", 1,1,0,0);
  h_ptpt_etaEta       = createHistogram(bn+TString("ptpt_etaEta"),      ac.nBins_eta,    ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}",    "#eta_{2}", "p_{T}xp_{T}", 1,1,0,0);
  h_npt_phiPhi        = createHistogram(bn+TString("npt_phiPhi"),       ac.nBins_phi,    ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "n x p_{T}",1,1,0,0);
  h_ptn_phiPhi        = createHistogram(bn+TString("ptn_phiPhi"),       ac.nBins_phi,    ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "p_{T} x n",1,1,0,0);
  h_ptpt_phiPhi       = createHistogram(bn+TString("ptpt_phiPhi"),      ac.nBins_phi,    ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "p_{T}xp_{T}",1,1,0,0);

  if (false && ac.fillQ3D)
    {
    h_n2_Q3D   = createHistogram(bn+TString("n2_Q3D"),ac.nBins_DeltaPlong, ac.min_DeltaPlong, ac.max_DeltaPlong,
                                 ac.nBins_DeltaPside, ac.min_DeltaPside, ac.max_DeltaPside,ac.nBins_DeltaPout,  ac.min_DeltaPout,  ac.max_DeltaPout,
                                 "Q_{long} (GeV/c)", "Q_{side}  (GeV/c)", "Q_{out}  (GeV/c)","Yield", 1,1,1,0);
    h_mInv_Lab = createHistogram(bn+TString("n2_mInvLab"),  100, static_cast<double>(0.0), 1.0, "m_{inv} (GeV/c^{2})", "yield", 1,1,1,0);
    h_beta     = createHistogram(bn+TString("n2_beta"),  100, static_cast<double>(0.0), 1.0, "#beta", "yield", 1,1,1,0);
    }
  if (ac.fillY)
    {
    h_n2_yY         = createHistogram(bn+TString("n2_yY"),        ac.nBins_y,    ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}","y_{2}", "N_{2}", 1,1,0,0);
#ifdef OPTIMIZEADDBINCONTENT
    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    bool defsumw2 = TH1::GetDefaultSumw2();
    TH1::SetDefaultSumw2(false);
#endif // OPTIMIZEADDBINCONTENT
    h_n2_phiYPhiY   = createHistogram(bn+TString("n2_phiYPhiY"),  ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}","y_{2}x#varphi_{2}", "N_{2}", 1,1,0,0);
    h_npt_phiYPhiY  = createHistogram(bn+TString("npt_phiYPhiY"), ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}","y_{2}x#varphi_{2}", "Nxp_{T}", 1,1,0,0);
    h_ptn_phiYPhiY  = createHistogram(bn+TString("ptn_phiYPhiY"), ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}","y_{2}x#varphi_{2}", "p_{T}xN", 1,1,0,0);
    h_ptpt_phiYPhiY = createHistogram(bn+TString("ptpt_phiYPhiY"),ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}","y_{2}x#varphi_{2}", "p_{T}xp_{T}", 1,1,0,0);
#ifdef OPTIMIZEADDBINCONTENT
    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    h_n2_phiYPhiY->SetBit(TH1::kIsNotW);
    h_n2_phiYPhiY->Sumw2(false);
    h_npt_phiYPhiY->SetBit(TH1::kIsNotW);
    h_npt_phiYPhiY->Sumw2(false);
    h_ptn_phiYPhiY->SetBit(TH1::kIsNotW);
    h_ptn_phiYPhiY->Sumw2(false);
    h_ptpt_phiYPhiY->SetBit(TH1::kIsNotW);
    h_ptpt_phiYPhiY->Sumw2(false);
    TH1::SetDefaultSumw2(defsumw2);
#endif // OPTIMIZEADDBINCONTENT
    h_npt_yY        = createHistogram(bn+TString("npt_yY"),       ac.nBins_y,    ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}","y_{2}", "n x p_{T}",   1,1,0,0);
    h_ptn_yY        = createHistogram(bn+TString("ptn_yY"),       ac.nBins_y,    ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}","y_{2}", "p_{T} x n",   1,1,0,0);
    h_ptpt_yY       = createHistogram(bn+TString("ptpt_yY"),      ac.nBins_y,    ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}","y_{2}", "p_{T}xp_{T}", 1,1,0,0);
    }
  if (reportDebug()) cout << "ParticlePairHistos::initialize() Completed." << endl;
}

void ParticlePairHistos::fill(Particle & particle1, Particle & particle2, double weight)
{
  //if (reportStart("ParticlePairHistos",getName(), "fill(Particle & particle1, Particle & particle2, double weight)"))
  //  ;
  double pt1   = particle1.pt;
//  double eta1  = particle1.eta;
//  double phi1  = particle1.phi; if (phi1<0) phi1 += TMath::TwoPi();

  double pt2   = particle2.pt;
//  double eta2  = particle2.eta;
//  double phi2  = particle2.phi; if (phi2<0) phi2 += TMath::TwoPi();

  ParticlePairAnalyzerConfiguration & ac = * (ParticlePairAnalyzerConfiguration*) getConfiguration();
  h_n2_ptPt    ->Fill(pt1,  pt2,  weight);
  // delayed fill h_n2_etaEta  ->Fill(eta1, eta2, weight);
  // delayed fill h_n2_phiPhi  ->Fill(phi1, phi2, weight);

  //if (reportDebug("ParticlePairHistos",getName(), "fill(Particle & particle1, Particle & particle2, double weight)")) cout << " -1- " << endl;

  // delayed fill h_ptn_etaEta ->Fill(eta1, eta2, weight*pt1);
  // delayed fill h_npt_etaEta ->Fill(eta1, eta2, weight*pt2);
  // delayed fill h_ptpt_etaEta->Fill(eta1, eta2, weight*pt1*pt2);
  // delayed fill h_ptn_phiPhi ->Fill(phi1, phi2, weight*pt1);
  // delayed fill h_npt_phiPhi ->Fill(phi1, phi2, weight*pt2);
  // delayed fill h_ptpt_phiPhi->Fill(phi1, phi2, weight*pt1*pt2);

  int iPhiEta1 = particle1.ixEtaPhi;
  int iPhiEta2 = particle2.ixEtaPhi;

  bool wrongix = iPhiEta1 < 0 || iPhiEta2 < 0;

  if (!wrongix) {


#ifdef OPTIMIZEADDBINCONTENT
    int binno = (iPhiEta2+1)*(h_n2_phiEtaPhiEta->GetNbinsX()+2)+(iPhiEta1+1);
    double nentries = h_n2_phiEtaPhiEta->GetEntries() + 1;
    h_n2_phiEtaPhiEta  ->AddBinContent(binno, weight);
    h_ptn_phiEtaPhiEta ->AddBinContent(binno, weight*pt1);
    h_npt_phiEtaPhiEta ->AddBinContent(binno, weight*pt2);
    h_ptpt_phiEtaPhiEta->AddBinContent(binno, weight*pt1*pt2);
    h_n2_phiEtaPhiEta->SetEntries(nentries);
    h_ptn_phiEtaPhiEta->SetEntries(nentries);
    h_npt_phiEtaPhiEta->SetEntries(nentries);
    h_ptpt_phiEtaPhiEta->SetEntries(nentries);
#else
    h_n2_phiEtaPhiEta  ->Fill(iPhiEta1+0.5, iPhiEta2+0.5, weight);
    h_ptn_phiEtaPhiEta ->Fill(iPhiEta1+0.5, iPhiEta2+0.5, weight*pt1);
    h_npt_phiEtaPhiEta ->Fill(iPhiEta1+0.5, iPhiEta2+0.5, weight*pt2);
    h_ptpt_phiEtaPhiEta->Fill(iPhiEta1+0.5, iPhiEta2+0.5, weight*pt1*pt2);
#endif // OPTIMIZEADDBINCONTENT
  }

  //      if (ac.fill6D)
  //      {
  //        int iPt1 = int( (pt1-ac.min_pt)/ac.range_pt );
  //        int iPt2 = int( (pt2-ac.min_pt)/ac.range_pt );
  //        double iPtPhiEta1 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iEta1 + ac.nBins_pt*iPhi1 + iPt1);
  //        double iPtPhiEta2 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iEta2 + ac.nBins_pt*iPhi2 + iPt2);
  //        h_n2_ptPhiEtaPtPhiEta->Fill(iPtPhiEta1, iPtPhiEta2, weight);
  //      }

  if (ac.fillY)
    {
//    double y1    = particle1.y;
//    double y2    = particle2.y;
    // delayed fill h_n2_yY  ->Fill(y1, y2, weight);
    // delayed fill h_ptn_yY ->Fill(y1, y2, weight*pt1);
    // delayed fill h_npt_yY ->Fill(y1, y2, weight*pt2);
    // delayed fill h_ptpt_yY->Fill(y1, y2, weight*pt1*pt2);
    int iPhiY1 = particle1.ixYPhi;
    int iPhiY2 = particle2.ixYPhi;
    wrongix = iPhiY1 < 0 || iPhiY2 < 0;
    if (!wrongix) {
#ifdef OPTIMIZEADDBINCONTENT
      int binno = (iPhiY2+1)*(h_n2_phiYPhiY->GetNbinsX()+2)+(iPhiY1+1);
      double nentries = h_n2_phiEtaPhiEta->GetEntries() + 1;
      h_n2_phiYPhiY  ->AddBinContent(binno, weight);
      h_ptn_phiYPhiY ->AddBinContent(binno, weight*pt1);
      h_npt_phiYPhiY ->AddBinContent(binno, weight*pt2);
      h_ptpt_phiYPhiY->AddBinContent(binno, weight*pt1*pt2);
      h_n2_phiYPhiY->SetEntries(nentries);
      h_ptn_phiYPhiY->SetEntries(nentries);
      h_npt_phiYPhiY->SetEntries(nentries);
      h_ptpt_phiYPhiY->SetEntries(nentries);
#else
      h_n2_phiYPhiY  ->Fill(iPhiY1+0.5, iPhiY2+0.5, weight);
      h_ptn_phiYPhiY ->Fill(iPhiY1+0.5, iPhiY2+0.5, weight*pt1);
      h_npt_phiYPhiY ->Fill(iPhiY1+0.5, iPhiY2+0.5, weight*pt2);
      h_ptpt_phiYPhiY->Fill(iPhiY1+0.5, iPhiY2+0.5, weight*pt1*pt2);
#endif // OPTIMIZEADDBINCONTENT
    }

    //        if (ac.fill6D)
    //        {
    //          int iPt1 = int( (pt1-ac.min_pt)/ac.range_pt );
    //          int iPt2 = int( (pt2-ac.min_pt)/ac.range_pt );
    //          double iPtPhiY1 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iY1 + ac.nBins_pt*iPhi1 + iPt1);
    //          double iPtPhiY2 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iY2 + ac.nBins_pt*iPhi2 + iPt2);
    //          h_n2_ptPhiYPtPhiY->Fill(iPtPhiY1, iPtPhiY2, weight);
    //        }
    }

//  if (false && ac.fillQ3D)
//    {
//    double Qlong, Qside, Qout, Qinv;
//    //double Qlong1, Qside1, Qout1;
//
//    double px1   = particle1.px;
//    double py1   = particle1.py;
//    double pz1   = particle1.pz;
//    double e1    = particle1.e;
//
//    double px2   = particle2.px;
//    double py2   = particle2.py;
//    double pz2   = particle2.pz;
//    double e2    = particle2.e;
//
//    double pt,s,Mlong,roots;
//    double ptot[4],q[4];
//    //const int g[4]={1,-1,-1,-1};
//    //int alpha;
//    Qinv=0.0;
//    s=0.0;
//    ptot[0] = e1  + e2;
//    ptot[1] = px1 + px2;
//    ptot[2] = py1 + py2;
//    ptot[3] = pz1 + pz2;
//    q[0] = e1  - e2;
//    q[1] = px1 - px2;
//    q[2] = py1 - py2;
//    q[3] = pz1 - pz2;
//    s = ptot[0]*ptot[0] - ptot[1]*ptot[1] - ptot[2]*ptot[2] - ptot[3]*ptot[3];
//    Qinv = -(q[0]*q[0] - q[1]*q[1] - q[2]*q[2] - q[3]*q[3]);
//    pt=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
//    Mlong=sqrt(s+pt*pt);
//    roots=sqrt(s);
//
//    if (pt>0)
//      {
//      Qside = (ptot[1]*q[2]-ptot[2]*q[1])/pt;
//      Qlong = (ptot[0]*q[3]-ptot[3]*q[0])/Mlong;
//      Qout  = (roots/Mlong)*(ptot[1]*q[1]+ptot[2]*q[2])/pt;
//      }
//    else
//      {
//      Qlong = q[3];
//      Qside = q[2];
//      Qout  = q[1];
//      }
//
//    h_n2_Q3D->Fill(Qlong, Qside, Qout,weight);
//    }

}

// complete filling the addicional histograms by projecting the 
// higher dimensional ones
void ParticlePairHistos::completeFill()
{
  ParticlePairAnalyzerConfiguration & ac = * (ParticlePairAnalyzerConfiguration*) getConfiguration();

  int nbinseta = h_n2_etaEta->GetNbinsX();
  int nbinsphi = h_n2_phiPhi->GetNbinsX();

  project_n2XYXY_n2XX(h_n2_phiEtaPhiEta,h_n2_etaEta,nbinseta,nbinsphi);
  project_n2XYXY_n2YY(h_n2_phiEtaPhiEta,h_n2_phiPhi,nbinseta,nbinsphi);

  project_n2XYXY_n2XX(h_ptn_phiEtaPhiEta,h_ptn_etaEta,nbinseta,nbinsphi);
  project_n2XYXY_n2XX(h_npt_phiEtaPhiEta,h_npt_etaEta,nbinseta,nbinsphi);
  project_n2XYXY_n2XX(h_ptpt_phiEtaPhiEta,h_ptpt_etaEta,nbinseta,nbinsphi);

  project_n2XYXY_n2YY(h_ptn_phiEtaPhiEta,h_ptn_phiPhi,nbinseta,nbinsphi);
  project_n2XYXY_n2YY(h_npt_phiEtaPhiEta,h_npt_phiPhi,nbinseta,nbinsphi);
  project_n2XYXY_n2YY(h_ptpt_phiEtaPhiEta,h_ptpt_phiPhi,nbinseta,nbinsphi);

  if (ac.fillY)
    {
      int nbinsy = h_n2_yY->GetNbinsX();

      project_n2XYXY_n2XX(h_n2_phiYPhiY,h_n2_yY,nbinsy,nbinsphi);
      project_n2XYXY_n2XX(h_ptn_phiYPhiY,h_ptn_yY,nbinsy,nbinsphi);
      project_n2XYXY_n2XX(h_npt_phiYPhiY,h_npt_yY,nbinsy,nbinsphi);
      project_n2XYXY_n2XX(h_ptpt_phiYPhiY,h_ptpt_yY,nbinsy,nbinsphi);
    }
}

//  void fill(TLorentzVector & p1, TLorentzVector & p2, double weight)
//  {
//
//    double pt1   = p1.Pt();
//    double eta1  = p1.Eta();
//    double phi1  = p1.Phi();
//    double y1    = p1.Rapidity();
//    double pt2   = p2.Pt();
//    double eta2  = p2.Eta();
//    double phi2  = p2.Phi();
//    double y2    = p2.Rapidity();
//
//    if (phi1<0) phi1 += TMath::TwoPi();  // sets the range to [0,2pi[
//    if (phi2<0) phi2 += TMath::TwoPi();
//
//    ParticlePairAnalyzerConfiguration & ac = *analysisConfiguration;
//    h_n2_ptPt    ->Fill(pt1,  pt2,  weight);
//    h_n2_etaEta  ->Fill(eta1, eta2, weight);
//    h_n2_phiPhi  ->Fill(phi1, phi2, weight);
//
//    h_ptn_etaEta ->Fill(eta1, eta2, weight*pt1);
//    h_npt_etaEta ->Fill(eta1, eta2, weight*pt2);
//    h_ptpt_etaEta->Fill(eta1, eta2, weight*pt1*pt2);
//    h_ptn_phiPhi ->Fill(phi1, phi2, weight*pt1);
//    h_npt_phiPhi ->Fill(phi1, phi2, weight*pt2);
//    h_ptpt_phiPhi->Fill(phi1, phi2, weight*pt1*pt2);
//
//    int iEta1 = int( double(ac.nBins_eta)*(eta1-ac.min_eta)/ac.range_eta );
//    int iEta2 = int( double(ac.nBins_eta)*(eta2-ac.min_eta)/ac.range_eta );
//    int iPhi1 = int( double(ac.nBins_phi)*(phi1-ac.min_phi)/ac.range_phi );
//    int iPhi2 = int( double(ac.nBins_phi)*(phi2-ac.min_phi)/ac.range_phi );
//
//    //if (debugLevel>0)  cout << "ParticlePairHistos::fill() iEta1:" << iEta1  << " iEta2:" << iEta2 << " iPhi1:" << iPhi1 << " iPhi2:" << iPhi2 << endl;
//
//    if (iEta1>=0 && iEta2>=0 && iPhi1>=0 && iPhi2>=0 &&
//        iEta1<ac.nBins_eta && iEta2< ac.nBins_eta && iPhi1<ac.nBins_phi && iPhi2<ac.nBins_phi)
//    {
//      double iPhiEta1 = 0.1 + double(ac.nBins_phi*iEta1 + iPhi1);
//      double iPhiEta2 = 0.1 + double(ac.nBins_phi*iEta2 + iPhi2);
//      h_n2_phiEtaPhiEta  ->Fill(iPhiEta1, iPhiEta2, weight);
//      h_ptn_phiEtaPhiEta ->Fill(iPhiEta1, iPhiEta2, weight*pt2);
//      h_npt_phiEtaPhiEta ->Fill(iPhiEta1, iPhiEta2, weight*pt1);
//      h_ptpt_phiEtaPhiEta->Fill(iPhiEta1, iPhiEta2, weight*pt1*pt2);
//
//      if (ac.fill6D)
//      {
//        int iPt1 = int( (pt1-ac.min_pt)/ac.range_pt );
//        int iPt2 = int( (pt2-ac.min_pt)/ac.range_pt );
//        double iPtPhiEta1 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iEta1 + ac.nBins_pt*iPhi1 + iPt1);
//        double iPtPhiEta2 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iEta2 + ac.nBins_pt*iPhi2 + iPt2);
//        h_n2_ptPhiEtaPtPhiEta->Fill(iPtPhiEta1, iPtPhiEta2, weight);
//      }
//    }
//
//    if (ac.fillY)
//    {
//      h_n2_yY  ->Fill(y1, y2, weight);
//      h_ptn_yY ->Fill(y1, y2, weight*pt1);
//      h_npt_yY ->Fill(y1, y2, weight*pt2);
//      h_ptpt_yY->Fill(y1, y2, weight*pt1*pt2);
//
//
//      int iY1   = int( double(ac.nBins_eta)*(y1-ac.min_eta)/ac.range_eta );
//      int iY2   = int( double(ac.nBins_eta)*(y2-ac.min_eta)/ac.range_eta );
//
//      if (iY1>=0 && iY2>=0 && iPhi1>=0 && iPhi2>=0 &&
//          iY1<ac.nBins_eta && iY2< ac.nBins_eta && iPhi1<ac.nBins_phi && iPhi2<ac.nBins_phi)
//      {
//        double iPhiY1 = 0.1 + double(ac.nBins_phi*iY1 + iPhi1);
//        double iPhiY2 = 0.1 + double(ac.nBins_phi*iY2 + iPhi2);
//        h_n2_phiYPhiY  ->Fill(iPhiY1, iPhiY2, weight);
//        h_ptn_phiYPhiY ->Fill(iPhiY1, iPhiY2, weight*pt2);
//        h_npt_phiYPhiY ->Fill(iPhiY1, iPhiY2, weight*pt1);
//        h_ptpt_phiYPhiY->Fill(iPhiY1, iPhiY2, weight*pt1*pt2);
//
//        if (ac.fill6D)
//        {
//          int iPt1 = int( (pt1-ac.min_pt)/ac.range_pt );
//          int iPt2 = int( (pt2-ac.min_pt)/ac.range_pt );
//          double iPtPhiY1 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iY1 + ac.nBins_pt*iPhi1 + iPt1);
//          double iPtPhiY2 = 0.1 + double(ac.nBins_phi*ac.nBins_pt*iY2 + ac.nBins_pt*iPhi2 + iPt2);
//          h_n2_ptPhiYPtPhiY->Fill(iPtPhiY1, iPtPhiY2, weight);
//        }
//      }
//    }
//
//    if (ac.fillQ3D)
//    {
//      double Qlong, Qside, Qout, Qinv;
//      double Qlong1, Qside1, Qout1, Qinv1;
//
//      ///
//      double pt,s,Mlong,roots;
//      double ptot[4],q[4];
//      const int g[4]={1,-1,-1,-1};
//      int alpha;
//      Qinv=0.0;
//      s=0.0;
//      ptot[0] = p1.E()  + p2.E();
//      ptot[1] = p1.Px() + p2.Px();
//      ptot[2] = p1.Py() + p2.Py();
//      ptot[3] = p1.Pz() + p2.Pz();
//      q[0] = p1.E()  - p2.E();
//      q[1] = p1.Px() - p2.Px();
//      q[2] = p1.Py() - p2.Py();
//      q[3] = p1.Pz() - p2.Pz();
//      s = ptot[0]*ptot[0] - ptot[1]*ptot[1] - ptot[2]*ptot[2] - ptot[3]*ptot[3];
//      Qinv = -(q[0]*q[0] - q[1]*q[1] - q[2]*q[2] - q[3]*q[3]);
//      pt=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
//      Mlong=sqrt(s+pt*pt);
//      roots=sqrt(s);
//
//      Qinv1  = sqrt(Qinv);
//
//      if (pt>0)
//      {
//        Qside = (ptot[1]*q[2]-ptot[2]*q[1])/pt;
//        Qlong = (ptot[0]*q[3]-ptot[3]*q[0])/Mlong;
//        Qout  = (roots/Mlong)*(ptot[1]*q[1]+ptot[2]*q[2])/pt;
//      }
//      else
//      {
//        Qlong = q[3];
//        Qside = q[2];
//        Qout  = q[1];
//      }
//
//      h_n2_Q3D->Fill(Qlong, Qside, Qout,weight);
//    }
//
//  }

//________________________________________________________________________
void ParticlePairHistos::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    cout << "-Fatal- Attempting to load ParticleHistos from an invalid file pointer" << endl;
    return;
    }
  ParticlePairAnalyzerConfiguration & ac = *(ParticlePairAnalyzerConfiguration*) getConfiguration();
   TString bn = getHistoBaseName();
  h_n2_ptPt           = loadH2(inputFile, bn+TString("n2_ptPt"));
  h_n2_etaEta         = loadH2(inputFile, bn+TString("n2_etaEta"));
  h_n2_phiPhi         = loadH2(inputFile, bn+TString("n2_phiPhi"));
  h_n2_phiEtaPhiEta   = loadH2(inputFile, bn+TString("n2_phiEtaPhiEta"));
  h_npt_phiEtaPhiEta  = loadH2(inputFile, bn+TString("npt_phiEtaPhiEta"));
  h_ptn_phiEtaPhiEta  = loadH2(inputFile, bn+TString("ptn_phiEtaPhiEta"));
  h_ptpt_phiEtaPhiEta = loadH2(inputFile, bn+TString("ptpt_phiEtaPhiEta"));
  h_npt_etaEta        = loadH2(inputFile, bn+TString("npt_etaEta"));
  h_ptn_etaEta        = loadH2(inputFile, bn+TString("ptn_etaEta"));
  h_ptpt_etaEta       = loadH2(inputFile, bn+TString("ptpt_etaEta"));
  h_npt_phiPhi        = loadH2(inputFile, bn+TString("npt_phiPhi"));
  h_ptn_phiPhi        = loadH2(inputFile, bn+TString("ptn_phiPhi"));
  h_ptpt_phiPhi       = loadH2(inputFile, bn+TString("ptpt_phiPhi"));

  if (ac.fillY)
    {
    h_n2_yY         = loadH2(inputFile, bn+TString("n2_yY"));
    h_n2_phiYPhiY   = loadH2(inputFile, bn+TString("n2_phiYPhiY"));
    h_npt_phiYPhiY  = loadH2(inputFile, bn+TString("npt_phiYPhiY"));
    h_ptn_phiYPhiY  = loadH2(inputFile, bn+TString("ptn_phiYPhiY"));
    h_ptpt_phiYPhiY = loadH2(inputFile, bn+TString("ptpt_phiYPhiY"));
    h_npt_yY        = loadH2(inputFile, bn+TString("npt_yY"));
    h_ptn_yY        = loadH2(inputFile, bn+TString("ptn_yY"));
    h_ptpt_yY       = loadH2(inputFile, bn+TString("ptpt_yY"));
    }
  if (ac.fill6D)
    {
    h_n2_ptPhiEtaPtPhiEta = loadH2(inputFile, bn+TString("n2_ptPhiEtaPtPhiEta"));
    if (ac.fillY)
      {
      h_n2_ptPhiYPtPhiY = loadH2(inputFile, bn+TString("n2_ptPhiYPtPhiY"));
      }
    }
  return;
}
