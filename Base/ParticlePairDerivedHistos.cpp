//
//  DerivedParticleHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include "ParticlePairDerivedHistos.hpp"
ClassImp(ParticlePairDerivedHistos);

ParticlePairDerivedHistos::ParticlePairDerivedHistos(const TString & name,
                                                     ParticlePairAnalyzerConfiguration * configuration,
                                                     LogLevel  debugLevel)
 :
 Histograms(name,configuration,150,debugLevel)
 {
  // no ops
 }


 ParticlePairDerivedHistos::~ParticlePairDerivedHistos()
 {
 // no ops
 }


 //////////////////////////////////////////////////////////////////////////////////////////////////////
 // Load histograms
 //////////////////////////////////////////////////////////////////////////////////////////////////////
 void ParticlePairDerivedHistos::loadHistograms(TFile * inputFile)
 {
   if (reportDebug()) cout << "ParticlePairDerivedHistos::loadHistograms(...)" << endl;

   ParticlePairAnalyzerConfiguration & ac = * (ParticlePairAnalyzerConfiguration*) getConfiguration();
  TString bn = getHistoBaseName();
   h_n1n1_phiEtaPhiEta    = loadH2(inputFile, bn+TString("n1n1_phiEtaPhiEta"));
   h_n1n1_etaEta          = loadH2(inputFile, bn+TString("n1n1_etaEta"));
   h_n1n1_phiPhi          = loadH2(inputFile, bn+TString("n1n1_phiPhi"));
   h_n1n1_ptPt            = loadH2(inputFile, bn+TString("n1n1_ptPt"));
   h_pt1pt1_phiEtaPhiEta  = loadH2(inputFile, bn+TString("pt1pt1_phiEtaPhiEta"));
   h_pt1pt1_etaEta        = loadH2(inputFile, bn+TString("pt1pt1_etaEta"));
   h_pt1pt1_phiPhi        = loadH2(inputFile, bn+TString("pt1pt1_phiPhi"));
   h_R2_phiEtaPhiEta      = loadH2(inputFile, bn+TString("R2_phiEtaPhiEta"));
   h_R2_etaEta            = loadH2(inputFile, bn+TString("R2_etaEta"));
   h_R2_phiPhi            = loadH2(inputFile, bn+TString("R2_phiPhi"));
   h_R2_ptPt              = loadH2(inputFile, bn+TString("R2_ptPt"));
   h_DptDpt_phiEtaPhiEta  = loadH2(inputFile, bn+TString("DptDpt_phiEtaPhiEta"));
   h_DptDpt_etaEta        = loadH2(inputFile, bn+TString("DptDpt_etaEta"));
   h_DptDpt_phiPhi        = loadH2(inputFile, bn+TString("DptDpt_phiPhi"));

   h_G2_phiEtaPhiEta      = loadH2(inputFile, bn+TString("G2_phiEtaPhiEta"));
   h_G2_etaEta            = loadH2(inputFile, bn+TString("G2_etaEta"));
   h_G2_phiPhi            = loadH2(inputFile, bn+TString("G2_phiPhi"));
   h_n1n1_DetaDphi        = loadH2(inputFile, bn+TString("n1n1_DetaDphi"));
   h_pt1pt1_DetaDphi      = loadH2(inputFile, bn+TString("pt1pt1_DetaDphi"));
   h_n2_DetaDphi          = loadH2(inputFile, bn+TString("n2_DetaDphi"));
   h_ptn_DetaDphi         = loadH2(inputFile, bn+TString("ptn_DetaDphi"));
   h_npt_DetaDphi         = loadH2(inputFile, bn+TString("npt_DetaDphi"));
   h_ptpt_DetaDphi        = loadH2(inputFile, bn+TString("ptpt_DetaDphi"));
   h_R2_DetaDphi          = loadH2(inputFile, bn+TString("R2_DetaDphi"));
   h_DptDpt_DetaDphi      = loadH2(inputFile, bn+TString("DptDpt_DetaDphi"));
   h_P2_DetaDphi          = loadH2(inputFile, bn+TString("P2_DetaDphi"));
   h_G2_DetaDphi          = loadH2(inputFile, bn+TString("G2_DetaDphi"));
   h_R2_DetaDphi_shft     = loadH2(inputFile, bn+TString("R2_DetaDphi_shft"));
   h_DptDpt_DetaDphi_shft = loadH2(inputFile, bn+TString("DptDpt_DetaDphi_shft"));
   h_P2_DetaDphi_shft     = loadH2(inputFile, bn+TString("P2_DetaDphi_shft"));
   h_G2_DetaDphi_shft     = loadH2(inputFile, bn+TString("G2_DetaDphi_shft"));
   h_P2_etaEta            = loadH2(inputFile, bn+TString("P2_etaEta"));
   h_P2_phiPhi            = loadH2(inputFile, bn+TString("P2_phiPhi"));

  h_bf12_phiEtaPhiEta     = loadH2(inputFile, bn+TString("bf12_phiEtaPhiEta"));
  h_bf12_etaEta           = loadH2(inputFile, bn+TString("bf12_etaEta"));
  h_bf12_phiPhi           = loadH2(inputFile, bn+TString("bf12_phiPhi"));
  h_bf12_DetaDphi         = loadH2(inputFile, bn+TString("bf12_DetaDphi"));
  h_bf12_DetaDphi_shft    = loadH2(inputFile, bn+TString("bf12_DetaDphi_shft"));
  h_bf21_phiEtaPhiEta     = loadH2(inputFile, bn+TString("bf21_phiEtaPhiEta"));
  h_bf21_etaEta           = loadH2(inputFile, bn+TString("bf21_etaEta"));
  h_bf21_phiPhi           = loadH2(inputFile, bn+TString("bf21_phiPhi"));
  h_bf21_DetaDphi         = loadH2(inputFile, bn+TString("bf21_DetaDphi"));
  h_bf21_DetaDphi_shft    = loadH2(inputFile, bn+TString("bf21_DetaDphi_shft"));


   if (ac.fillY)
   {
     h_n1n1_phiYPhiY      = loadH2(inputFile, bn+TString("n1n1_phiYPhiY"));
     h_n1n1_yY            = loadH2(inputFile, bn+TString("n1n1_yY"));
     h_pt1pt1_phiYPhiY    = loadH2(inputFile, bn+TString(""));
     h_pt1pt1_yY          = loadH2(inputFile, bn+TString("pt1pt1_yY"));
     h_R2_phiYPhiY        = loadH2(inputFile, bn+TString("R2_phiYPhiY"));
     h_R2_yY              = loadH2(inputFile, bn+TString("R2_yY"));
     h_DptDpt_phiYPhiY    = loadH2(inputFile, bn+TString("DptDpt_phiYPhiY"));
     h_DptDpt_yY          = loadH2(inputFile, bn+TString("DptDpt_yY"));
     h_G2_phiYPhiY        = loadH2(inputFile, bn+TString("G2_phiYPhiY"));
     h_G2_yY              = loadH2(inputFile, bn+TString("G2_yY"));
     h_n1n1_DyDphi        = loadH2(inputFile, bn+TString("n1n1_DyDphi"));
     h_pt1pt1_DyDphi      = loadH2(inputFile, bn+TString("pt1pt1_DyDphi"));
     h_n2_DyDphi          = loadH2(inputFile, bn+TString("n2_DyDphi"));
     h_ptn_DyDphi         = loadH2(inputFile, bn+TString("ptn_DyDphi"));
     h_npt_DyDphi         = loadH2(inputFile, bn+TString("npt_DyDphi"));
     h_ptpt_DyDphi        = loadH2(inputFile, bn+TString("ptpt_DyDphi"));
     h_R2_DyDphi          = loadH2(inputFile, bn+TString("R2_DyDphi"));
     h_DptDpt_DyDphi      = loadH2(inputFile, bn+TString("DptDpt_DyDphi"));
     h_P2_DetaDphi        = loadH2(inputFile, bn+TString("P2_DetaDphi"));
     h_G2_DyDphi          = loadH2(inputFile, bn+TString("G2_DyDphi"));
     h_R2_DyDphi_shft     = loadH2(inputFile, bn+TString("R2_DyDphi_shft"));
     h_DptDpt_DyDphi_shft = loadH2(inputFile, bn+TString("DptDpt_DyDphi_shft"));
     h_P2_DyDphi_shft     = loadH2(inputFile, bn+TString("P2_DyDphi_shft"));
     h_G2_DyDphi_shft     = loadH2(inputFile, bn+TString("G2_DyDphi_shft"));
     h_P2_yY              = loadH2(inputFile, bn+TString("P2_yY"));

     h_bf12_phiYPhiY      = loadH2(inputFile, bn+TString("bf12_phiYPhiY"));
     h_bf12_yY            = loadH2(inputFile, bn+TString("bf12_yY"));
     h_bf12_DyDphi        = loadH2(inputFile, bn+TString("bf12_DyDphi"));
     h_bf12_DyDphi_shft   = loadH2(inputFile, bn+TString("bf12_DyDphi_shft"));
     h_bf21_phiYPhiY      = loadH2(inputFile, bn+TString("bf21_phiYPhiY"));
     h_bf21_yY            = loadH2(inputFile, bn+TString("bf21_yY"));
     h_bf21_DyDphi        = loadH2(inputFile, bn+TString("bf21_DyDphi"));
     h_bf21_DyDphi_shft   = loadH2(inputFile, bn+TString("bf21_DyDphi_shft"));
   }
   if (reportDebug()) cout << "ParticlePairDerivedHistos::loadHistograms() Completed." << endl;
 }

 //////////////////////////////////////////////////////////////////////////////////////////////////////
 // Initialize the histograms
 //////////////////////////////////////////////////////////////////////////////////////////////////////
 void ParticlePairDerivedHistos::createHistograms()
 {
   if (reportDebug()) cout << "ParticlePairDerivedHistos::initialize() starting"<< endl;
   TString name;
   ijNormalization = false;
   ParticlePairAnalyzerConfiguration & ac = * (ParticlePairAnalyzerConfiguration*) getConfiguration();
    TString bn = getHistoBaseName();
   ac.nBins_phiEta    = ac.nBins_eta*ac.nBins_phi;
   ac.nBins_Dphi_shft = ac.nBins_phi/4;
   ac.nBins_Dphi      = ac.nBins_phi;
   ac.width_Dphi      = TMath::TwoPi()/ac.nBins_Dphi;
   ac.min_Dphi        = -ac.width_Dphi/2.;
   ac.max_Dphi        = TMath::TwoPi() - ac.width_Dphi/2.;
   ac.min_Dphi_shft   = ac.min_Dphi - ac.width_Dphi*double(ac.nBins_Dphi/4);
   ac.max_Dphi_shft   = ac.max_Dphi - ac.width_Dphi*double(ac.nBins_Dphi/4);
   ac.nBins_Deta      = 2*ac.nBins_eta - 1;
   ac.width_Deta      = (ac.max_eta - ac.min_eta)/ac.nBins_Deta;
   ac.min_Deta        = ac.min_eta - ac.max_eta;
   ac.max_Deta        = ac.max_eta - ac.min_eta;

   ac.nBins_Dy        = 2*ac.nBins_y - 1;
   ac.width_Dy        = (ac.max_y - ac.min_y)/ac.nBins_Dy;
   ac.min_Dy          = ac.min_y - ac.max_y;
   ac.max_Dy          = ac.max_y - ac.min_y;


   h_n1n1_phiEtaPhiEta    = createHistogram(bn+TString("n1n1_phiEtaPhiEta"),   ac.nBins_phiEta, static_cast<double>(0.),double(ac.nBins_phiEta), ac.nBins_phiEta, 0.,double(ac.nBins_phiEta),  "#eta_{1} x #varphi_{1}", "#eta_{2} x #varphi_{2}", "<n_{1}><n_{1}>", 0,1,0,0);
   h_n1n1_etaEta          = createHistogram(bn+TString("n1n1_etaEta"),         ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,  "#eta_{1}", "#eta_{2}", "<n_{1}><n_{2}>", 0,1,0,0);
   h_n1n1_phiPhi          = createHistogram(bn+TString("n1n1_phiPhi"),         ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,  "#varphi_{1}", "#varphi_{2}", "<n_{1}><n_{2}>", 0,1,0,0);
   h_n1n1_ptPt            = createHistogram(bn+TString("n1n1_ptPt"),           ac.nBins_pt,    ac.min_pt,    ac.max_pt,    ac.nBins_pt,    ac.min_pt,   ac.max_pt,   "p_{T,1}", "p_{T,2}",  "<n_{1}><n_{2}>", 0,1,0,0);
   h_n1n1_DetaDphi        = createHistogram(bn+TString("n1n1_DetaDphi"),       ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi, ac.max_Dphi, "#Delta#eta","#Delta#varphi", "<n_{1}><n_{1}>", 0,1,0,0);

   h_pt1pt1_phiEtaPhiEta  = createHistogram(bn+TString("pt1pt1_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.),double(ac.nBins_phiEta), ac.nBins_phiEta, 0.,double(ac.nBins_phiEta),       "#eta_{1} x #varphi_{1}", "#eta_{2} x #varphi_{2}", "pt1pt1", 0,1,0,0);
   h_pt1pt1_etaEta        = createHistogram(bn+TString("pt1pt1_etaEta"),       ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,            "#eta_{1}",   "#eta_{2}",      "pt1pt1", 0,1,0,0);
   h_pt1pt1_phiPhi        = createHistogram(bn+TString("pt1pt1_phiPhi"),       ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,            "#varphi_{1}","#varphi_{2}",   "pt1pt1", 0,1,0,0);
   h_pt1pt1_DetaDphi      = createHistogram(bn+TString("pt1pt1_DetaDphi"),     ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta", "#Delta#varphi", "pt1pt1", 0,1,0,0);

   h_ptn_DetaDphi         = createHistogram(bn+TString("ptn_DetaDphi"),        ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta", "#Delta#varphi", "ptn", 0,1,0,0);
   h_npt_DetaDphi         = createHistogram(bn+TString("npt_DetaDphi"),        ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta", "#Delta#varphi", "npt", 0,1,0,0);
   h_ptpt_DetaDphi        = createHistogram(bn+TString("ptpt_DetaDphi"),       ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta", "#Delta#varphi", "ptpt", 0,1,0,0);
   h_ptpt_DetaDphi_shft   = createHistogram(bn+TString("ptpt_DetaDphi_shft"),  ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "ptpt", 0,1,0,0);

   h_n2_DetaDphi          = createHistogram(bn+TString("n2_DetaDphi"),         ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta", "#Delta#varphi", "<n_{2}>", 0,1,0,0);
   h_n2_DetaDphi_shft     = createHistogram(bn+TString("n2_DetaDphi_shft"),    ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "n_{2}", 0,1,0,0);

   h_R2_phiEtaPhiEta      = createHistogram(bn+TString("R2_phiEtaPhiEta"),     ac.nBins_phiEta, static_cast<double>(0.),double(ac.nBins_phiEta), ac.nBins_phiEta, 0.,double(ac.nBins_phiEta),       "#eta_{1} x #varphi_{1}", "#eta_{2} x #varphi_{2}", "R_{2}", 0,1,0,0);
   h_R2_etaEta            = createHistogram(bn+TString("R2_etaEta"),           ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,            "#eta_{1}",    "#eta_{2}",      "R_{2}", 0,1,0,0);
   h_R2_phiPhi            = createHistogram(bn+TString("R2_phiPhi"),           ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,            "#varphi_{1}", "#varphi_{2}",   "R_{2}", 0,1,0,0);
   h_R2_ptPt              = createHistogram(bn+TString("R2_ptPt"),             ac.nBins_pt,    ac.min_pt,    ac.max_pt,    ac.nBins_pt,    ac.min_pt,   ac.max_pt,             "p_{T,1}",     "p_{T,2}",       "R_{2}", 0,1,0,0);
   h_R2_DetaDphi          = createHistogram(bn+TString("R2_DetaDphi"),         ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta",  "#Delta#varphi", "R_{2}", 0,1,0,0);
   h_R2_DetaDphi_shft     = createHistogram(bn+TString("R2_DetaDphi_shft"),    ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "R_{2}", 0,1,0,0);

   h_DptDpt_phiEtaPhiEta  = createHistogram(bn+TString("DptDpt_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.),double(ac.nBins_phiEta), ac.nBins_phiEta, 0.,double(ac.nBins_phiEta),       "#eta_{1} x #varphi_{1}", "#eta_{2} x #varphi_{2}", "<#Delta p_{T} #Delta p_{T}>", 0,1,0,0);
   h_DptDpt_etaEta        = createHistogram(bn+TString("DptDpt_etaEta"),       ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,            "#eta_{1}",    "#eta_{2}",      "<#Delta p_{T} #Delta p_{T}>", 0,1,0,0);
   h_DptDpt_phiPhi        = createHistogram(bn+TString("DptDpt_phiPhi"),       ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,            "#varphi_{1}", "#varphi_{2}",   "<#Delta p_{T} #Delta p_{T}>", 0,1,0,0);
   h_DptDpt_DetaDphi      = createHistogram(bn+TString("DptDpt_DetaDphi"),     ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta",  "#Delta#varphi", "<#Delta p_{T}#Delta p_{T}>", 0,1,0,0);
   h_DptDpt_DetaDphi_shft = createHistogram(bn+TString("DptDpt_DetaDphi_shft"),ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "<#Delta p_{T} #Delta p_{T}>", 0,1,0,0);

   h_P2_etaEta            = createHistogram(bn+TString("P2_etaEta"),           ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,            "#eta_{1}",    "#eta_{2}",      "P_{2}", 0,1,0,0);
   h_P2_phiPhi            = createHistogram(bn+TString("P2_phiPhi"),           ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,            "#varphi_{1}", "#varphi_{2}",   "P_{2}", 0,1,0,0);
   h_P2_DetaDphi          = createHistogram(bn+TString("P2_DetaDphi"),         ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta",  "#Delta#varphi", "P_{2}", 0,1,0,0);
   h_P2_DetaDphi_shft     = createHistogram(bn+TString("P2_DetaDphi_shft"),    ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "P_{2}", 0,1,0,0);

   h_G2_phiEtaPhiEta      = createHistogram(bn+TString("G2_phiEtaPhiEta"),     ac.nBins_phiEta, static_cast<double>(0.),double(ac.nBins_phiEta), ac.nBins_phiEta, 0.,double(ac.nBins_phiEta),       "#eta x #varphi", "#eta x #varphi","G_{2}", 0,1,0,0);
   h_G2_etaEta            = createHistogram(bn+TString("G2_etaEta"),           ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,            "#eta_{1}", "#eta_{2}",         "G_{2}", 0,1,0,0);
   h_G2_phiPhi            = createHistogram(bn+TString("G2_phiPhi"),           ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,            "#varphi_{1}", "#varphi_{2}",   "G_{2}", 0,1,0,0);
   h_G2_DetaDphi          = createHistogram(bn+TString("G2_DetaDphi"),         ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta","#Delta#varphi",   "G_{2}", 0,1,0,0);
   h_G2_DetaDphi_shft     = createHistogram(bn+TString("G2_DetaDphi_shft"),    ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi",  "G_{2}", 0,1,0,0);

   h_bf12_phiEtaPhiEta    = createHistogram(bn+TString("bf12_phiEtaPhiEta"),   ac.nBins_phiEta, static_cast<double>(0.),double(ac.nBins_phiEta), ac.nBins_phiEta, 0.,double(ac.nBins_phiEta),       "#eta x #varphi", "#eta x #varphi","n_{2}^{12}/n_{1}^{2}", false, true, false, false);
   h_bf12_etaEta          = createHistogram(bn+TString("bf12_etaEta"),         ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,            "#eta_{1}", "#eta_{2}",         "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
   h_bf12_phiPhi          = createHistogram(bn+TString("bf12_phiPhi"),         ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,            "#varphi_{1}", "#varphi_{2}",   "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
   h_bf12_DetaDphi        = createHistogram(bn+TString("bf12_DetaDphi"),       ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta","#Delta#varphi",   "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
   h_bf12_DetaDphi_shft   = createHistogram(bn+TString("bf12_DetaDphi_shft"),  ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi",  "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
   h_bf21_phiEtaPhiEta    = createHistogram(bn+TString("bf21_phiEtaPhiEta"),   ac.nBins_phiEta, static_cast<double>(0.),double(ac.nBins_phiEta), ac.nBins_phiEta, 0.,double(ac.nBins_phiEta),       "#eta x #varphi", "#eta x #varphi","n_{2}^{21}/n_{1}^{1}", false, true, false, false);
   h_bf21_etaEta          = createHistogram(bn+TString("bf21_etaEta"),         ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,  ac.max_eta,            "#eta_{1}", "#eta_{2}",         "n_{2}^{21}/n_{1}^{1}", false, true, false, false);
   h_bf21_phiPhi          = createHistogram(bn+TString("bf21_phiPhi"),         ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,  ac.max_phi,            "#varphi_{1}", "#varphi_{2}",   "n_{2}^{21}/n_{1}^{1}", false, true, false, false);
   h_bf21_DetaDphi        = createHistogram(bn+TString("bf21_DetaDphi"),       ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta#eta","#Delta#varphi",   "n_{2}^{21}/n_{1}^{1}", false, true, false, false);
   h_bf21_DetaDphi_shft   = createHistogram(bn+TString("bf21_DetaDphi_shft"),  ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi",  "n_{2}^{21}/n_{1}^{1}", false, true, false, false);


   if (ac.fillY)
   {
     h_n1n1_phiYPhiY     = createHistogram(bn+TString("n1n1_phiYPhiY"),     ac.nBins_phiY, static_cast<double>(0.),double(ac.nBins_phiY), ac.nBins_phiY, 0.,double(ac.nBins_phiY),        "y_{1} x #varphi_{1}", "y_{2} x #varphi_{2}", "<n_{1}><n_{1}>", 0,1,0,0);
     h_n1n1_yY           = createHistogram(bn+TString("n1n1_yY"),           ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,   ac.min_y,  ac.max_y,                 "y_{1}",     "y_{2}", "<n_{1}><n_{1}>", 0,1,0,0);
     h_n1n1_DyDphi       = createHistogram(bn+TString("n1n1_DyDphi"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,     "#Delta y",  "#Delta#varphi", "<n_{1}><n_{1}>", 0,1,0,0);

     h_pt1pt1_phiYPhiY   = createHistogram(bn+TString("pt1pt1_phiYPhiY"),   ac.nBins_phiY, static_cast<double>(0.),double(ac.nBins_phiY), ac.nBins_phiY, 0.,double(ac.nBins_phiY),         "y_{1} x #varphi_{1}", "y_{2} x #varphi_{2}", "pt1pt1", 0,1,0,0);
     h_pt1pt1_yY         = createHistogram(bn+TString("pt1pt1_yY"),         ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,   ac.min_y,  ac.max_y,                  "y_{1}",    "y_{2}", "pt1pt1", 0,1,0,0);
     h_pt1pt1_DyDphi     = createHistogram(bn+TString("pt1pt1_DyDphi"),     ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "pt1pt1", 0,1,0,0);

     h_ptn_DyDphi        = createHistogram(bn+TString("ptn_DyDphi"),        ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "ptn", 0,1,0,0);
     h_npt_DyDphi        = createHistogram(bn+TString("npt_DyDphi"),        ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "npt", 0,1,0,0);
     h_ptpt_DyDphi       = createHistogram(bn+TString("ptpt_DyDphi"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "ptpt", 0,1,0,0);

     h_n2_DyDphi         = createHistogram(bn+TString("n2_DyDphi"),         ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "<n_{2}>", 0,1,0,0);
     h_n2_DyDphi_shft    = createHistogram(bn+TString("n2_DyDphi_shft"),    ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "n_{2}", 0,1,0,0);

     h_R2_phiYPhiY       = createHistogram(bn+TString("R2_phiYPhiY"),       ac.nBins_phiY, static_cast<double>(0.),double(ac.nBins_phiY), ac.nBins_phiY, 0.,double(ac.nBins_phiY),         "y_{1} x #varphi_{1}", "y_{2} x #varphi_{2}", "R_{2}", 0,1,0,0);
     h_R2_yY             = createHistogram(bn+TString("R2_yY"),             ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,   ac.min_y,  ac.max_y,                  "y_{1}",    "y_{2}", "R_{2}", 0,1,0,0);
     h_R2_DyDphi_shft    = createHistogram(bn+TString("R2_DyDphi_shft"),    ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "R_{2}", 0,1,0,0);
     h_R2_DyDphi         = createHistogram(bn+TString("R2_DyDphi"),         ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "R_{2}", 0,1,0,0);

     h_DptDpt_phiYPhiY   = createHistogram(bn+TString("DptDpt_phiYPhiY"),   ac.nBins_phiY, static_cast<double>(0.),double(ac.nBins_phiY), ac.nBins_phiY, 0.,double(ac.nBins_phiY),         "y_{1} x #varphi_{1}", "y_{2} x #varphi_{2}", "<#Delta p_{T} #Delta p_{T}>", 0,1,0,0);
     h_DptDpt_yY         = createHistogram(bn+TString("DptDpt_yY"),         ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,   ac.min_y,  ac.max_y,                  "y_{1}",    "y_{2}", "<#Delta p_{T} #Delta p_{T}>", 0,1,0,0);
     h_DptDpt_DyDphi     = createHistogram(bn+TString("DptDpt_DyDphi"),     ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "<#Delta p_{T}#Delta p_{T}>", 0,1,0,0);
     h_DptDpt_DyDphi_shft= createHistogram(bn+TString("DptDpt_DyDphi_shft"),ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "<#Delta p_{T}#Delta p_{T}>", 0,1,0,0);

     h_P2_yY             = createHistogram(bn+TString("P2_yY"),             ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,     ac.min_y,         ac.max_y,         "y_{1}",    "y_{2}",         "P_{2}", 0,1,0,0);
     h_P2_DyDphi         = createHistogram(bn+TString("P2_DyDphi"),         ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "P_{2}", 0,1,0,0);
     h_P2_DyDphi_shft    = createHistogram(bn+TString("P2_DyDphi_shft"),    ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "P_{2}", 0,1,0,0);

     h_G2_phiYPhiY       = createHistogram(bn+TString("G2_phiYPhiY"),       ac.nBins_phiY, static_cast<double>(0.),double(ac.nBins_phiY), ac.nBins_phiY, 0.,double(ac.nBins_phiY),         "y_{1} x #varphi_{1}", "y_{2} x #varphi_{2}", "G_{2}", 0,1,0,0);
     h_G2_yY             = createHistogram(bn+TString("G2_yY"),             ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,   ac.min_y,  ac.max_y,                  "y_{1}",    "y_{2}", "G_{2}", 0,1,0,0);
     h_G2_DyDphi         = createHistogram(bn+TString("G2_DyDphi"),         ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "G_{2}", 0,1,0,0);
     h_G2_DyDphi_shft    = createHistogram(bn+TString("G2_DyDphi_shft"),    ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "G_{2}", 0,1,0,0);


     h_bf12_phiYPhiY     = createHistogram(bn+TString("bf12_phiYPhiY"),     ac.nBins_phiY, static_cast<double>(0.),double(ac.nBins_phiY), ac.nBins_phiY, 0.,double(ac.nBins_phiY),         "y_{1} x #varphi_{1}", "y_{2} x #varphi_{2}", "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
     h_bf12_yY           = createHistogram(bn+TString("bf12_yY"),           ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,   ac.min_y,  ac.max_y,                  "y_{1}",    "y_{2}", "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
     h_bf12_DyDphi       = createHistogram(bn+TString("bf12_DyDphi"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
     h_bf12_DyDphi_shft  = createHistogram(bn+TString("bf12_DyDphi_shft"),  ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "n_{2}^{12}/n_{1}^{2}", false, true, false, false);
     h_bf21_phiYPhiY     = createHistogram(bn+TString("bf21_phiYPhiY"),     ac.nBins_phiY, static_cast<double>(0.),double(ac.nBins_phiY), ac.nBins_phiY, 0.,double(ac.nBins_phiY),         "y_{1} x #varphi_{1}", "y_{2} x #varphi_{2}", "n_{2}^{21}/n_{1}^{1}", false, true, false, false);
     h_bf21_yY           = createHistogram(bn+TString("bf21_yY"),           ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,   ac.min_y,  ac.max_y,                  "y_{1}",    "y_{2}", "n_{2}^{21}/n_{1}^{1}", false, true, false, false);
     h_bf21_DyDphi       = createHistogram(bn+TString("bf21_DyDphi"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi,      ac.max_Dphi,      "#Delta y", "#Delta#varphi", "n_{2}^{21}/n_{1}^{1}", false, true, false, false);
     h_bf21_DyDphi_shft  = createHistogram(bn+TString("bf21_DyDphi_shft"),  ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "n_{2}^{21}/n_{1}^{1}", false, true, false, false);

//     h_n1n1_Q3D          = createHistogram(bn+TString("n1n1_Q3D"),
//                                            ac.nBins_DeltaPlong, ac.min_DeltaPlong, ac.max_DeltaPlong,
//                                            ac.nBins_DeltaPside, ac.min_DeltaPside, ac.max_DeltaPside,
//                                            ac.nBins_DeltaPout,  ac.min_DeltaPout,  ac.max_DeltaPout,
//                                            "Q_{long} (GeV/c)", "Q_{side}  (GeV/c)", "Q_{out}  (GeV/c)","Yield", 0,1,0,0);
//     h_R2_Q3D             = createHistogram(bn+TString("R2_Q3D"),
//                                            ac.nBins_DeltaPlong, ac.min_DeltaPlong, ac.max_DeltaPlong,
//                                            ac.nBins_DeltaPside, ac.min_DeltaPside, ac.max_DeltaPside,
//                                            ac.nBins_DeltaPout,  ac.min_DeltaPout,  ac.max_DeltaPout,
//                                            "Q_{long} (GeV/c)", "Q_{side}  (GeV/c)", "Q_{out}  (GeV/c)","Yield", 0,1,0,0);
   }

   if (reportDebug()) cout << "ParticlePairDerivedHistos::initialize() completed."<< endl;
 }


 double ParticlePairDerivedHistos::calculateWeightCorrection(TFile * fRaw,
                                  TFile * fCal)
 {
   if (reportDebug()) cout << "ParticlePairDerivedHistos::calculateWeightCorrection() started."<< endl;

   TH1 * zV = (TH1*) fRaw->Get("zV");
   if (!zV)
   {
     cout << "calculateWeightCorrection(...) !zv.ABORT." << endl;
     return 0.;
   }

   int nZ = zV->GetNbinsX();
   double * zWeight = new double[nZ];
   double sumWeight = 0;
   double v;
   for (int iZ=0; iZ<nZ; ++iZ)
   {
     v = zV->GetBinContent(iZ+1);
     zWeight[iZ] = v;
     sumWeight += v;
   }
   for (int iZ=0; iZ<nZ; ++iZ)
   {
     zWeight[iZ] /= sumWeight;
     cout << "iZ:" << iZ << "   " << zWeight[iZ] << endl;
   }

   TH3 * correction_1;
   TH3 * correction_2;
   correction_1 = (TH3*) fCal->Get("correction_p");
   correction_2 = (TH3*) fCal->Get("correction_p");
   int index;
   int nEtaBins = 20;
   int nPhiBins = 72;
   int nEtaPhiBins = nEtaBins*nPhiBins;
   int nDetaBins  = 2*nEtaBins - 1;
   int nDphiBins  = nPhiBins;
   int nDetaDphiBins = nDetaBins*nDphiBins;
   double sumW1 = 0;
   double sumW2 = 0;
   double sumW12 = 0;
   double sumWW = 0;
   double v1, v2;
   cout << "Create arrays..." << endl;
   double *w1 = new double[nEtaPhiBins];
   double *w2 = new double[nEtaPhiBins];
   double *numerator = new double[nDetaDphiBins];
   double *denominator = new double[nDetaDphiBins];
   double *ww = new double[nEtaPhiBins*nEtaPhiBins];

   for (int iZ=0; iZ<nZ; ++iZ)
   {
     for (int iEtaPhi=0; iEtaPhi<nEtaPhiBins; ++iEtaPhi)
     {
       w1[iEtaPhi] = 0.;
       w2[iEtaPhi] = 0.;
       for (int iP=1; iP<=9; ++iP)
       {
         v1 = correction_1->GetBinContent(iZ+1,iEtaPhi+1,iP);
         v2 = correction_2->GetBinContent(iZ+1,iEtaPhi+1,iP);
         w1[iEtaPhi] += v1;
         w2[iEtaPhi] += v2;
       }
       w1[iEtaPhi] /= 9.;
       w2[iEtaPhi] /= 9.;
     }

     double sum1 = 0;
     double sum2 = 0;
     double sum12 = 0;
     for (int iEtaPhi=0; iEtaPhi<nEtaPhiBins; ++iEtaPhi)
     {
       sum1 += w1[iEtaPhi];
       sum2 += w2[iEtaPhi];
     }

     // cross product
     for (int iEtaPhi=0; iEtaPhi<nEtaPhiBins; ++iEtaPhi)
     {
       for (int jEtaPhi=0; jEtaPhi<nEtaPhiBins; ++jEtaPhi)
       {
         v1 = w1[iEtaPhi]*w2[jEtaPhi];
         ww[iEtaPhi*nEtaPhiBins + jEtaPhi] = v1;
         sum12 += v1;
       }
     }
     sum1 /= double(nEtaPhiBins);
     sum2 /= double(nEtaPhiBins);
     sum12 /= double(nEtaPhiBins*nEtaPhiBins);
     sumW1 += zWeight[iZ]*sum1;
     sumW2 += zWeight[iZ]*sum2;
     sumWW += zWeight[iZ]*sum12;

     for (int k=0;k<nDetaDphiBins;++k)
     {
       numerator[k]  = 0;
       denominator[k] = 0;
     }
     //reduce
     int j;
     int i = 0;
     for (int iEta=0;iEta<nEtaBins; ++iEta)
     {
       for (int iPhi=0;iPhi<nPhiBins; ++iPhi)
       {
         j=1;
         for (int jEta=0;jEta<nEtaBins; ++jEta)
         {
           for (int jPhi=0;jPhi<nPhiBins; ++jPhi)
           {
             int dPhi = iPhi-jPhi; if (dPhi<0) dPhi += nPhiBins; dPhi+=1;
             int dEta = iEta-jEta + nEtaBins;
             v1  = ww[i*nEtaPhiBins + j];
             index = (dEta-1)*nPhiBins + dPhi-1;
             numerator[index]  += v1;
             denominator[index] += 1.;
             //cout << " iEta:" << iEta << " iPhi:" << iPhi << " jEta:" << jEta << " jPhi:" << jPhi << " v1:" << v1 << " ev1:" << ev1 << endl;
             ++j;
           }
         }
         ++i;
       }
     }
     for (int dEta=0;dEta<nDetaBins;++dEta)
     {
       for (int dPhi=0;dPhi<nDphiBins;++dPhi)
       {
         index = dEta*nDphiBins + dPhi;
         v1  = numerator[index];
         v2  = denominator[index];
         if (v2<=0) cout << "miserable idiot!!!!!!!!!!!!!!!!!!!!!!!" << endl;
         numerator[index] = v1/v2;
       }
     }
     double sum = 0.;
     for (int i=0;i<nDetaDphiBins;++i) sum += numerator[i];
     sum /= double(nDetaDphiBins);
     sumW12 += zWeight[iZ]*sum;
     cout << "iZ:" << iZ << " " << zWeight[iZ] << " sum1:" << sum1 << " sum2:" << sum2 << " sum12:" << sum12 << " sum:" << sum << endl;

   }
   fCal->Close();

   delete[] w1;
   delete[] w2;
   delete[] numerator;
   delete[] denominator;
   delete[] ww;
   delete[] zWeight;

   cout << "==========================" << endl;
   cout << " W1:" << sumW1 << endl;
   cout << " W2:" << sumW2 << endl;
   cout << "W12:" << sumW12 << endl;
   cout << " WW:" << sumWW << endl;
   double ratio = sumWW/sumW1/sumW2;
   cout << "WW/W1/W2: " << ratio << endl;

   if (reportDebug()) cout << "ParticlePairDerivedHistos::calculateWeightCorrection() completed."<< endl;

   return 1.; //for now

 }



 // Histograms from ParticlePairDerivedHistos must be normalized "per event" before calling this function
 void ParticlePairDerivedHistos::calculateDerivedHistograms(ParticleHistos * part1Histos, ParticleHistos * part2Histos, ParticlePairHistos * pairHistos, double bincorrection)
 {
   if (reportDebug()) cout << "ParticlePairDerivedHistos::calculateDerivedHistograms() started."<< endl;
 ParticlePairAnalyzerConfiguration & ac = * (ParticlePairAnalyzerConfiguration*) getConfiguration();
   /////////////////////////////////////////////////////////////////////////////////////////////////////////
   // R2 related histograms
   /////////////////////////////////////////////////////////////////////////////////////////////////////////

   calculateN1N1_H2H2H2( part1Histos->h_n1_phiEta,   part2Histos->h_n1_phiEta,   h_n1n1_phiEtaPhiEta,   1.0, 1.0);
   calculateN1N1_H1H1H2( part1Histos->h_n1_pt,       part2Histos->h_n1_pt,       h_n1n1_ptPt,           1.0, 1.0);
   calculateN1N1_H1H1H2( part1Histos->h_n1_eta,      part2Histos->h_n1_eta,      h_n1n1_etaEta,         1.0, 1.0);
   calculateN1N1_H1H1H2( part1Histos->h_n1_phi,      part2Histos->h_n1_phi,      h_n1n1_phiPhi,         1.0, 1.0);

   calculateR2_H2H2H2(   pairHistos->h_n2_phiEtaPhiEta, h_n1n1_phiEtaPhiEta,  h_R2_phiEtaPhiEta,  ijNormalization, bincorrection, 1.0);
   calculateR2_H2H2H2(   pairHistos->h_n2_ptPt,         h_n1n1_ptPt,          h_R2_ptPt,          ijNormalization, bincorrection, 1.0);
   calculateR2_H2H2H2(   pairHistos->h_n2_etaEta,       h_n1n1_etaEta,        h_R2_etaEta,        ijNormalization, bincorrection, 1.0);
   calculateR2_H2H2H2(   pairHistos->h_n2_phiPhi,       h_n1n1_phiPhi,        h_R2_phiPhi,        ijNormalization, bincorrection, 1.0);

   reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_n2_phiEtaPhiEta,     h_n2_DetaDphi,     ac.nBins_eta, ac.nBins_phi);
   reduce_n2xEtaPhi_n2DetaDphi( h_n1n1_phiEtaPhiEta,               h_n1n1_DetaDphi,   ac.nBins_eta, ac.nBins_phi);
   reduce_n2xEtaPhi_n2DetaDphi( h_R2_phiEtaPhiEta,                 h_R2_DetaDphi,     ac.nBins_eta, ac.nBins_phi);
   symmetrizeDeltaEtaDeltaPhi(h_n2_DetaDphi,ijNormalization);
   symmetrizeDeltaEtaDeltaPhi(h_R2_DetaDphi,ijNormalization);
   shiftY(*h_n2_DetaDphi,      *h_n2_DetaDphi_shft,      ac.nBins_Dphi_shft);
   shiftY(*h_R2_DetaDphi,      *h_R2_DetaDphi_shft,      ac.nBins_Dphi_shft);

   /////////////////////////////////////////////////////////////////////////////////////////////////////////
   // Q3D related Histograms
   /////////////////////////////////////////////////////////////////////////////////////////////////////////
// if (ac.fillQ3D)
//   {
//     if (reportDebug()) cout << "ParticlePairDerivedHistos::calculateDerivedHistograms() Q3D Histograms projection."<< endl;
//
//     symmetrize3D(pairHistos->h_n2_Q3D);
//
//     h_n2_Q3D_xy = (TH2 *) pairHistos->h_n2_Q3D->Project3D("xy");     append(h_n2_Q3D_xy, 0,1,0,0);
//     h_n2_Q3D_xz = (TH2 *) pairHistos->h_n2_Q3D->Project3D("xz");     append(h_n2_Q3D_xz, 0,1,0,0);
//     h_n2_Q3D_yz = (TH2 *) pairHistos->h_n2_Q3D->Project3D("yz");     append(h_n2_Q3D_yz, 0,1,0,0);
//
//     if (ac.fillY)
//     {
//       calculateN1N1H2H2_Q3D_MCY(part1Histos->h_n1_ptY, part2Histos->h_n1_ptY, h_n1n1_Q3D,  1.0, 1.0);
//       symmetrize3D(h_n1n1_Q3D);
//     }
//     else
//     {
//       calculateN1N1H2H2_Q3D_MCEta(part1Histos->h_n1_ptEta, part2Histos->h_n1_ptEta, h_n1n1_Q3D,  1.0, 1.0);
//       symmetrize3D(h_n1n1_Q3D);
//     }
//     h_n1n1_Q3D_xy = (TH2 *) h_n1n1_Q3D->Project3D("xy"); append(h_n1n1_Q3D_xy, 0,1,0,0);
//     h_n1n1_Q3D_xz = (TH2 *) h_n1n1_Q3D->Project3D("xz"); append(h_n1n1_Q3D_xz, 0,1,0,0);
//     h_n1n1_Q3D_yz = (TH2 *) h_n1n1_Q3D->Project3D("yz"); append(h_n1n1_Q3D_yz, 0,1,0,0);
//
//     calculateR2_Q3D(pairHistos->h_n2_Q3D,h_n1n1_Q3D,h_R2_Q3D, 1.0, 1.0);
//     h_R2_Q3D_xy = (TH2 *) h_R2_Q3D->Project3D("xy");     append(h_R2_Q3D_xy,  0,1,0,0);
//     h_R2_Q3D_xz = (TH2 *) h_R2_Q3D->Project3D("xz");     append(h_R2_Q3D_xz,  0,1,0,0);
//     h_R2_Q3D_yz = (TH2 *) h_R2_Q3D->Project3D("yz");     append(h_R2_Q3D_yz,  0,1,0,0);
//     //h_R2_Q3D_x  = h_R2_Q3D->Project3D("x");              append(h_R2_Q3D_x,   notScaled, notSaved, plotted, notPrinted);
//     //h_R2_Q3D_y  = h_R2_Q3D->Project3D("y");              append(h_R2_Q3D_y,   notScaled, notSaved, plotted, notPrinted);
//     //h_R2_Q3D_z  = h_R2_Q3D->Project3D("z");              append(h_R2_Q3D_z,   notScaled, notSaved, plotted, notPrinted);
//
//     if (reportDebug()) cout << "ParticlePairDerivedHistos::calculateDerivedHistograms() Q3D Histograms projection completed!!!!!"<< endl;
//   }


   /////////////////////////////////////////////////////////////////////////////////////////////////////////
   // DptDpt related Histograms
   /////////////////////////////////////////////////////////////////////////////////////////////////////////

   calculateN1N1_H1H1H2( part1Histos->h_pt_phi,      part2Histos->h_pt_phi,      h_pt1pt1_phiPhi,       1.0, 1.0);
   calculateN1N1_H1H1H2( part1Histos->h_pt_eta,      part2Histos->h_pt_eta,      h_pt1pt1_etaEta,       1.0, 1.0);
   calculateN1N1_H2H2H2( part1Histos->h_pt_phiEta,   part2Histos->h_pt_phiEta,   h_pt1pt1_phiEtaPhiEta, 1.0, 1.0);


   calculateDptDpt(pairHistos->h_ptpt_phiEtaPhiEta, pairHistos->h_ptn_phiEtaPhiEta, pairHistos->h_npt_phiEtaPhiEta,  pairHistos->h_n2_phiEtaPhiEta,
                   part1Histos->h_pt_phiEta,
                   part2Histos->h_pt_phiEta,
                   h_DptDpt_phiEtaPhiEta,h_DptDpt_phiEtaPhiEta,
                   ijNormalization, ac.nBins_eta, ac.nBins_phi);

   calculateDptDpt(pairHistos->h_ptpt_etaEta, pairHistos->h_ptn_etaEta, pairHistos->h_npt_etaEta,  pairHistos->h_n2_etaEta,
                   part1Histos->h_pt_eta,
                   part2Histos->h_pt_eta,
                   h_DptDpt_etaEta,
                   ijNormalization, ac.nBins_eta);

   calculateDptDpt(pairHistos->h_ptpt_phiPhi, pairHistos->h_ptn_phiPhi, pairHistos->h_npt_phiPhi,  pairHistos->h_n2_phiPhi,
                   part1Histos->h_pt_phi,
                   part2Histos->h_pt_phi,
                   h_DptDpt_phiPhi,
                   ijNormalization, ac.nBins_phi);

   reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_ptpt_phiEtaPhiEta,   h_ptpt_DetaDphi,   ac.nBins_eta, ac.nBins_phi);
   reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_ptn_phiEtaPhiEta,    h_ptn_DetaDphi,    ac.nBins_eta, ac.nBins_phi);
   reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_npt_phiEtaPhiEta,    h_npt_DetaDphi,    ac.nBins_eta, ac.nBins_phi);
   reduce_n2xEtaPhi_n2DetaDphi( h_pt1pt1_phiEtaPhiEta,             h_pt1pt1_DetaDphi, ac.nBins_eta, ac.nBins_phi);
   reduce_n2xEtaPhi_n2DetaDphi( h_DptDpt_phiEtaPhiEta,             h_DptDpt_DetaDphi, ac.nBins_eta, ac.nBins_phi);

   /// calculate average pt base on h_pt_eta histograms
   double avgPt1Eta = avgValue(part1Histos->h_pt_eta);
   double avgPt2Eta = avgValue(part2Histos->h_pt_eta);
   if (reportInfo()) cout << "calculateDerivedHistograms(...) avgPt1Eta = " << avgPt1Eta << endl;
   if (reportInfo()) cout << "calculateDerivedHistograms(...) avgPt2Eta = " << avgPt2Eta << endl;
   h_P2_DetaDphi->Add(h_DptDpt_DetaDphi,h_DptDpt_DetaDphi,1.0,0.0);
   h_P2_DetaDphi->Scale(1.0/avgPt1Eta/avgPt2Eta);

   symmetrizeDeltaEtaDeltaPhi(h_DptDpt_DetaDphi,ijNormalization);
   symmetrizeDeltaEtaDeltaPhi(h_P2_DetaDphi,ijNormalization);
   shiftY(*h_DptDpt_DetaDphi,  *h_DptDpt_DetaDphi_shft,  ac.nBins_Dphi_shft);
   shiftY(*h_P2_DetaDphi,      *h_P2_DetaDphi_shft,      ac.nBins_Dphi_shft);
   shiftY(*h_ptpt_DetaDphi,     *h_ptpt_DetaDphi_shft,    ac.nBins_Dphi_shft);


 /////////////////////////////////////////////////////////////////////////////////////////////////////////
 // G2 related histograms
 /////////////////////////////////////////////////////////////////////////////////////////////////////////
 // will go here...

   calculateG2_H2H2H2H2( pairHistos->h_ptpt_phiEtaPhiEta, h_n1n1_phiEtaPhiEta,  h_pt1pt1_phiEtaPhiEta, h_G2_phiEtaPhiEta, ijNormalization, bincorrection, 1.0);
   calculateG2_H2H2H2H2( pairHistos->h_ptpt_etaEta, h_n1n1_etaEta,  h_pt1pt1_etaEta, h_G2_etaEta, ijNormalization, bincorrection, 1.0);
   calculateG2_H2H2H2H2( pairHistos->h_ptpt_phiPhi, h_n1n1_phiPhi,  h_pt1pt1_phiPhi, h_G2_phiPhi, ijNormalization, bincorrection, 1.0);
   reduce_n2xEtaPhi_n2DetaDphi( h_G2_phiEtaPhiEta,                 h_G2_DetaDphi,     ac.nBins_eta, ac.nBins_phi);
   symmetrizeDeltaEtaDeltaPhi(h_G2_DetaDphi,ijNormalization);
   shiftY(*h_G2_DetaDphi,      *h_G2_DetaDphi_shft,      ac.nBins_Dphi_shft);

   /*
    symmetrizeDeltaEtaDeltaPhi( h_R2_DetaDphi,     0);
    symmetrizeDeltaEtaDeltaPhi( h_ptpt_DetaDphi,   ijNormalization);
    symmetrizeDeltaEtaDeltaPhi( h_ptn_DetaDphi,    ijNormalization);
    symmetrizeDeltaEtaDeltaPhi( h_npt_DetaDphi,    ijNormalization);
    symmetrizeDeltaEtaDeltaPhi( h_DptDpt_DetaDphi, ijNormalization);
    symmetrizeDeltaEtaDeltaPhi( h_ptpt_DetaDphi,   0);
    symmetrizeDeltaEtaDeltaPhi( h_pt1pt1_DetaDphi, 0);
    symmetrizeDeltaEtaDeltaPhi( h_DptDpt_DetaDphi, 0);
    symmetrizeDeltaEtaDeltaPhi( h_G2_DetaDphi,     0);
    symmetrizeXX( h_R2_ptPt,     ijNormalization);
    symmetrizeXX( h_R2_etaEta,   ijNormalization);
    symmetrizeXX( h_R2_phiPhi,   ijNormalization);
    */

   /* balance functions components not R2 based */
   calculateBf( pairHistos->h_n2_phiEtaPhiEta,part1Histos->h_n1_phiEta,part2Histos->h_n1_phiEta,h_bf12_phiEtaPhiEta,h_bf21_phiEtaPhiEta );
   calculateBf( pairHistos->h_n2_etaEta,part1Histos->h_n1_phiEta,part2Histos->h_n1_phiEta,h_bf12_etaEta,h_bf21_etaEta );
   calculateBf( pairHistos->h_n2_phiPhi,part1Histos->h_n1_phiEta,part2Histos->h_n1_phiEta,h_bf12_phiPhi,h_bf21_phiPhi );
   reduce_n2xEtaPhi_n2DetaDphi( h_bf12_phiEtaPhiEta, h_bf12_DetaDphi, ac.nBins_eta, ac.nBins_phi);
   reduce_n2xEtaPhi_n2DetaDphi( h_bf21_phiEtaPhiEta, h_bf21_DetaDphi, ac.nBins_eta, ac.nBins_phi);
   shiftY(*h_bf12_DetaDphi, *h_bf12_DetaDphi_shft, ac.nBins_Dphi_shft);
   shiftY(*h_bf21_DetaDphi, *h_bf21_DetaDphi_shft, ac.nBins_Dphi_shft);

   if (ac.fillY)
   {
     calculateN1N1_H2H2H2( part1Histos->h_n1_phiY,     part2Histos->h_n1_phiY, h_n1n1_phiYPhiY,   1.0, 1.0);
     calculateN1N1_H1H1H2( part1Histos->h_n1_y,        part2Histos->h_n1_y,    h_n1n1_yY,         1.0, 1.0);
     calculateR2_H2H2H2(   pairHistos->h_n2_phiYPhiY,  h_n1n1_phiYPhiY,        h_R2_phiYPhiY,     ijNormalization, bincorrection, 1.0);
     calculateR2_H2H2H2(   pairHistos->h_n2_yY,        h_n1n1_yY,              h_R2_yY,           ijNormalization, bincorrection, 1.0);

     calculateN1N1_H1H1H2( part1Histos->h_pt_y,    part2Histos->h_pt_y,        h_pt1pt1_yY,       1.0, 1.0);
     calculateN1N1_H2H2H2( part1Histos->h_pt_phiY, part2Histos->h_pt_phiY,     h_pt1pt1_phiYPhiY, 1.0, 1.0);
     calculateDptDpt(pairHistos->h_ptpt_phiYPhiY,  pairHistos->h_ptn_phiYPhiY, pairHistos->h_npt_phiYPhiY,  pairHistos->h_n2_phiYPhiY,
                     part1Histos->h_pt_phiY,
                     part2Histos->h_pt_phiY,
                     h_DptDpt_phiYPhiY,h_DptDpt_phiYPhiY,
                     ijNormalization, ac.nBins_y, ac.nBins_phi);
     calculateDptDpt(pairHistos->h_ptpt_yY, pairHistos->h_ptn_yY, pairHistos->h_npt_yY,  pairHistos->h_n2_yY,
                     part1Histos->h_pt_y,
                     part2Histos->h_pt_y,
                     h_DptDpt_yY,
                     ijNormalization, ac.nBins_y);
     calculateG2_H2H2H2H2( pairHistos->h_ptpt_phiYPhiY, h_n1n1_phiYPhiY,  h_pt1pt1_phiYPhiY, h_G2_phiYPhiY, ijNormalization, bincorrection, 1.0);
     calculateG2_H2H2H2H2( pairHistos->h_ptpt_yY, h_n1n1_yY,  h_pt1pt1_yY, h_G2_yY, ijNormalization, bincorrection, 1.0);
     calculateG2_H2H2H2H2( pairHistos->h_ptpt_phiPhi, h_n1n1_phiPhi,  h_pt1pt1_phiPhi, h_G2_phiPhi, ijNormalization, bincorrection, 1.0);

     reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_n2_phiYPhiY,   h_n2_DyDphi,     ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_ptpt_phiYPhiY, h_ptpt_DyDphi,   ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_ptn_phiYPhiY,  h_ptn_DyDphi,    ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( pairHistos->h_npt_phiYPhiY,  h_npt_DyDphi,    ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( h_n1n1_phiYPhiY,             h_n1n1_DyDphi,   ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( h_R2_phiYPhiY,               h_R2_DyDphi,     ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( h_pt1pt1_phiYPhiY,           h_pt1pt1_DyDphi, ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( h_DptDpt_phiYPhiY,           h_DptDpt_DyDphi, ac.nBins_y, ac.nBins_phi);
     reduce_n2xEtaPhi_n2DetaDphi( h_G2_phiYPhiY,               h_G2_DyDphi,     ac.nBins_y, ac.nBins_phi);

     /// calculate average pt base on h_pt_eta histograms
     double avgPt1Y = avgValue(part1Histos->h_pt_y);
     double avgPt2Y = avgValue(part2Histos->h_pt_y);
     if (reportInfo()) cout << "calculateDerivedHistograms(...) avgPt1Y = " << avgPt1Y << endl;
     if (reportInfo()) cout << "calculateDerivedHistograms(...) avgPt2Y = " << avgPt2Y << endl;
     h_P2_DyDphi->Add(h_DptDpt_DyDphi,h_DptDpt_DyDphi,1.0,0.0);
     h_P2_DyDphi->Scale(1.0/avgPt1Y/avgPt2Y);

     symmetrizeDeltaEtaDeltaPhi(h_n2_DyDphi,ijNormalization);
     symmetrizeDeltaEtaDeltaPhi(h_R2_DyDphi,ijNormalization);
     symmetrizeDeltaEtaDeltaPhi(h_DptDpt_DyDphi,ijNormalization);
     symmetrizeDeltaEtaDeltaPhi(h_P2_DyDphi,ijNormalization);
     symmetrizeDeltaEtaDeltaPhi(h_G2_DyDphi,ijNormalization);

     shiftY(*h_n2_DyDphi,      *h_n2_DyDphi_shft,      ac.nBins_Dphi_shft);
     shiftY(*h_R2_DyDphi,      *h_R2_DyDphi_shft,      ac.nBins_Dphi_shft);
     shiftY(*h_DptDpt_DyDphi,  *h_DptDpt_DyDphi_shft,  ac.nBins_Dphi_shft);
     shiftY(*h_P2_DyDphi,      *h_P2_DyDphi_shft,      ac.nBins_Dphi_shft);
     shiftY(*h_G2_DyDphi,      *h_G2_DyDphi_shft,      ac.nBins_Dphi_shft);

    /* balance functions components not R2 based */
    calculateBf( pairHistos->h_n2_phiYPhiY,part1Histos->h_n1_phiY,part2Histos->h_n1_phiY,h_bf12_phiYPhiY,h_bf21_phiYPhiY );
    calculateBf( pairHistos->h_n2_yY,part1Histos->h_n1_phiY,part2Histos->h_n1_phiY,h_bf12_yY,h_bf21_yY );
    reduce_n2xEtaPhi_n2DetaDphi( h_bf12_phiYPhiY, h_bf12_DyDphi, ac.nBins_y, ac.nBins_phi);
    reduce_n2xEtaPhi_n2DetaDphi( h_bf21_phiYPhiY, h_bf21_DyDphi, ac.nBins_y, ac.nBins_phi);
    shiftY(*h_bf12_DyDphi, *h_bf12_DyDphi_shft, ac.nBins_Dphi_shft);
    shiftY(*h_bf21_DyDphi, *h_bf21_DyDphi_shft, ac.nBins_Dphi_shft);
   }

   if (reportDebug()) cout << "ParticlePairDerivedHistos::calculateDerivedHistograms() completed."<< endl;
 }
