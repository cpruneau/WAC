//
//  ParticlePairCombinedHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 12/18/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "ParticlePairCombinedHistos.hpp"
ClassImp(ParticlePairCombinedHistos)


ParticlePairCombinedHistos::ParticlePairCombinedHistos(const TString & name,
                                                       AnalysisConfiguration * configuration,
                                                       LogLevel  debugLevel)
  :
  Histograms(name,configuration,100,debugLevel)
  {
    initialize();
  }


  ParticlePairCombinedHistos::ParticlePairCombinedHistos(TFile * inputFile,
                             const TString & collectionName,
                             AnalysisConfiguration * analysisConfiguration,
                             LogLevel  debugLevel)
  :
  Histograms(collectionName,analysisConfiguration,100,debugLevel)
  {
     loadHistograms(inputFile);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // DTOR
  //////////////////////////////////////////////////////////////////////////////////////////////////////
ParticlePairCombinedHistos::~ParticlePairCombinedHistos()
  {
    if (reportDebug()) cout << "ParticlePairCombinedHistos::DTOR() No ops." << endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // load histograms
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticlePairCombinedHistos::loadHistograms(TFile * inputFile)
  {
    if (reportDebug()) cout << "ParticlePairCombinedHistos::loadHistograms(...) strted" << endl;
    AnalysisConfiguration & ac = *getConfiguration();
    TString bn = getHistoBaseName();
    h_R2_etaEta        = loadH2(inputFile,bn+TString("R2_etaEta"));
    h_R2_phiPhi        = loadH2(inputFile,bn+TString("R2_phiPhi"));
    h_R2_ptPt          = loadH2(inputFile,bn+TString("R2_ptPt"));
    h_R2_DetaDphi_shft = loadH2(inputFile,bn+TString("R2_DetaDphi_shft"));

    h_DptDpt_etaEta        = loadH2(inputFile,bn+TString("DptDpt_etaEta"));
    h_DptDpt_phiPhi        = loadH2(inputFile,bn+TString("DptDpt_phiPhi"));
    h_DptDpt_DetaDphi_shft = loadH2(inputFile,bn+TString("DptDpt_DetaDphi_shft"));

    h_P2_DetaDphi_shft     = loadH2(inputFile,bn+TString("P2_DetaDphi_shft"));

    h_G2_etaEta            = loadH2(inputFile,bn+TString("G2_etaEta"));
    h_G2_phiPhi            = loadH2(inputFile,bn+TString("G2_phiPhi"));
    h_G2_DetaDphi_shft     = loadH2(inputFile,bn+TString("G2_DetaDphi_shft"));

    if (ac.fillY)
    {
      h_R2_yY              = loadH2(inputFile,bn+TString("R2_yY"));
      h_R2_DyDphi_shft     = loadH2(inputFile,bn+TString("R2_DyDphi_shft"));
      h_DptDpt_yY          = loadH2(inputFile,bn+TString("DptDpt_yY"));
      h_DptDpt_DyDphi_shft = loadH2(inputFile,bn+TString("DptDpt_DyDphi_shft"));
      h_P2_DyDphi_shft     = loadH2(inputFile,bn+TString("P2_DyDphi_shft"));
      h_G2_DyDphi_shft     = loadH2(inputFile,bn+TString("G2_DyDphi_shft"));
    }
    if (reportDebug()) cout << "ParticlePairCombinedHistos::loadHistograms(...) completed" << endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // initialize histograms of this class
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticlePairCombinedHistos::createHistograms()
  {
    if (reportDebug()) cout << "ParticlePairCombinedHistos::initialize() started" << endl;
    TString name;

    AnalysisConfiguration & ac = *getConfiguration();
    TString bn = getHistoBaseName();
    ac.nBins_phiEta     = ac.nBins_eta*ac.nBins_phi;
    ac.nBins_Dphi_shft  = ac.nBins_phi/4;
    ac.nBins_Dphi       = ac.nBins_phi;
    ac.width_Dphi       = TMath::TwoPi()/ac.nBins_Dphi;
    ac.min_Dphi         = -ac.width_Dphi/2.;
    ac.max_Dphi         = TMath::TwoPi() - ac.width_Dphi/2.;
    ac.min_Dphi_shft    = ac.min_Dphi - ac.width_Dphi*double(ac.nBins_Dphi/4);
    ac.max_Dphi_shft    = ac.max_Dphi - ac.width_Dphi*double(ac.nBins_Dphi/4);
    ac.nBins_Deta       = 2*ac.nBins_eta - 1;
    ac.width_Deta       = (ac.max_eta - ac.min_eta)/ac.nBins_Deta;
    ac.min_Deta         = ac.min_eta - ac.max_eta;
    ac.max_Deta         = ac.max_eta - ac.min_eta;

    ac.nBins_Dy         = 2*ac.nBins_y - 1;
    ac.width_Dy         = (ac.max_y - ac.min_y)/ac.nBins_Dy;
    ac.min_Dy           = ac.min_y - ac.max_y;
    ac.max_Dy           = ac.max_y - ac.min_y;


    h_R2_etaEta            = createHistogram(bn+TString("R2_etaEta"),            ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,       ac.max_eta,       "#eta_{1}",    "#eta_{2}",      "R_{2}", true,true,true,true);
    h_R2_phiPhi            = createHistogram(bn+TString("R2_phiPhi"),            ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,       ac.max_phi,       "#varphi_{1}", "#varphi_{2}",   "R_{2}", true,true,true,true);
    h_R2_ptPt              = createHistogram(bn+TString("R2_ptPt"),              ac.nBins_pt,    ac.min_pt,    ac.max_pt,    ac.nBins_pt,    ac.min_pt,        ac.max_pt,        "p_{T,1}",     "p_{T,2}",       "R_{2}", true,true,true,true);
    h_R2_DetaDphi_shft     = createHistogram(bn+TString("R2_DetaDphi_shft"),     ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "R_{2}", true,true,true,true);

    h_DptDpt_etaEta        = createHistogram(bn+TString("DptDpt_etaEta"),        ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,       ac.max_eta,       "#eta_{1}",    "#eta_{2}",      "<#Delta p_{T} #Delta p_{T}>", true,true,true,true);
    h_DptDpt_phiPhi        = createHistogram(bn+TString("DptDpt_phiPhi"),        ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,       ac.max_phi,       "#varphi_{1}", "#varphi_{2}",   "<#Delta p_{T} #Delta p_{T}>", true,true,true,true);
    h_DptDpt_DetaDphi_shft = createHistogram(bn+TString("DptDpt_DetaDphi_shft"), ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "<#Delta p_{T} #Delta p_{T}>", true,true,true,true);

    h_P2_etaEta            = createHistogram(bn+TString("P2_etaEta"),            ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,       ac.max_eta,       "#eta_{1}",    "#eta_{2}",      "P_{2}", true,true,true,true);
    h_P2_phiPhi            = createHistogram(bn+TString("P2_phiPhi"),            ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,       ac.max_phi,       "#varphi_{1}", "#varphi_{2}",   "P_{2}", true,true,true,true);
    h_P2_DetaDphi_shft     = createHistogram(bn+TString("P2_DetaDphi_shft"),     ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "P_{2}", true,true,true,true);

    h_G2_etaEta            = createHistogram(bn+TString("G2_etaEta"),            ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,       ac.max_eta,       "#eta_{1}",    "#eta_{2}",      "G_{2}", true,true,true,true);
    h_G2_phiPhi            = createHistogram(bn+TString("G2_phiPhi"),            ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,       ac.max_phi,       "#varphi_{1}", "#varphi_{2}",   "G_{2}", true,true,true,true);
    h_G2_DetaDphi_shft     = createHistogram(bn+TString("G2_DetaDphi_shft"),     ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "G_{2}", true,true,true,true);

    h_BF_etaEta            = createHistogram(bn+TString("BF_etaEta"),            ac.nBins_eta,   ac.min_eta,   ac.max_eta,   ac.nBins_eta,   ac.min_eta,       ac.max_eta,       "#eta_{1}",    "#eta_{2}",      "BF_{2}", true,true,true,true);
    h_BF_phiPhi            = createHistogram(bn+TString("BF_phiPhi"),            ac.nBins_phi,   ac.min_phi,   ac.max_phi,   ac.nBins_phi,   ac.min_phi,       ac.max_phi,       "#varphi_{1}", "#varphi_{2}",   "BF_{2}", true,true,true,true);
    h_BF_DetaDphi_shft     = createHistogram(bn+TString("BF_DetaDphi_shft"),     ac.nBins_Deta,  ac.min_Deta,  ac.max_Deta,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta",  "#Delta#varphi", "BF_{2}", true,true,true,true);


    if (ac.fillY)
    {
      h_R2_yY              = createHistogram(bn+TString("R2_yY"),                ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,     ac.min_y,         ac.max_y,         "y_{1}",      "y_{2}", "R_{2}", true,true,true,true);
      h_R2_DyDphi_shft     = createHistogram(bn+TString("R2_DyDphi_shft"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y",   "#Delta#varphi", "R_{2}", true,true,true,true);

      h_DptDpt_yY          = createHistogram(bn+TString("DptDpt_yY"),            ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,     ac.min_y,         ac.max_y,         "y_{1}",      "y_{2}",         "<#Delta p_{T} #Delta p_{T}>", true,true,true,true);
      h_DptDpt_DyDphi_shft = createHistogram(bn+TString("DptDpt_DyDphi_shft"),   ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y",   "#Delta#varphi", "<#Delta p_{T} #Delta p_{T}>", true,true,true,true);

      h_P2_yY              = createHistogram(bn+TString("P2_yY"),                ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,     ac.min_y,         ac.max_y,         "y_{1}",      "y_{2}",         "P_{2}", true,true,true,true);
      h_P2_DyDphi_shft     = createHistogram(bn+TString("P2_DyDphi_shft"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y",   "#Delta#varphi", "P_{2}", true,true,true,true);

      h_G2_yY              = createHistogram(bn+TString("G2_yY"),                ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,     ac.min_y,         ac.max_y,         "y_{1}",      "y_{2}",         "G_{2}", true,true,true,true);
      h_G2_DyDphi_shft     = createHistogram(bn+TString("G2_DyDphi_shft"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y",   "#Delta#varphi", "G_{2}", true,true,true,true);

      h_BF_yY              = createHistogram(bn+TString("BF_yY"),                ac.nBins_y,   ac.min_y,   ac.max_y,   ac.nBins_y,     ac.min_y,         ac.max_y,         "y_{1}",      "y_{2}",         "BF_{2}", true,true,true,true);
      h_BF_DyDphi_shft     = createHistogram(bn+TString("BF_DyDphi_shft"),       ac.nBins_Dy,  ac.min_Dy,  ac.max_Dy,  ac.nBins_Dphi,  ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y",   "#Delta#varphi", "BF_{2}", true,true,true,true);
    }

  if (ac.fillQ3D)
    {
      h_R2_Q3D = createHistogram(bn+TString("R2_Q3D"),
                                 ac.nBins_DeltaPlong, ac.min_DeltaPlong, ac.max_DeltaPlong,
                                 ac.nBins_DeltaPside, ac.min_DeltaPside, ac.max_DeltaPside,
                                 ac.nBins_DeltaPout,  ac.min_DeltaPout,  ac.max_DeltaPout,
                                 "Q_{long} (GeV/c)", "Q_{side}  (GeV/c)", "Q_{out}  (GeV/c)","R_{2}", true,true,false,false);

      h_RR2_Q3D = createHistogram(bn+TString("RR2_Q3D"),
                                  ac.nBins_DeltaPlong, ac.min_DeltaPlong, ac.max_DeltaPlong,
                                  ac.nBins_DeltaPside, ac.min_DeltaPside, ac.max_DeltaPside,
                                  ac.nBins_DeltaPout,  ac.min_DeltaPout,  ac.max_DeltaPout,
                                  "Q_{long} (GeV/c)", "Q_{side}  (GeV/c)", "Q_{out}  (GeV/c)","RR_{2}", true,true,false,false);

    }


    if (reportDebug()) cout << "ParticlePairCombinedHistos::initialize() completed" << endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // calculate combinations of histograms with given coefficients
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  void ParticlePairCombinedHistos::calculate(ParticlePairDerivedHistos * pp, ParticlePairDerivedHistos * mm, ParticlePairDerivedHistos * pm, double app, double amm, double apm)
  {
    if (reportDebug()) cout << "ParticlePairCombinedHistos::calculate() started" << endl;
    AnalysisConfiguration & ac = *getConfiguration();

    h_R2_etaEta->Add(pp->h_R2_etaEta, mm->h_R2_etaEta, app, amm);
    h_R2_etaEta->Add(pm->h_R2_etaEta, apm);

    h_R2_phiPhi->Add(pp->h_R2_phiPhi, mm->h_R2_phiPhi, app, amm);
    h_R2_phiPhi->Add(pm->h_R2_phiPhi, apm);

    h_R2_ptPt->Add(pp->h_R2_ptPt, mm->h_R2_ptPt, app, amm);
    h_R2_ptPt->Add(pm->h_R2_ptPt, apm);

    h_DptDpt_etaEta->Add(pp->h_DptDpt_etaEta, mm->h_DptDpt_etaEta, app, amm);
    h_DptDpt_etaEta->Add(pm->h_DptDpt_etaEta, apm);

    h_DptDpt_phiPhi->Add(pp->h_DptDpt_phiPhi, mm->h_DptDpt_phiPhi, app, amm);
    h_DptDpt_phiPhi->Add(pm->h_DptDpt_phiPhi, apm);

    h_G2_etaEta->Add(pp->h_G2_etaEta, mm->h_G2_etaEta, app, amm);
    h_G2_etaEta->Add(pm->h_G2_etaEta, apm);

    h_G2_phiPhi->Add(pp->h_G2_phiPhi, mm->h_G2_phiPhi, app, amm);
    h_G2_phiPhi->Add(pm->h_G2_phiPhi, apm);

    h_R2_DetaDphi_shft->Add(pp->h_R2_DetaDphi_shft, mm->h_R2_DetaDphi_shft, app, amm);
    h_R2_DetaDphi_shft->Add(pm->h_R2_DetaDphi_shft, apm);

    h_DptDpt_DetaDphi_shft->Add(pp->h_DptDpt_DetaDphi_shft, mm->h_DptDpt_DetaDphi_shft, app, amm);
    h_DptDpt_DetaDphi_shft->Add(pm->h_DptDpt_DetaDphi_shft, apm);

    h_P2_DetaDphi_shft->Add(pp->h_P2_DetaDphi_shft, mm->h_P2_DetaDphi_shft, app, amm);
    h_P2_DetaDphi_shft->Add(pm->h_P2_DetaDphi_shft, apm);

    h_G2_DetaDphi_shft->Add(pp->h_G2_DetaDphi_shft, mm->h_G2_DetaDphi_shft, app, amm);
    h_G2_DetaDphi_shft->Add(pm->h_G2_DetaDphi_shft, apm);

    h_BF_DetaDphi_shft->Add(pp->h_bf12_DetaDphi_shft, mm->h_bf12_DetaDphi_shft, app, amm);
    h_BF_DetaDphi_shft->Add(pm->h_bf12_DetaDphi_shft, apm/2.0);
    h_BF_DetaDphi_shft->Add(pm->h_bf21_DetaDphi_shft, apm/2.0);

    h_BF_etaEta->Add(pp->h_bf12_etaEta, mm->h_bf12_etaEta, app, amm);
    h_BF_etaEta->Add(pm->h_bf12_etaEta, apm/2.0);
    h_BF_etaEta->Add(pm->h_bf21_etaEta, apm/2.0);

    h_BF_phiPhi->Add(pp->h_bf12_phiPhi, mm->h_bf12_phiPhi, app, amm);
    h_BF_phiPhi->Add(pm->h_bf12_phiPhi, apm/2.0);
    h_BF_phiPhi->Add(pm->h_bf21_phiPhi, apm/2.0);

    if (ac.fillY)
    {
      h_R2_yY->Add(pp->h_R2_yY, mm->h_R2_yY, app, amm);
      h_R2_yY->Add(pm->h_R2_yY, apm);

      h_DptDpt_yY->Add(pp->h_DptDpt_yY, mm->h_DptDpt_yY, app, amm);
      h_DptDpt_yY->Add(pm->h_DptDpt_yY, apm);

      h_G2_yY->Add(pp->h_G2_yY, mm->h_G2_yY, app, amm);
      h_G2_yY->Add(pm->h_G2_yY, apm);

      h_R2_DyDphi_shft->Add(pp->h_R2_DyDphi_shft, mm->h_R2_DyDphi_shft, app, amm);
      h_R2_DyDphi_shft->Add(pm->h_R2_DyDphi_shft, apm);

      h_DptDpt_DyDphi_shft->Add(pp->h_DptDpt_DyDphi_shft, mm->h_DptDpt_DyDphi_shft, app, amm);
      h_DptDpt_DyDphi_shft->Add(pm->h_DptDpt_DyDphi_shft, apm);

      h_P2_DyDphi_shft->Add(pp->h_P2_DyDphi_shft, mm->h_P2_DyDphi_shft, app, amm);
      h_P2_DyDphi_shft->Add(pm->h_P2_DyDphi_shft, apm);

      h_G2_DyDphi_shft->Add(pp->h_G2_DyDphi_shft, mm->h_G2_DyDphi_shft, app, amm);
      h_G2_DyDphi_shft->Add(pm->h_G2_DyDphi_shft, apm);

      h_BF_yY->Add(pp->h_bf12_yY, mm->h_bf12_yY, app, amm);
      h_BF_yY->Add(pm->h_bf12_yY, apm/2.0);
      h_BF_yY->Add(pm->h_bf21_yY, apm/2.0);

      h_BF_DyDphi_shft->Add(pp->h_bf12_DyDphi_shft, mm->h_bf12_DyDphi_shft, app, amm);
      h_BF_DyDphi_shft->Add(pm->h_bf12_DyDphi_shft, apm/2.0);
      h_BF_DyDphi_shft->Add(pm->h_bf21_DyDphi_shft, apm/2.0);
    }

    if (ac.fillQ3D)
    {
      h_R2_Q3D->Add(pp->h_R2_Q3D, mm->h_R2_Q3D, app, amm);
      h_R2_Q3D->Add(pm->h_R2_Q3D, apm);
      h_R2_Q3D_xy = (TH2 *) h_R2_Q3D->Project3D("xy");     add(h_R2_Q3D_xy, false,true,true,true);
      h_R2_Q3D_xz = (TH2 *) h_R2_Q3D->Project3D("xz");     add(h_R2_Q3D_xz, false,true,true,true);
      h_R2_Q3D_yz = (TH2 *) h_R2_Q3D->Project3D("yz");     add(h_R2_Q3D_yz, false,true,true,true);
    }

    if (reportDebug()) cout << "ParticlePairCombinedHistos::calculate() completed" << endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // calculate combinations of histograms with given coefficients
  //////////////////////////////////////////////////////////////////////////////////////////////////////
//  void ParticlePairCombinedHistos::calculateRR_Q3D(ClusterHistos * ch)
//  {
//    if (reportDebug()) "ParticlePairCombinedHistos::calculate() started" << endl;
//    AnalysisConfiguration & ac = *analysisConfiguration;
//
//    if (ac.fillQ3D)
//    {
//      h_RR2_Q3D->Divide(h_R2_Q3D, ch->h_cluster_Q3D,1.0,1.0);
//
//      h_RR2_Q3D_xy = (TH2 *) h_RR2_Q3D->Project3D("xy");     add(h_RR2_Q3D_xy, false,true,true,true);
//      h_RR2_Q3D_xz = (TH2 *) h_RR2_Q3D->Project3D("xz");     add(h_RR2_Q3D_xz, false,true,true,true);
//      h_RR2_Q3D_yz = (TH2 *) h_RR2_Q3D->Project3D("yz");     add(h_RR2_Q3D_yz, false,true,true,true);
//    }
//
//    if (reportDebug()) "ParticlePairCombinedHistos::calculate() completed" << endl;
//  }
