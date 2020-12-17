//
//  CollisionGeometryHistograms.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
// ===========================================================
#include "CollisionGeometryHistograms.hpp"
ClassImp(CollisionGeometryHistograms);

CollisionGeometryHistograms::CollisionGeometryHistograms(const TString & collectionName,
                                                         CollisionGeometryConfiguration * _configuration,
                                                         LogLevel  debugLevel)
:
Histograms(collectionName,nullptr,100,debugLevel),
configuration(_configuration),
nEventNcollGE1(0),
nEventNcollGE0(0),
h_crossSection(0),
h_b(0),
h_nPart(0),
h_nBinary(0),
h_nPartVsB_Prof(0),
h_nBinaryVsB_Prof(0),
h_nPartVsB(0),
h_nBinaryVsB(0),
h_xyNNIntVsB(0),
h_varXVsB_Prof(0),
h_varYVsB_Prof(0),
h_covXYVsB_Prof(0),
h_epsilonXVsB_Prof(0),
h_epsilonYVsB_Prof(0),
h_epsilonXYVsB_Prof(0),
h_psi2VsB_Prof(0),
h_psi2VsB(0)
{
  // no ops
}



// for now use the same boundaries for eta and y histogram
void CollisionGeometryHistograms::createHistograms()
{
  CollisionGeometryConfiguration & ac = *configuration;
  TString bn = ac.histoBaseName; bn += "_";

  //  1=events, 2=nBin>=1, 3=xSect
  h_crossSection       = createHistogram(bn+TString("Xsect"), 10, 0.0, 10.0, "label",  "Varia", 0, 1);
  double bMax = ac.maxB;
  h_crossSection->SetBinContent( 4, 3.1415927*bMax*bMax);
  h_crossSection->SetBinError  ( 4, 0.0);
  h_b                  = createHistogram(bn+TString("b"),                     ac.nBins_b,      ac.min_b,       ac.max_b,       "b (fm)",        "Counts", 0, 1);
  h_nPart              = createHistogram(bn+TString("nPart"),                 ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nBinary            = createHistogram(bn+TString("nBinary"),               ac.nBins_nBinary,ac.min_nBinary, ac.max_nBinary, "N_{Bin}",  "Counts", 0, 1);
  h_nPartVsB_Prof      = createProfile(bn+TString("nPartVsB_Prof"),           ac.nBins_b,     ac.min_b,  ac.max_b,  "b (fm)",  "<N_{Part}>",1);
  h_nBinaryVsB_Prof    = createProfile(bn+TString("nBinaryVsB_Prof"),         ac.nBins_b,     ac.min_b,  ac.max_b,  "b (fm)",  "<N_{Bin}>", 1);
  h_nPartSqVsB_Prof    = createProfile(bn+TString("nPartSqVsB_Prof"),         ac.nBins_b,     ac.min_b,  ac.max_b,  "b (fm)",  "<N_{Part}^{2}",1);
  h_nBinarySqVsB_Prof  = createProfile(bn+TString("nBinarySqVsB_Prof"),       ac.nBins_b,     ac.min_b,  ac.max_b,  "b (fm)",  "N_{Bin}^{2}", 1);
  h_nBinaryVsNPart_Prof   = createProfile(bn+TString("nBinaryVsNPart_Prof"),  ac.nBins_nPart, ac.min_nPart,   ac.max_nPart,  "N_{Part}", "<N_{Bin}>", 1);
  h_nBinarySqVsNPart_Prof = createProfile(bn+TString("nBinarySqVsNPart_Prof"),ac.nBins_nPart, ac.min_nPart,   ac.max_nPart,"N_{Part}", "<N_{Bin}^{2}>", 1);
  h_nPartVsB          = createHistogram(bn+TString("nPartVsB"),               ac.nBins_b,     ac.min_b,  ac.max_b,  ac.nBins_nPart,   ac.min_nPart,   ac.max_nPart,   "b (fm)", "N_{Part}", "Counts",0,1);
  h_nBinaryVsB        = createHistogram(bn+TString("nBinaryVsB"),             ac.nBins_b,     ac.min_b,  ac.max_b,  ac.nBins_nBinary, ac.min_nBinary, ac.max_nBinary, "b (fm)", "N_{Bin}",  "Counts",0,1);
  h_nPartRmsVsB       = createHistogram(bn+TString("nPartRmsVsB"),            ac.nBins_b,     ac.min_b,  ac.max_b,  "b (fm)", "RMS(N_{Part})",0,1);
  h_nBinaryRmsVsB     = createHistogram(bn+TString("nBinaryRmsVsB"),          ac.nBins_b,     ac.min_b,  ac.max_b,  "b (fm)", "RMS(N_{Bin})",0,1);
  h_nBinaryRmsVsNPart = createHistogram(bn+TString("nBinaryRmsVsNPart"),      ac.nBins_nPart, ac.min_nPart,   ac.max_nPart, "N_{Part}", "RMS(N_{Bin})",0,1);

  h_nPartOmegaVsB        = createHistogram(bn+TString("nPartOmegaVsB"),       ac.nBins_b,     ac.min_b,     ac.max_b,     "b (fm)",     "#omega(N_{part})");
  h_nBinaryOmegaVsB      = createHistogram(bn+TString("nBinaryOmegaVsB"),     ac.nBins_b,     ac.min_b,     ac.max_b,     "b (fm)",     "#omega(N_{Bin})");
  h_nBinaryOmegaVsNPart  = createHistogram(bn+TString("nBinaryOmegaVsNPart"), ac.nBins_nPart, ac.min_nPart, ac.max_nPart, "n_{Part}",   "#omega(N_{Bin})");
  h_nPartR2VsB           = createHistogram(bn+TString("nPartR2VsB"),          ac.nBins_b,     ac.min_b,     ac.max_b,     "b (fm)",     "R_{2}(N_{part})");
  h_nBinaryR2VsB         = createHistogram(bn+TString("nBinaryR2VsB"),        ac.nBins_b,     ac.min_b,     ac.max_b,     "b (fm)",     "R_{2}(N_{bin})");
  h_nBinaryR2VsNPart     = createHistogram(bn+TString("nBinaryR2VsNPart"),    ac.nBins_nPart, ac.min_nPart, ac.max_nPart, "n_{Part}",   "R_{2}(N_{bin})");

 
  // assume a maximum b of 17 fm for simplicity and slice it 10
  // equal fractions as bdb...

  bValues10[0] = 0.0;
  bValues10[1] = 3.39927;
  bValues10[2] = 4.79929;
  bValues10[3] = 6.79951;
  bValues10[4] = 8.39898;
  bValues10[5] = 9.79733;
  bValues10[6] = 10.9965;
  bValues10[7] = 11.9972;
  bValues10[8] = 12.9963;
  bValues10[9] = 13.7987;
  bValues10[10] = 14.7986;
  bValues10[11] = 20.0;

  int nSlices10 = 11;
  int nSlices100 = 100;
  //double xSlice = 16.5*16.5/10.0;  761.632
  double xSlice = 761.632/3.1415927; // PbPb total xsection in fm^2 divided Pi

  for (int i=0;i<=nSlices100;i++)
  {
  bValues100[i] = sqrt( xSlice* double(i)/100.0 );
  }

  h_nPartVsXsect_Prof     = createProfile(bn+TString("nPartVsXsect_Prof"),     nSlices10, bValues10, "b (fm)",  "<N_{Part}>",1);
  h_nPartSqVsXsect_Prof   = createProfile(bn+TString("nPartSqVsXsect_Prof"),   nSlices10, bValues10, "b (fm)",  "<N_{Part}^2>",1);
  h_nBinaryVsXsect_Prof   = createProfile(bn+TString("nBinaryVsXsect_Prof"),   nSlices10, bValues10, "b (fm)",  "<N_{Bin}>",1);
  h_nBinarySqVsXsect_Prof = createProfile(bn+TString("nBinarySqVsXsect_Prof"), nSlices10, bValues10, "b (fm)",  "<N_{Bin}^2>",1);
  h_nPartRmsVsXsect       = createHistogram(bn+TString("nPartRmsVsXsect"),     nSlices10, bValues10, "b (fm)",  "RMS(N_{Part})",    0,1);
  h_nBinaryRmsVsXsect     = createHistogram(bn+TString("nBinaryRmsVsXsect"),   nSlices10, bValues10, "b (fm)",  "RMS(N_{Bin})",     0,1);
  h_nPartOmegaVsXsect     = createHistogram(bn+TString("nPartOmegaVsXsect"),   nSlices10, bValues10, "b (fm)",  "#omega(N_{part})", 0,1);
  h_nBinaryOmegaVsXsect   = createHistogram(bn+TString("nBinaryOmegaVsXsect"), nSlices10, bValues10, "b (fm)",  "#omega(N_{Bin})",  0,1);
  h_nPartR2VsXsect        = createHistogram(bn+TString("nPartR2VsXsect"),      nSlices10, bValues10, "b (fm)",  "R_{2}(N_{part})",  0,1);
  h_nBinaryR2VsXsect      = createHistogram(bn+TString("nBinaryR2VsXsect"),    nSlices10, bValues10, "b (fm)",  "R_{2}(N_{bin})",   0,1);

  h_nPartVsXsect_Prof100     = createProfile(bn+TString("nPartVsXsect_Prof100"),     nSlices100, bValues100, "b (fm)",  "<N_{Part}>",1);
  h_nPartSqVsXsect_Prof100   = createProfile(bn+TString("nPartSqVsXsect_Prof100"),   nSlices100, bValues100, "b (fm)",  "<N_{Part}^2>",1);
  h_nBinaryVsXsect_Prof100   = createProfile(bn+TString("nBinaryVsXsect_Prof100"),   nSlices100, bValues100, "b (fm)",  "<N_{Bin}>",1);
  h_nBinarySqVsXsect_Prof100 = createProfile(bn+TString("nBinarySqVsXsect_Prof100"), nSlices100, bValues100, "b (fm)",  "<N_{Bin}^2>",1);
  h_nPartRmsVsXsect100       = createHistogram(bn+TString("nPartRmsVsXsect100"),     nSlices100, bValues100, "b (fm)",  "RMS(N_{Part})",    0,1);
  h_nBinaryRmsVsXsect100     = createHistogram(bn+TString("nBinaryRmsVsXsect100"),   nSlices100, bValues100, "b (fm)",  "RMS(N_{Bin})",     0,1);
  h_nPartOmegaVsXsect100     = createHistogram(bn+TString("nPartOmegaVsXsect100"),   nSlices100, bValues100, "b (fm)",  "#omega(N_{part})", 0,1);
  h_nBinaryOmegaVsXsect100   = createHistogram(bn+TString("nBinaryOmegaVsXsect100"), nSlices100, bValues100, "b (fm)",  "#omega(N_{Bin})",  0,1);
  h_nPartR2VsXsect100        = createHistogram(bn+TString("nPartR2VsXsect100"),      nSlices100, bValues100, "b (fm)",  "R_{2}(N_{part})",  0,1);
  h_nBinaryR2VsXsect100      = createHistogram(bn+TString("nBinaryR2VsXsect100"),    nSlices100, bValues100, "b (fm)",  "R_{2}(N_{bin})",   0,1);

  h_nPart_0_5       = createHistogram(bn+TString("nPart_0_5"),     ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_5_10      = createHistogram(bn+TString("nPart_5_10"),    ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_10_20     = createHistogram(bn+TString("nPart_10_20"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_20_30     = createHistogram(bn+TString("nPart_20_30"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_30_40     = createHistogram(bn+TString("nPart_30_40"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_40_50     = createHistogram(bn+TString("nPart_40_50"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_50_60     = createHistogram(bn+TString("nPart_50_60"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_60_70     = createHistogram(bn+TString("nPart_60_70"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_70_80     = createHistogram(bn+TString("nPart_70_80"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_80_90     = createHistogram(bn+TString("nPart_80_90"),   ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);
  h_nPart_90_100    = createHistogram(bn+TString("nPart_90_100"),  ac.nBins_nPart,  ac.min_nPart,   ac.max_nPart,   "N_{Part}", "Counts", 0, 1);

  h_nBinary_0_5       = createHistogram(bn+TString("nBinary_0_5"),     ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_5_10      = createHistogram(bn+TString("nBinary_5_10"),    ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_10_20     = createHistogram(bn+TString("nBinary_10_20"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_20_30     = createHistogram(bn+TString("nBinary_20_30"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_30_40     = createHistogram(bn+TString("nBinary_30_40"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_40_50     = createHistogram(bn+TString("nBinary_40_50"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_50_60     = createHistogram(bn+TString("nBinary_50_60"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_60_70     = createHistogram(bn+TString("nBinary_60_70"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_70_80     = createHistogram(bn+TString("nBinary_70_80"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_80_90     = createHistogram(bn+TString("nBinary_80_90"),   ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);
  h_nBinary_90_100    = createHistogram(bn+TString("nBinary_90_100"),  ac.nBins_nBinary,  ac.min_nBinary,   ac.max_nBinary,   "N_{Bin}", "Counts", 0, 1);

// =============
  // plots to make
  // to compare to Loizides PHYSICAL REVIEW C 94, 024914 (2016)
  // Calculated total cross sections for PbPb and pPb collisions as a function of σNN
  // Ncoll/(Npart/2) vs Npart
  // S(fm^2) vs Npart
  // eccentricity epsilon2 vs npart
  // triangularity epsilon3 vs npart
  

  h_xyNNIntVsB       = createHistogram(bn+TString("xyNNIntVsB"),     ac.nBins_b,  ac.min_b,  ac.max_b,  40, -10.0, 10.0, 40, -10.0, 10.0, "b", "x", "y", "Counts",0,1);

  h_varXVsB_Prof     = createProfile(bn+"varXVsB_Prof",     ac.nBins_b,  ac.min_b,  ac.max_b,  "b",  "Var[x]", 1);
  h_varYVsB_Prof     = createProfile(bn+"varYVsB_Prof",     ac.nBins_b,  ac.min_b,  ac.max_b,  "b",  "Var[y]", 1);
  h_covXYVsB_Prof    = createProfile(bn+"covXYVsB_Prof",    ac.nBins_b,  ac.min_b,  ac.max_b,  "b",  "Cov[x,y]", 1);
  h_epsilonXVsB_Prof = createProfile(bn+"epsilonXVsB_Prof", ac.nBins_b,  ac.min_b,  ac.max_b,  "b",  "#epsilon_{x}", 1);
  h_epsilonYVsB_Prof = createProfile(bn+"epsilonYVsB_Prof", ac.nBins_b,  ac.min_b,  ac.max_b,  "b",  "#epsilon_{y}", 1);
  h_epsilonXYVsB_Prof  = createProfile(bn+"epsilonXYVsB_Prof",  ac.nBins_b,  ac.min_b,  ac.max_b,  "b",  "|#epsilon|", 1);
  h_epsilonXYVsB       = createHistogram(bn+TString("epsilonXYVsB"), ac.nBins_b,  ac.min_b,  ac.max_b,  40, -1.0, 1.0,  40, -1.0, 1.0, "b", "#epsilon_{x}", "#epsilon_{y}", "Counts",0,1);;

  h_psi2VsB_Prof  = createProfile(bn+"psi2VsB_Prof",         ac.nBins_b,  ac.min_b,  ac.max_b,  "b",  "#psi_{2}", 1);
  h_psi2VsB       = createHistogram(bn+TString("psi2VsB"),   ac.nBins_b,  ac.min_b,  ac.max_b,  40, -TMath::Pi(), TMath::Pi(), "b", "#psi_2", "Counts",0,1);
}

//________________________________________________________________________
void CollisionGeometryHistograms::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportFatal()) cout << "-Fatal- Attempting to load CollisionGeometryHistograms from an invalid file pointer" << endl;
    return;
    }
  CollisionGeometryConfiguration & ac = *configuration;
  TString bn = ac.histoBaseName; bn += "_";
  h_crossSection      = loadH1(inputFile,      bn+TString("Xsect") );
  h_b                 = loadH1(inputFile,      bn+TString("b") );
  h_nPart             = loadH1(inputFile,      bn+TString("nPart") );
  h_nBinary           = loadH1(inputFile,      bn+TString("nBinary") );
  h_nPartVsB_Prof     = loadProfile(inputFile, bn+TString("nPartVsB_Prof"));
  h_nPartSqVsB_Prof   = loadProfile(inputFile, bn+TString("nPartSqVsB_Prof"));
  h_nBinaryVsB_Prof   = loadProfile(inputFile, bn+TString("nBinaryVsB_Prof"));
  h_nBinarySqVsB_Prof = loadProfile(inputFile, bn+TString("nBinarySqVsB_Prof"));
  h_nPartVsB          = loadH2(inputFile,      bn+TString("nPartVsB"));
  //h_nPartVsB          = loadH2(inputFile,      bn+TString("nPartVsB"));
  h_nBinaryVsB        = loadH2(inputFile,      bn+TString("nBinaryVsB"));

  h_nBinaryVsNPart_Prof   = loadProfile(inputFile, bn+TString("nBinaryVsNPart_Prof"));
  h_nBinarySqVsNPart_Prof = loadProfile(inputFile, bn+TString("nBinarySqVsNPart_Prof"));
  h_nPartRmsVsB           = loadH1(inputFile,      bn+TString("nPartRmsVsB"));
  h_nBinaryRmsVsB         = loadH1(inputFile,      bn+TString("nBinaryRmsVsB"));
  h_nBinaryRmsVsNPart     = loadH1(inputFile,      bn+TString("nBinaryRmsVsNPart"));

  h_nPartOmegaVsB         = loadH1(inputFile,      bn+TString("nPartOmegaVsB"));
  h_nBinaryOmegaVsB       = loadH1(inputFile,      bn+TString("nPartOmegaVsB"));
  h_nBinaryOmegaVsNPart   = loadH1(inputFile,      bn+TString("nBinaryOmegaVsNPart"));
  h_nPartR2VsB            = loadH1(inputFile,      bn+TString("nPartR2VsB"));
  h_nBinaryR2VsB          = loadH1(inputFile,      bn+TString("nBinaryR2VsB"));
  h_nBinaryR2VsNPart      = loadH1(inputFile,      bn+TString("nBinaryR2VsNPart"));


  h_nPartVsXsect_Prof     = loadProfile(inputFile, bn+TString("nPartVsXsect_Prof"));
  h_nPartSqVsXsect_Prof   = loadProfile(inputFile, bn+TString("nPartSqVsXsect_Prof"));
  h_nBinaryVsXsect_Prof   = loadProfile(inputFile, bn+TString("nBinaryVsXsect_Prof"));
  h_nBinarySqVsXsect_Prof = loadProfile(inputFile, bn+TString("nBinarySqVsXsect_Prof"));
  h_nPartRmsVsXsect       = loadH1(inputFile,      bn+TString("nPartRmsVsXsect"));
  h_nBinaryRmsVsXsect     = loadH1(inputFile,      bn+TString("nBinaryRmsVsXsect"));
  h_nPartOmegaVsXsect     = loadH1(inputFile,      bn+TString("nPartOmegaVsXsect"));
  h_nBinaryOmegaVsXsect   = loadH1(inputFile,      bn+TString("nBinaryOmegaVsXsect"));
  h_nPartR2VsXsect        = loadH1(inputFile,      bn+TString("nPartR2VsXsect"));
  h_nBinaryR2VsXsect      = loadH1(inputFile,      bn+TString("nBinaryR2VsXsect"));

  h_nPartVsXsect_Prof100     = loadProfile(inputFile, bn+TString("nPartVsXsect_Prof100"));
  h_nPartSqVsXsect_Prof100   = loadProfile(inputFile, bn+TString("nPartSqVsXsect_Prof100"));
  h_nBinaryVsXsect_Prof100   = loadProfile(inputFile, bn+TString("nBinaryVsXsect_Prof100"));
  h_nBinarySqVsXsect_Prof100 = loadProfile(inputFile, bn+TString("nBinarySqVsXsect_Prof100"));
  h_nPartRmsVsXsect100       = loadH1(inputFile,      bn+TString("nPartRmsVsXsect100"));
  h_nBinaryRmsVsXsect100     = loadH1(inputFile,      bn+TString("nBinaryRmsVsXsect100"));
  h_nPartOmegaVsXsect100     = loadH1(inputFile,      bn+TString("nPartOmegaVsXsect100"));
  h_nBinaryOmegaVsXsect100   = loadH1(inputFile,      bn+TString("nBinaryOmegaVsXsect100"));
  h_nPartR2VsXsect100        = loadH1(inputFile,      bn+TString("nPartR2VsXsect100"));
  h_nBinaryR2VsXsect100      = loadH1(inputFile,      bn+TString("nBinaryR2VsXsect100"));

  h_nPart_0_5       = loadH1(inputFile,      bn+TString("nPart_0_5"));
  h_nPart_5_10      = loadH1(inputFile,      bn+TString("nPart_5_10"));
  h_nPart_10_20     = loadH1(inputFile,      bn+TString("nPart_10_20"));
  h_nPart_20_30     = loadH1(inputFile,      bn+TString("nPart_20_30"));
  h_nPart_30_40     = loadH1(inputFile,      bn+TString("nPart_30_40"));
  h_nPart_40_50     = loadH1(inputFile,      bn+TString("nPart_40_50"));
  h_nPart_50_60     = loadH1(inputFile,      bn+TString("nPart_50_60"));
  h_nPart_60_70     = loadH1(inputFile,      bn+TString("nPart_60_70"));
  h_nPart_70_80     = loadH1(inputFile,      bn+TString("nPart_70_80"));
  h_nPart_80_90     = loadH1(inputFile,      bn+TString("nPart_80_90"));
  h_nPart_90_100    = loadH1(inputFile,      bn+TString("nPart_90_100"));

  h_nBinary_0_5       = loadH1(inputFile,      bn+TString("nBinary_0_5"));
  h_nBinary_5_10      = loadH1(inputFile,      bn+TString("nBinary_5_10"));
  h_nBinary_10_20     = loadH1(inputFile,      bn+TString("nBinary_10_20"));
  h_nBinary_20_30     = loadH1(inputFile,      bn+TString("nBinary_20_30"));
  h_nBinary_30_40     = loadH1(inputFile,      bn+TString("nBinary_30_40"));
  h_nBinary_40_50     = loadH1(inputFile,      bn+TString("nBinary_40_50"));
  h_nBinary_50_60     = loadH1(inputFile,      bn+TString("nBinary_50_60"));
  h_nBinary_60_70     = loadH1(inputFile,      bn+TString("nBinary_60_70"));
  h_nBinary_70_80     = loadH1(inputFile,      bn+TString("nBinary_70_80"));
  h_nBinary_80_90     = loadH1(inputFile,      bn+TString("nBinary_80_90"));
  h_nBinary_90_100    = loadH1(inputFile,      bn+TString("nBinary_90_100"));

  h_xyNNIntVsB         = loadH3(inputFile,     bn+TString("xyNNIntVsB"));

  h_varXVsB_Prof       = loadProfile(inputFile,bn+TString("varXVsB_Prof"));
  h_varYVsB_Prof       = loadProfile(inputFile,bn+TString("varYVsB_Prof"));
  h_covXYVsB_Prof      = loadProfile(inputFile,bn+TString("covXYVsB_Prof"));
  h_epsilonXVsB_Prof   = loadProfile(inputFile,bn+TString("epsilonXVsB_Prof"));
  h_epsilonYVsB_Prof   = loadProfile(inputFile,bn+TString("epsilonYVsB_Prof"));
  h_epsilonXYVsB_Prof  = loadProfile(inputFile,bn+TString("epsilonXYVsB_Prof"));
  h_epsilonXYVsB       = loadH3(inputFile,     bn+TString("epsilonXYVsB"));
  h_psi2VsB_Prof       = loadProfile(inputFile,bn+TString("psi2VsB_Prof"));
  h_psi2VsB            = loadH2(inputFile,     bn+TString("psi2VsB"));
}


void CollisionGeometryHistograms::fill(CollisionGeometry * collisionGeometry, double weight)
{
  double impactPar = collisionGeometry->b;
  double nPart     = collisionGeometry->nParticipant;
  double nBinary   = collisionGeometry->nBinary;

  nEventNcollGE0++; // events counted....
  h_crossSection->Fill(0.5);
  if (nBinary<1) return;
  h_crossSection->Fill(1.5);
  nEventNcollGE1++; // at least one binary collision

  h_b        ->Fill(impactPar,   weight);
  h_nPart    ->Fill(nPart,       weight);
  h_nBinary  ->Fill(nBinary,     weight);

  h_nPartVsB_Prof    ->Fill(impactPar, nPart,   weight);
  h_nPartSqVsB_Prof  ->Fill(impactPar, nPart*nPart,   weight);
  h_nBinaryVsB_Prof  ->Fill(impactPar, nBinary, weight);
  h_nBinarySqVsB_Prof->Fill(impactPar, nBinary*nBinary, weight);
  h_nPartVsB        ->Fill(impactPar, nPart,   weight);
  h_nBinaryVsB      ->Fill(impactPar, nBinary, weight);
  h_nBinaryVsNPart_Prof   ->Fill(nPart, nBinary,   weight);
  h_nBinarySqVsNPart_Prof ->Fill(nPart, nBinary*nBinary,   weight);

  h_nPartVsXsect_Prof     ->Fill(impactPar, nPart,   weight);
  h_nPartSqVsXsect_Prof   ->Fill(impactPar, nPart*nPart,   weight);
  h_nBinaryVsXsect_Prof   ->Fill(impactPar, nBinary, weight);
  h_nBinarySqVsXsect_Prof ->Fill(impactPar, nBinary*nBinary, weight);

  h_nPartVsXsect_Prof100     ->Fill(impactPar, nPart,   weight);
  h_nPartSqVsXsect_Prof100   ->Fill(impactPar, nPart*nPart,   weight);
  h_nBinaryVsXsect_Prof100   ->Fill(impactPar, nBinary, weight);
  h_nBinarySqVsXsect_Prof100 ->Fill(impactPar, nBinary*nBinary, weight);

  if (impactPar<bValues10[1] && impactPar>=bValues10[0])  h_nPart_0_5->Fill(nPart,       weight);
  if (impactPar<bValues10[2] && impactPar>=bValues10[1])  h_nPart_5_10->Fill(nPart,       weight);
  if (impactPar<bValues10[3] && impactPar>=bValues10[2])  h_nPart_10_20->Fill(nPart,       weight);
  if (impactPar<bValues10[4] && impactPar>=bValues10[3])  h_nPart_20_30->Fill(nPart,       weight);
  if (impactPar<bValues10[5] && impactPar>=bValues10[4])  h_nPart_30_40->Fill(nPart,       weight);
  if (impactPar<bValues10[6] && impactPar>=bValues10[5])  h_nPart_40_50->Fill(nPart,       weight);
  if (impactPar<bValues10[7] && impactPar>=bValues10[6])  h_nPart_50_60->Fill(nPart,       weight);
  if (impactPar<bValues10[8] && impactPar>=bValues10[7])  h_nPart_60_70->Fill(nPart,       weight);
  if (impactPar<bValues10[9] && impactPar>=bValues10[8])  h_nPart_70_80->Fill(nPart,       weight);
  if (impactPar<bValues10[10] && impactPar>=bValues10[9])  h_nPart_80_90->Fill(nPart,       weight);
  if (impactPar<bValues10[11] && impactPar>=bValues10[10])  h_nPart_90_100->Fill(nPart,       weight);


  if (impactPar<bValues10[1] && impactPar>=bValues10[0])  h_nBinary_0_5->Fill(nBinary,       weight);
  if (impactPar<bValues10[2] && impactPar>=bValues10[1])  h_nBinary_5_10->Fill(nBinary,       weight);
  if (impactPar<bValues10[3] && impactPar>=bValues10[2])  h_nBinary_10_20->Fill(nBinary,       weight);
  if (impactPar<bValues10[4] && impactPar>=bValues10[3])  h_nBinary_20_30->Fill(nBinary,       weight);
  if (impactPar<bValues10[5] && impactPar>=bValues10[4])  h_nBinary_30_40->Fill(nBinary,       weight);
  if (impactPar<bValues10[6] && impactPar>=bValues10[5])  h_nBinary_40_50->Fill(nBinary,       weight);
  if (impactPar<bValues10[7] && impactPar>=bValues10[6])  h_nBinary_50_60->Fill(nBinary,       weight);
  if (impactPar<bValues10[8] && impactPar>=bValues10[7])  h_nBinary_60_70->Fill(nBinary,       weight);
  if (impactPar<bValues10[9] && impactPar>=bValues10[8])  h_nBinary_70_80->Fill(nBinary,       weight);
  if (impactPar<bValues10[10] && impactPar>=bValues10[9])  h_nBinary_80_90->Fill(nBinary,       weight);
  if (impactPar<bValues10[11] && impactPar>=bValues10[10])  h_nBinary_90_100->Fill(nBinary,       weight);


  calculateMomentsFromBinaryColl(collisionGeometry);
  calculateMomentsFromParticipants(collisionGeometry);




//  h_varXVsB_Prof       ->Fill(impactPar, varX,   weight);
//  h_varYVsB_Prof       ->Fill(impactPar, varY,   weight);
//  h_covXYVsB_Prof      ->Fill(impactPar, varXY,  weight);
//  h_epsilonXVsB_Prof   ->Fill(impactPar, epsX,   weight);
//  h_epsilonYVsB_Prof   ->Fill(impactPar, epsY,   weight);
//  h_epsilonXYVsB_Prof  ->Fill(impactPar, epsMod, weight);
//  h_epsilonXYVsB       ->Fill(impactPar, epsX, epsY, weight);
//  h_psi2VsB_Prof       ->Fill(impactPar, psi2,   weight);
//  h_psi2VsB            ->Fill(impactPar, psi2,   weight);

}

void CollisionGeometryHistograms::calculateDerivedHistograms()
{
  double nEvents      = h_crossSection->GetBinContent(1);
  double nEventsWColl = h_crossSection->GetBinContent(2); double eNEventWColl = h_crossSection->GetBinError(2);
  double area = h_crossSection->GetBinContent(4);
  double ratio1  = nEvents>0      ? nEventsWColl/nEvents : 1.0e-30;
  double eRatio1 = nEvents>0      ? eNEventWColl/nEvents : 1.0e-30;
  double ratio2  = nEventNcollGE0 ? double(nEventNcollGE1)/double(nEventNcollGE0) : 1.0e-30;
  double eRatio2 = nEventNcollGE0 ? sqrt(double(nEventNcollGE1))/double(nEventNcollGE0) : 1.0e-30;
  double xSect1  = ratio1*area;
  double eXSect1 = eRatio1*area;
  double xSect2  = ratio2*area;
  double eXSect2 = eRatio2*area;

  h_crossSection->SetBinContent(5,ratio1); h_crossSection->SetBinError(5,eRatio1);
  h_crossSection->SetBinContent(6,ratio2); h_crossSection->SetBinError(6,eRatio2);
  h_crossSection->SetBinContent(7,xSect1); h_crossSection->SetBinError(7,eXSect1);
  h_crossSection->SetBinContent(8,xSect2); h_crossSection->SetBinError(8,eXSect2);

  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()                         nEvents: " <<  nEvents << endl;
  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()                    nEventsWColl: " <<  nEventsWColl << endl;
  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()            nEventsWColl/nEvents: " <<  ratio1 << " +-" << eRatio1 << endl;
  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()                  nEventNcollGE0: " <<  nEventNcollGE0 << endl;
  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()                  nEventNcollGE1: " <<  nEventNcollGE1 << endl;
  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()   nEventNcollGE1/nEventNcollGE0: " <<  ratio2 << " +-" << eRatio2 << endl;
  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()            xSect1 = ratio1*area: " <<  xSect1 << " +-" << eXSect1 << endl;
  cout << "  CollisionGeometryHistograms::calculateDerivedHistograms()            xSect2 = ratio2*area: " <<  xSect2 << " +-" << eXSect2 << endl;

  calculateRms(h_nPartVsB_Prof,       h_nPartSqVsB_Prof,       h_nPartRmsVsB,       h_nPartOmegaVsB,       h_nPartR2VsB);
  calculateRms(h_nBinaryVsB_Prof,     h_nBinarySqVsB_Prof,     h_nBinaryRmsVsB,     h_nBinaryOmegaVsB,     h_nBinaryR2VsB);
  calculateRms(h_nBinaryVsNPart_Prof, h_nBinarySqVsNPart_Prof, h_nBinaryRmsVsNPart, h_nBinaryOmegaVsNPart, h_nBinaryR2VsNPart);
  calculateRms(h_nPartVsXsect_Prof,   h_nPartSqVsXsect_Prof,   h_nPartRmsVsXsect,   h_nPartOmegaVsXsect,   h_nPartR2VsXsect);
  calculateRms(h_nBinaryVsXsect_Prof, h_nBinarySqVsXsect_Prof, h_nBinaryRmsVsXsect, h_nBinaryOmegaVsXsect, h_nBinaryR2VsXsect);
}

void CollisionGeometryHistograms::calculateRms(TProfile * h1, TProfile * h1Sq,   TH1* h1Rms, TH1* h1Omega, TH1* h1R2)
{
  double v, ev, v2, ev2, rms, erms, omega, eomega, R2, eR2;
  int nBins = h1->GetNbinsX();
  for (int iBin=1; iBin<nBins; iBin++)
  {
  v = h1->GetBinContent(iBin);
  ev = 0.0;
  v2 = h1Sq->GetBinContent(iBin);
  ev2 = 0.0;
  rms = v2 - v*v;
  if (rms>0) rms = sqrt(rms);
  erms = 0.0;


  if (v>0)
    {
    omega  = rms*rms/v;
    eomega = 0.0;
    R2 = (v2-v)/v/v - 1.0;
    eR2 = 0;
    }
  else
    {
    omega  = 0.0;
    eomega = 0.0;
    R2     = 0.0;
    eR2    = 0;
    }

  //cout << " iBin: " << iBin << " v:" << v << " v2:" << v2 << " rms:" << rms << " omega:" << omega << " R2:" << R2 << endl;

  h1Rms->SetBinContent(iBin, rms);
  h1Rms->SetBinError  (iBin, erms);

  h1Omega->SetBinContent(iBin, omega);
  h1Omega->SetBinError(iBin, eomega);

  h1R2->SetBinContent(iBin, R2);
  h1R2->SetBinError(iBin, eR2);
  }
}


void CollisionGeometryHistograms::calculateInelCrossSection()
{
  
}


// calculate the moments from the binary colls and save them
// in the geometry event...
// calculate eccentricities from n=1 to n=10;
void CollisionGeometryHistograms::calculateMomentsFromBinaryColl(CollisionGeometry * collisionGeometry)
{
  double sum;
  double meanX, meanY;
  double varX, varY, varXY;
  double w, n;
  double r, rw, phi;
  double cphiN[10], sphiN[10], rN[10], psiN[10], eccN[10];
  double area;
  double x, x2, y, y2, xy;
  double epsX;
  double epsY;
  double epsDenom, epsMod;
  double psi2;

  // mean and variances
  int nBinary   = collisionGeometry->nBinary;
  sum    = 0.0;
  meanX  = 0.0;
  meanY  = 0.0;
  varX   = 0.0;
  varY   = 0.0;
  varXY  = 0.0;
  epsDenom  = -999.0;
  epsX      = -999.0;
  epsY      = -999.0;
  epsMod    = -999.0;
  psi2      = -999.0;
  area      = -999.0;
  // compute 1st and 2nd moments
  for (int iBinary=0; iBinary<nBinary; iBinary++)
    {
    w  = 1.0;
    x  = collisionGeometry->x[iBinary];  x2 = x*x;
    y  = collisionGeometry->y[iBinary];  y2 = y*y;
    xy = x*y;
    sum += w;
    meanX += w*x;
    meanY += w*y;
    varX  += w*x2;
    varY  += w*y2;
    varXY += w*xy;
    }
  if (sum>0)
    {
    meanX /= sum;
    meanY /= sum;
    varX  /= sum; varX  -= meanX*meanX;
    varY  /= sum; varY  -= meanY*meanY;
    varXY /= sum; varXY -= meanX*meanY;
    epsDenom = varX + varY;
    if (epsDenom>0)
      {
      epsX     = (varY - varX)/epsDenom;
      epsY     = 2*varXY/epsDenom;
      epsMod   = sqrt(epsX*epsX + epsY*epsY);
      psi2     = atan2(epsY,epsX);
      }
    }

  // compute n>=0 moments with weights...
  for (int iN=1; iN<10; iN++)
  {
  cphiN[iN] = 0.0; sphiN[iN] = 0.0; rN[iN] = 0.0; psiN[iN] = 0.0; eccN[iN] = 0.0;
  }
  sum = 0.0;
  for (int iBinary=0; iBinary<nBinary; iBinary++)
  {
  sum  += 1.0;
  x  = collisionGeometry->x[iBinary] - meanX;
  y  = collisionGeometry->y[iBinary] - meanY;
  r  = TMath::Sqrt(x*x+y*y);
  phi = TMath::ATan2(y,x);
  for (int iN=1; iN<10; iN++)
    {
    n = iN+1;
    w = (iN==0) ? 3 : n;
    rw = TMath::Power(r,w);
    cphiN[iN] += rw*TMath::Cos(n*phi);
    sphiN[iN] += rw*TMath::Sin(n*phi);
    rN[iN]    += rw;
    }
  }
  for (int iN=1; iN<10; iN++)
  {
  psiN[iN] = TMath::ATan2(sphiN[iN],cphiN[iN]) + TMath::Pi()/double(iN);
  eccN[iN] = TMath::Sqrt(sphiN[iN]*sphiN[iN]  + cphiN[iN]*cphiN[iN]) / rN[iN];
  }

  // copy values to the geometry event
  collisionGeometry->binaryMoments. meanX  =  meanX;
  collisionGeometry->binaryMoments. meanY  =  meanY;
  collisionGeometry->binaryMoments. varX   =  varX;
  collisionGeometry->binaryMoments. varY   =  varY;
  collisionGeometry->binaryMoments. varXY  =  varXY;
  collisionGeometry->binaryMoments. epsX   =  epsX;
  collisionGeometry->binaryMoments. epsY   =  epsY;
  collisionGeometry->binaryMoments. epsMod =  epsMod;
  collisionGeometry->binaryMoments. psi2   =  psi2;
  collisionGeometry->binaryMoments. area   =  area;
  for (int iN=0;iN<10;iN++)
  {
  collisionGeometry->binaryMoments. cphiN[iN]  =  cphiN[iN];
  collisionGeometry->binaryMoments. sphiN[iN]  =  sphiN[iN];
  collisionGeometry->binaryMoments. sphiN[iN]  =  sphiN[iN];
  collisionGeometry->binaryMoments. psiN[iN]   =  psiN[iN];
  collisionGeometry->binaryMoments. eccN[iN]   =  eccN[iN];
  }
}

// calculate the moments from the participants and save them
// in the geometry event...
// calculate eccentricities from n=1 to n=10;
void CollisionGeometryHistograms::calculateMomentsFromParticipants(CollisionGeometry * collisionGeometry)
{
  double sum;
  double meanX, meanY;
  double varX, varY, varXY;
  double w, n;
  double r, rw, phi;
  double cphiN[10], sphiN[10], rN[10], psiN[10], eccN[10];
  double area;
  double x, x2, y, y2, xy;
  double epsX;
  double epsY;
  double epsDenom, epsMod;
  double psi2;

  // mean and variances
  int nA = collisionGeometry->nucleusA->nNucleons;
  int nB = collisionGeometry->nucleusA->nNucleons;
  sum    = 0.0;
  meanX  = 0.0;
  meanY  = 0.0;
  varX   = 0.0;
  varY   = 0.0;
  varXY  = 0.0;
  epsDenom  = -999.0;
  epsX      = -999.0;
  epsY      = -999.0;
  epsMod    = -999.0;
  psi2      = -999.0;
  area      = -999.0;
  // compute 1st and 2nd moments
  // based on the position of the nucleons at impact
  for (int iA=0; iA<nA; iA++)
    {
    w  = 1.0;
    x  = collisionGeometry->nucleusA->getNucleon(iA)->x;  x2 = x*x;
    y  = collisionGeometry->nucleusA->getNucleon(iA)->y;  y2 = y*y;
    xy = x*y;
    sum += w;
    meanX += w*x;
    meanY += w*y;
    varX  += w*x2;
    varY  += w*y2;
    varXY += w*xy;
    }
  for (int iB=0; iB<nB; iB++)
    {
    w  = 1.0;
    x  = collisionGeometry->nucleusB->getNucleon(iB)->x;  x2 = x*x;
    y  = collisionGeometry->nucleusB->getNucleon(iB)->y;  y2 = y*y;
    xy = x*y;
    sum += w;
    meanX += w*x;
    meanY += w*y;
    varX  += w*x2;
    varY  += w*y2;
    varXY += w*xy;
    }


  if (sum>0)
    {
    meanX /= sum;
    meanY /= sum;
    varX  /= sum; varX  -= meanX*meanX;
    varY  /= sum; varY  -= meanY*meanY;
    varXY /= sum; varXY -= meanX*meanY;
    epsDenom = varX + varY;
    if (epsDenom>0)
      {
      epsX     = (varY - varX)/epsDenom;
      epsY     = 2*varXY/epsDenom;
      epsMod   = sqrt(epsX*epsX + epsY*epsY);
      psi2     = atan2(epsY,epsX);
      }
    }

  // compute n>=0 moments with weights...
  for (int iN=1; iN<10; iN++)
  {
  cphiN[iN] = 0.0; sphiN[iN] = 0.0; rN[iN] = 0.0; psiN[iN] = 0.0; eccN[iN] = 0.0;
  }
  sum = 0.0;
  for (int iA=0; iA<nA; iA++)
  {
  sum  += 1.0;
  x  = collisionGeometry->nucleusA->getNucleon(iA)->x - meanX;
  y  = collisionGeometry->nucleusA->getNucleon(iA)->y - meanY;
  r  = TMath::Sqrt(x*x+y*y);
  phi = TMath::ATan2(y,x);
  for (int iN=1; iN<10; iN++)
    {
    n = iN+1;
    w = (iN==0) ? 3 : n;
    rw = TMath::Power(r,w);
    cphiN[iN] += rw*TMath::Cos(n*phi);
    sphiN[iN] += rw*TMath::Sin(n*phi);
    rN[iN]    += rw;
    }
  }
  for (int iB=0; iB<nA; iB++)
  {
  sum  += 1.0;
  x  = collisionGeometry->nucleusB->getNucleon(iB)->x - meanX;
  y  = collisionGeometry->nucleusB->getNucleon(iB)->y - meanY;
  r  = TMath::Sqrt(x*x+y*y);
  phi = TMath::ATan2(y,x);
  for (int iN=1; iN<10; iN++)
    {
    n = iN+1;
    w = (iN==0) ? 3 : n;
    rw = TMath::Power(r,w);
    cphiN[iN] += rw*TMath::Cos(n*phi);
    sphiN[iN] += rw*TMath::Sin(n*phi);
    rN[iN]    += rw;
    }
  }
  for (int iN=1; iN<10; iN++)
  {
  psiN[iN] = TMath::ATan2(sphiN[iN],cphiN[iN]) + TMath::Pi()/double(iN);
  eccN[iN] = TMath::Sqrt(sphiN[iN]*sphiN[iN]  + cphiN[iN]*cphiN[iN]) / rN[iN];
  }

  // copy values to the geometry event
  collisionGeometry->participantMoments. meanX  =  meanX;
  collisionGeometry->participantMoments. meanY  =  meanY;
  collisionGeometry->participantMoments. varX   =  varX;
  collisionGeometry->participantMoments. varY   =  varY;
  collisionGeometry->participantMoments. varXY  =  varXY;
  collisionGeometry->participantMoments. epsX   =  epsX;
  collisionGeometry->participantMoments. epsY   =  epsY;
  collisionGeometry->participantMoments. epsMod =  epsMod;
  collisionGeometry->participantMoments. psi2   =  psi2;
  collisionGeometry->participantMoments. area   =  area;
  for (int iN=0;iN<10;iN++)
  {
  collisionGeometry->participantMoments. cphiN[iN]  =  cphiN[iN];
  collisionGeometry->participantMoments. sphiN[iN]  =  sphiN[iN];
  collisionGeometry->participantMoments. sphiN[iN]  =  sphiN[iN];
  collisionGeometry->participantMoments. psiN[iN]   =  psiN[iN];
  collisionGeometry->participantMoments. eccN[iN]   =  eccN[iN];
  }
}
