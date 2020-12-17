// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_CollisionGeometryHistograms
#define WAC_CollisionGeometryHistograms
#include "Histograms.hpp"
#include "CollisionGeometry.hpp"
#include "CollisionGeometryConfiguration.hpp"

class CollisionGeometryHistograms : public Histograms
{
public:

  CollisionGeometryHistograms(const TString & collectionName,
                              CollisionGeometryConfiguration * configuration,
                              LogLevel  debugLevel);
  virtual ~CollisionGeometryHistograms() { }
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void fill(CollisionGeometry * collisionGeometry, double weight);
  virtual void calculateDerivedHistograms();
  virtual void calculateRms(TProfile * h1, TProfile * h1Sq, TH1* h1Rms, TH1* h1Omega, TH1* h1R2);
  virtual void calculateInelCrossSection();
  virtual void calculateMomentsFromBinaryColl(CollisionGeometry * collisionGeometry);
  virtual void calculateMomentsFromParticipants(CollisionGeometry * collisionGeometry);
  
  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////

  CollisionGeometryConfiguration * configuration;

  long nEventNcollGE1;
  long nEventNcollGE0;

  TH1      * h_crossSection;

  TH1      * h_b;
  TH1      * h_nPart;
  TH1      * h_nBinary;
  TProfile * h_nPartVsB_Prof;
  TProfile * h_nPartSqVsB_Prof;
  TProfile * h_nBinaryVsB_Prof;
  TProfile * h_nBinarySqVsB_Prof;
  TProfile * h_nBinaryVsNPart_Prof;
  TProfile * h_nBinarySqVsNPart_Prof;
  TH2      * h_nPartVsB;
  TH2      * h_nBinaryVsB;

  TH1      * h_nPartRmsVsB;
  TH1      * h_nBinaryRmsVsB;
  TH1      * h_nBinaryRmsVsNPart;
  TH1      * h_nPartOmegaVsB;
  TH1      * h_nBinaryOmegaVsB;
  TH1      * h_nBinaryOmegaVsNPart;
  TH1      * h_nPartR2VsB;
  TH1      * h_nBinaryR2VsB;
  TH1      * h_nBinaryR2VsNPart;

  TProfile * h_nPartVsXsect_Prof;
  TProfile * h_nPartSqVsXsect_Prof;
  TProfile * h_nBinaryVsXsect_Prof;
  TProfile * h_nBinarySqVsXsect_Prof;
  TH1      * h_nPartRmsVsXsect;
  TH1      * h_nBinaryRmsVsXsect;
  TH1      * h_nPartOmegaVsXsect;
  TH1      * h_nBinaryOmegaVsXsect;
  TH1      * h_nPartR2VsXsect;
  TH1      * h_nBinaryR2VsXsect;

  TProfile * h_nPartVsXsect_Prof100;
  TProfile * h_nPartSqVsXsect_Prof100;
  TProfile * h_nBinaryVsXsect_Prof100;
  TProfile * h_nBinarySqVsXsect_Prof100;
  TH1      * h_nPartRmsVsXsect100;
  TH1      * h_nBinaryRmsVsXsect100;
  TH1      * h_nPartOmegaVsXsect100;
  TH1      * h_nBinaryOmegaVsXsect100;
  TH1      * h_nPartR2VsXsect100;
  TH1      * h_nBinaryR2VsXsect100;

  double bValues10[12];
  double bValues100[101];

  TH1      * h_nPart_0_5;
  TH1      * h_nPart_5_10;
  TH1      * h_nPart_10_20;
  TH1      * h_nPart_20_30;
  TH1      * h_nPart_30_40;
  TH1      * h_nPart_40_50;
  TH1      * h_nPart_50_60;
  TH1      * h_nPart_60_70;
  TH1      * h_nPart_70_80;
  TH1      * h_nPart_80_90;
  TH1      * h_nPart_90_100;

  TH1      * h_nBinary_0_5;
  TH1      * h_nBinary_5_10;
  TH1      * h_nBinary_10_20;
  TH1      * h_nBinary_20_30;
  TH1      * h_nBinary_30_40;
  TH1      * h_nBinary_40_50;
  TH1      * h_nBinary_50_60;
  TH1      * h_nBinary_60_70;
  TH1      * h_nBinary_70_80;
  TH1      * h_nBinary_80_90;
  TH1      * h_nBinary_90_100;


//  TGraph      * g_nPartRmsVsAvgNpart;
//  TGraph      * g_nBinaryRmsVsAvgNpart;
//  TGraph      * g_nPartOmegaVsAvgNpart;
//  TGraph      * g_nBinaryOmegaVsAvgNpart;
//  TGraph      * g_nPartR2VsAvgNpart;
//  TGraph      * g_nBinaryR2VsAvgNpart;

  TH3      * h_xyNNIntVsB;
  TProfile * h_varXVsB_Prof;
  TProfile * h_varYVsB_Prof;
  TProfile * h_covXYVsB_Prof;
  TProfile * h_epsilonXVsB_Prof;
  TProfile * h_epsilonYVsB_Prof;
  TProfile * h_epsilonXYVsB_Prof;
  TH3      * h_epsilonXYVsB;
  TProfile * h_psi2VsB_Prof;
  TH2      * h_psi2VsB;

  ClassDef(CollisionGeometryHistograms,0)

};

#endif /* WAC_CollisionGeometryHistograms  */



