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
  CollisionGeometryHistograms(TFile * inputFile,
                              const TString & collectionName,
                              CollisionGeometryConfiguration * configuration,
                              LogLevel  debugLevel);
  virtual ~CollisionGeometryHistograms() { }
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void fill(CollisionGeometry * collisionGeometry, double weight);


  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////

  CollisionGeometryConfiguration * configuration;

  TH1      * h_b;
  TH1      * h_nPart;
  TH1      * h_nBinary;
  TProfile * h_nPartVsB_Prof;
  TProfile * h_nBinaryVsB_Prof;
  TH2      * h_nPartVsB;
  TH2      * h_nBinaryVsB;
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



