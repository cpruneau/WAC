// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_ParticleAnalyzerConfiguration
#define WAC_ParticleAnalyzerConfiguration
#include "TaskConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class ParticleAnalyzerConfiguration : public TaskConfiguration
{
public:
  
  ParticleAnalyzerConfiguration(const TString & name,
                                const TString & type,
                                const TString & version);
  ParticleAnalyzerConfiguration(const ParticleAnalyzerConfiguration & source);
  virtual ~ParticleAnalyzerConfiguration(){}

  ParticleAnalyzerConfiguration & operator=(const ParticleAnalyzerConfiguration & source);

  virtual void validate();
  virtual void printConfiguration(ostream & os);
  int getIxEtaPhi(double eta, double phi);
  int getIxYPhi(double y, double phi);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  TString histoBaseName;
  // single particle
  int nBins_pt;  double min_pt;  double max_pt;  double range_pt;
  int nBins_eta; double min_eta; double max_eta; double range_eta;
  int nBins_y;   double min_y;   double max_y;   double range_y;
  int nBins_phi; double min_phi; double max_phi; double range_phi;
  int nBins_phiEta;
  int nBins_phiEtaPt;
  int nBins_phiY;
  int nBins_phiYPt;

  bool fill3D;
  bool fillY;

  ClassDef(ParticleAnalyzerConfiguration,0)
};

#endif /* WAC_ParticleAnalyzerConfiguration */
