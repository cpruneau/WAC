// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_ParticlePairAnalyzerConfiguration
#define WAC_ParticlePairAnalyzerConfiguration
#include "ParticleAnalyzerConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class ParticlePairAnalyzerConfiguration : public ParticleAnalyzerConfiguration
{
public:
  
  ParticlePairAnalyzerConfiguration(const TString & name,
                        const TString & type,
                        const TString & version);
  ParticlePairAnalyzerConfiguration(const ParticlePairAnalyzerConfiguration & source);
  virtual ~ParticlePairAnalyzerConfiguration(){}

  ParticlePairAnalyzerConfiguration & operator=(const ParticlePairAnalyzerConfiguration & source);

  virtual void validate();
  virtual void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////

  // pair multiplicity bin correction
  float binCorrPM;
  float binCorrMM;
  float binCorrPP;

  // pairs q3d
  int nBins_DeltaPlong;  double min_DeltaPlong;  double max_DeltaPlong; double range_DeltaPlong;
  int nBins_DeltaPside;  double min_DeltaPside;  double max_DeltaPside; double range_DeltaPside;
  int nBins_DeltaPout;   double min_DeltaPout;   double max_DeltaPout;  double range_DeltaPout;
  int nBins_Dphi;        double min_Dphi;        double max_Dphi;       double width_Dphi;
  int nBins_Dphi_shft;   double min_Dphi_shft;   double max_Dphi_shft;
  int nBins_Deta;        double min_Deta;        double max_Deta;       double width_Deta;
  int nBins_Dy;          double min_Dy;          double max_Dy;         double width_Dy;

  bool fillPairs;
  bool fill6D;
  bool fillQ3D;

ClassDef(ParticlePairAnalyzerConfiguration,0)
};

#endif /* WAC_ParticlePairAnalyzerConfiguration */
