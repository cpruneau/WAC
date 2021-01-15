// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_GlobalAnalyzerConfiguration
#define WAC_GlobalAnalyzerConfiguration
#include "TaskConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class GlobalAnalyzerConfiguration : public TaskConfiguration
{
public:
  
  GlobalAnalyzerConfiguration(const TString & name,
                              const TString & type,
                              const TString & version);
  GlobalAnalyzerConfiguration(const GlobalAnalyzerConfiguration & source);
  virtual ~GlobalAnalyzerConfiguration(){}

  GlobalAnalyzerConfiguration & operator=(const GlobalAnalyzerConfiguration & source);

  virtual void validate();
  virtual void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  int nBins_n;  int nBins_n2;  double min_n;  double max_n;  double range_n;
  int nBins_e;  int nBins_e2;  double min_e;  double max_e;  double range_e;
  int nBins_q;  int nBins_q2;  double min_q;  double max_q;  double range_q;
  int nBins_b;  int nBins_b2;  double min_b;  double max_b;  double range_b;

  bool countParticles;
  bool setEvent;
  bool fillCorrelationHistos;

  ClassDef(GlobalAnalyzerConfiguration,0)
};

#endif /* WAC_GlobalAnalyzerConfiguration */
