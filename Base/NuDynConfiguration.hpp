// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_NuDynConfiguration
#define WAC_NuDynConfiguration
#include "TMath.h"
#include "TaskConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class NuDynConfiguration : public TaskConfiguration
{
public:

  enum MultiplicityType   { Centrality, TotalMultiplicity, AcceptedMultiplicity};

  NuDynConfiguration(const TString & name,
                        const TString & type,
                        const TString & version);
  NuDynConfiguration(const NuDynConfiguration & source);
  virtual ~NuDynConfiguration();
  NuDynConfiguration & operator=(const NuDynConfiguration & source);
  void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////

  MultiplicityType multiplicityType;
  int nBins_mult;  double min_mult; double max_mult;

ClassDef(NuDynConfiguration,0)
};

#endif /* WAC_NuDynConfiguration */
