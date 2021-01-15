// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_ParticleFilterAliceV0
#define WAC_ParticleFilterAliceV0
#include "TString.h"
#include "ParticleFilter.hpp"

//////////////////////////////////////////////////////////////////////////////////////////
// Single Particle Filter
//
// chargeRequested:
// case -1:    accepts negative only
// case  0:    accepts neutral only
// case  1:    accepts positive only
// case  999:  accepts all
//////////////////////////////////////////////////////////////////////////////////////////



class ParticleFilterAliceV0 : public ParticleFilter
{
public:

  enum V0Selection   { V0A, V0C, V0M };

  ParticleFilterAliceV0(V0Selection v0Selection=V0M,
                        SpeciesSelection pidRequested=AllSpecies,
                        ChargeSelection  chargeRequested=Charged,
                        double minPt=0.050,
                        double maxPt=1000.0);
  virtual ~ParticleFilterAliceV0();
  virtual bool accept(Particle & particle);

  virtual TString getName();
  virtual TString getTitle();
  virtual TString getLongName();
  virtual TString getLongTitle();

  //////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////////////////////////////////
  V0Selection v0Selection;
  double min_eta_v0A;
  double max_eta_v0A;
  double min_eta_v0C;
  double max_eta_v0C;

  ClassDef(ParticleFilterAliceV0,0)
};

#endif /* WAC_ParticleFilterAliceV0 */
