// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_ParticleFilter
#define WAC_ParticleFilter
#include "TString.h"
#include "Particle.hpp"

//////////////////////////////////////////////////////////////////////////////////////////
// Single Particle Filter
//
// chargeRequested:
// case -1:    accepts negative only
// case  0:    accepts neutral only
// case  1:    accepts positive only
// case  999:  accepts all
//////////////////////////////////////////////////////////////////////////////////////////



class ParticleFilter
{
public:

  enum ChargeSelection   { AllCharges, Negative, Positive, Charged, Neutral };
  enum SpeciesSelection  { AllSpecies, Photon, Lepton, Electron, Muon, Hadron, Pion, Kaon, Baryon, Proton, Lambda };

  ParticleFilter(SpeciesSelection pidRequested,
                  ChargeSelection  chargeRequested,
                 double minPt,
                 double maxPt,
                 double minEta,
                 double maxEta,
                 double minY,
                 double maxY);
  virtual ~ParticleFilter();
  virtual bool accept(Particle & particle);
  virtual bool acceptCharge(double charge);
  virtual bool acceptPid(double pid);
  inline bool acceptPt(double pt)
  {
    return (min_pt<pt) &&  (pt<= max_pt);
  }
  inline bool acceptEta(double eta)
  {
    return (min_eta<eta) && (eta<= max_eta);
  }
  inline bool acceptY(double y)
  {
    return (min_y<y) && (y<= max_y);
  }
  virtual TString getName();
  virtual TString getTitle();
  virtual TString getLongName();
  virtual TString getLongTitle();

  //////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////////////////////////////////

  SpeciesSelection pidRequested;
  ChargeSelection  chargeRequested;
  bool   filterOnPt;
  bool   filterOnEta;
  bool   filterOnY;
  double min_pt;
  double max_pt;
  double min_eta;
  double max_eta;
  double min_y;
  double max_y;


  ClassDef(ParticleFilter,0)
};

#endif /* WAC_ParticleFilter */
