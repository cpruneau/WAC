//
//  ParticleFilterAliceV0.cpp
//  MyMC
//
//  Created by Claude Pruneau on 12/12/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include <TMath.h>
#include "ParticleFilterAliceV0.hpp"

ClassImp(ParticleFilterAliceV0);

// ==========================================================================================
// CTOR
// Note: To turn off filtering on a specific variable, set the min large than the max.
// ==========================================================================================
ParticleFilterAliceV0::ParticleFilterAliceV0(V0Selection v0Requested,
                                             SpeciesSelection pidRequested,
                                             ChargeSelection  chargeRequested,
                                             double minPtRequested,
                                             double maxPtRequested)
:
ParticleFilter(pidRequested,
               chargeRequested,
               minPtRequested, maxPtRequested,
               1.0,    -1.0,
               1.0,    -1.0),
v0Selection(v0Requested),
min_eta_v0A(2.8),
max_eta_v0A(5.1),
min_eta_v0C(-3.7),
max_eta_v0C(-1.7)
{
  // no ops
}

// ==========================================================================================
// DTOR
// ==========================================================================================
ParticleFilterAliceV0::~ParticleFilterAliceV0()
{
  // no ops
}

// ==========================================================================================
// accept/reject the given particle based on filter parameter
// Filtering is based on
// Charge : enum ChargeSelection   { AllCharges, Negative, Positive, Charged, Neutral };
// Species: enum SpeciesSelection  { AllSpecies, Photon, Lepton, Electron, Muon, Hadron, Pion, Kaon, Baryon, Proton, Lambda };
// pt     : accept conditionally if min_pt < pt <= max_pt  OR  if min_pt >= max_pt
// eta    : accept conditionally if min_eta< eta<= max_eta OR  if min_eta>= max_eta
// y      : accept conditionally if   min_y< y  <= max_y OR    if min_y>  = max_y
// ==========================================================================================
bool ParticleFilterAliceV0::accept(Particle & particle)
{
  if (!acceptCharge(particle.charge))          return false;
  if (!acceptPid(TMath::Abs(particle.pid)))    return false;
  if (filterOnPt  && !acceptPt(particle.pt) )  return false;
  bool accepting = true;
  double eta = particle.eta;
  switch (v0Selection)
    {
      case V0A:   accepting = (min_eta_v0A<eta)&&(eta<= max_eta_v0A); break;
      case V0C:   accepting = (min_eta_v0C<eta)&&(eta<= max_eta_v0C); break;
      case V0M:   accepting = ((min_eta_v0A<eta)&&(eta<= max_eta_v0A)) || ((min_eta_v0C<eta)&&(eta<= max_eta_v0C));break;
    }
  return accepting;
}

// ==========================================================================================
// Creates a short filter name based on the PID and charge accepted
// ==========================================================================================
TString ParticleFilterAliceV0::getName()
{
  TString name;
  switch (v0Selection)
    {
      default: name = "WakyWaky"; break;
      case V0A:     name = "V0A"; break;
      case V0C:     name = "V0C"; break;
      case V0M:     name = "V0M"; break;
    }
  return name;
}

// ==========================================================================================
// Creates a short filter title  based on the PID and charge accepted
// ==========================================================================================
TString ParticleFilterAliceV0::getTitle()
{
  TString name;
  switch (v0Selection)
    {
      default: name = "WakyWaky"; break;
      case V0A:     name = "V0A"; break;
      case V0C:     name = "V0C"; break;
      case V0M:     name = "V0M"; break;
    }
  return name;
}



// ==========================================================================================
// Creates a long filter name based on the PID and charge accepted
// as well as the pT, eta, and y minimum and maximum accepted values.
// To avoid floating point values, all floats are multiplied by 1000.
// ==========================================================================================
TString ParticleFilterAliceV0::getLongName()
{
  TString name = getName();
  if (filterOnPt)
    {
    name += "PtGeq";
    name += int(1000.0*min_pt);
    name += "Lt";
    name += int(1000.0*max_pt);
    }
  if (filterOnEta)
    {
    name += "EtaGeq";
    name += int(1000.0*min_eta);
    name += "Lt";
    name += int(1000.0*max_eta);
    }
  if (filterOnY)
    {
    name += "YGeq";
    name += int(1000.0*min_y);
    name += "Lt";
    name += int(1000.0*max_y);
    }
  return name;
}


// ==========================================================================================
// Creates a long filter name based on the PID and charge accepted
// as well as the pT, eta, and y minimum and maximum accepted values.
// To avoid floating point values, all floats are multiplied by 1000.
// ==========================================================================================
TString ParticleFilterAliceV0::getLongTitle()
{
  TString name = getTitle();
  if (filterOnPt)   name += Form(" %g < p_{T} < %g;",min_pt,max_pt);
  if (filterOnEta)  name += Form(" %g < #eta < %g;",min_eta,max_eta);
  if (filterOnY)    name += Form(" %g < Y < %g",min_y,max_y);
  return name;
}
