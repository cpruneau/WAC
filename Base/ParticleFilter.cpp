//
//  ParticleFilter.cpp
//  MyMC
//
//  Created by Claude Pruneau on 12/12/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include <TMath.h>
#include "ParticleFilter.hpp"

ClassImp(ParticleFilter);

// ==========================================================================================
// CTOR
// Note: To turn off filtering on a specific variable, set the min large than the max.
// ==========================================================================================
ParticleFilter::ParticleFilter(ParticleFilter::SpeciesSelection pid,
                               ParticleFilter::ChargeSelection  charge,
                               double minPt,
                               double maxPt,
                               double minEta,
                               double maxEta,
                               double minY,
                               double maxY)
:
pidRequested     (pid),
chargeRequested  (charge),
filterOnPt       (maxPt  > minPt),
filterOnEta      (maxEta > minEta),
filterOnY        (maxY   > minY),
min_pt           (minPt),
max_pt           (maxPt),
min_eta          (minEta),
max_eta          (maxEta),
min_y            (minY),
max_y            (maxY)
{
  // no ops
}

// ==========================================================================================
// DTOR
// ==========================================================================================
ParticleFilter::~ParticleFilter()
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
bool ParticleFilter::accept(Particle & particle)
{
  bool accepting = true;
  double charge = particle.charge;
  switch (chargeRequested)
    {
      case AllCharges: accepting = true; break;
      case Negative:   accepting = (charge < 0); break;
      case Positive:   accepting = (charge > 0); break;
      case Charged:    accepting = (charge != 0); break;
      case Neutral:    accepting = (charge == 0); break;
    }
  if (!accepting) return false;
  double pid = TMath::Abs(particle.pid);
  switch (pidRequested)
    {
      case AllSpecies:   accepting = true; break;
      case Photon:   accepting = (pid == 22); break;
      case Lepton:   accepting = (pid == 11) || (pid == 12)|| (pid == 13)|| (pid == 14)|| (pid == 15)|| (pid == 16)  ; break;
      case Electron: accepting = (pid == 11); break;
      case Muon:     accepting = (pid == 13); break;
      case Hadron:  accepting = (pid == 211) || (pid == 321) || (pid == 2212); break;
      case Pion:    accepting = (pid == 111) || (pid == 211); break;
      case Kaon:    accepting = (pid == 321) || (pid == 311) || (pid == 310) || (pid == 130); break;
      case Baryon:  accepting = (pid == 2212) || (pid == 2112) || (pid == 3122); break;
      case Proton:  accepting = (pid == 2212); break;
      case Lambda:  accepting = (pid == 3122); break;
    }
  if (!accepting) return false;

  if (filterOnPt)
    {
    double pt  = particle.pt;
    accepting = (min_pt<pt) &&  (pt<= max_pt);
    }
  if (!accepting) return false;

  if (filterOnEta)
    {
    double eta = particle.eta;
    accepting = (min_eta<eta) &&  (eta<= max_eta);
    }
  if (!accepting) return false;

  if (filterOnY)
    {
    double y = particle.y;
    accepting = (min_y<y) &&  (y<= max_y);
    }

  return accepting;
}

// ==========================================================================================
// Creates a short filter name based on the PID and charge accepted
// ==========================================================================================
TString ParticleFilter::getName()
{
  TString name;
  switch (pidRequested)
    {
      default: name = "WakyWaky"; break;
      case AllSpecies: name = "All"; break;
      case Photon:     name = "Gam"; break;
      case Lepton:     name = "Lep"  ; break;
      case Electron:   name = "El"; break;
      case Muon:       name = "Mu"; break;
      case Hadron:     name = "H"; break;
      case Pion:       name = "Pi"; break;
      case Kaon:       name = "K"; break;
      case Baryon:     name = "B"; break;
      case Proton:     name = "P"; break;
      case Lambda:     name = "Lam"; break;
    }

  if (pidRequested == Photon) return name;

  switch (chargeRequested)
     {
       default: name = "WakyWakyWakyWaky"; break;
       case AllCharges:   name += "All"; break; // all charges
       case Negative: name += "M"; break; // minus
       case Positive: name += "P"; break;
       case Charged:  name += "C"; break;
       case Neutral:  name += "0"; break;
     }
  return name;
}

// ==========================================================================================
// Creates a short filter title  based on the PID and charge accepted
// ==========================================================================================
TString ParticleFilter::getTitle()
{
  TString name;
  switch (pidRequested)
    {
      case AllSpecies:name = "l,h"; break;
      case Photon:   name = "#gamma"; break;
      case Lepton:   name = "l"  ; break;
      case Electron: name = "e"; break;
      case Muon:     name = "#muon"; break;
      case Hadron:   name = "h"; break;
      case Pion:     name = "#pi"; break;
      case Kaon:     name = "K"; break;
      case Baryon:   name = "Baryon"; break;
      case Proton:   name = "p"; break;
      case Lambda:   name = "#Lambda"; break;
    }

  if (pidRequested == Photon) return name;


  switch (chargeRequested)
    {
      case AllCharges:name += "^{0,#pm}"; break; // all charges
      case Negative: name += "^{-}"; break; // minus
      case Positive: name += "^{+}"; break;
      case Charged:  name += "^{#pm}"; break;
      case Neutral:  name += "^{0}"; break;
    }
  return name;
}

// ==========================================================================================
// Creates a long filter name based on the PID and charge accepted
// as well as the pT, eta, and y minimum and maximum accepted values.
// To avoid floating point values, all floats are multiplied by 1000.
// ==========================================================================================
TString ParticleFilter::getLongName()
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
TString ParticleFilter::getLongTitle()
{
  TString name = getTitle();
  if (filterOnPt)   name += Form(" %g < p_{T} < %g;",min_pt,max_pt);
  if (filterOnEta)  name += Form(" %g < #eta < %g;",min_eta,max_eta);
  if (filterOnY)    name += Form(" %g < Y < %g",min_y,max_y);
  return name;
}
