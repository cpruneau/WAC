// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_PythiaConfiguration
#define WAC_PythiaConfiguration
#include "TaskConfiguration.hpp"

/*
pythia8->ReadString("Init:showChangedSettings = on");      // list changed settings
pythia8->ReadString("Init:showChangedParticleData = off"); // list changed particle data
pythia8->ReadString("Next:numberCount = 100000");            // print message every n events
pythia8->ReadString("Next:numberShowInfo = 1");            // print event information n times
pythia8->ReadString("Next:numberShowProcess = 0");         // print process record n times
pythia8->ReadString("Next:numberShowEvent = 0");

pythia8->ReadString("SoftQCD:all = off");                   // Allow total sigma = elastic/SD/DD/ND
                                                           // Optionally only study one or a few processes at a time.
pythia8->ReadString("SoftQCD:elastic = on");                 Elastic
pythia8->ReadString("SoftQCD:singleDiffractive = on");      Single diffractive
pythia8->ReadString("SoftQCD:doubleDiffractive = on");      Double diffractive
pythia8->ReadString("SoftQCD:centralDiffractive = on");     Central diffractive
pythia8->ReadString("SoftQCD:nonDiffractive = on");         Nondiffractive (inelastic)
pythia8->ReadString("SoftQCD:inelastic = on");              All inelastic
                                                            Optionally switch on hand-set cross section and Coulomb term.
                                                            Note: these values are illustrative only, not to be taken seriously.
pythia8->ReadString("SigmaTotal:setOwn = on");               set your own total/elastic/diffr. sigma
pythia8->ReadString("SigmaTotal:sigmaTot = 106.0");          total cross section (excluding Coulomb)
pythia8->ReadString("SigmaTotal:sigmaXX = 5.0");         reduce double diffractive cross section
                                                             pythia8->ReadString("SigmaElastic:setOwn = on");         suboption for Coulomb term in elastic
                                                             pythia8->ReadString("SigmaElastic:tAbsMin = 4e-5");         divergence requires lower cutoff

pythia8->ReadString("HardQCD:all = on");
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////
// PYTHIA Generator Configuration
////////////////////////////////////////////////////////////////////////////////////////////////////////
class PythiaConfiguration : public TaskConfiguration
{
public:
  
  PythiaConfiguration(int    beam,
                      int    target,
                      double energy,
                      int    nOptions,
                      TString ** options);
  virtual ~PythiaConfiguration(){}
  void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  int    beam;  // PDG Code   proton is 2212
  int    target;
  double energy; // in GeV
  int    nOptions;
  TString ** options;



ClassDef(PythiaConfiguration,0)
};

#endif /* WAC_PythiaConfiguration */
