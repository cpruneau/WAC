// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_PythiaEventGenerator
#define WAC_PythiaEventGenerator
#include "TParticle.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TFile.h"
#include "TTree.h"
#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"

class PythiaEventGenerator : public Task
{
public:


  PythiaEventGenerator(const TString & name,
                       TaskConfiguration * configuration,
                       Event * event,
                       EventFilter * ef,
                       ParticleFilter * pf,
                       LogLevel selectedLevel);
  virtual ~PythiaEventGenerator();
  virtual void initialize();
  virtual void finalize();
  void execute();

  TPythia8* pythia8; // = new TPythia8();

  // For WAC analyses
  int nMax; //  = 10000;
  TClonesArray* particles; // = new TClonesArray("TParticle", nMax);
  EventFilter * eventFilter;
  ParticleFilter * particleFilter;

  // For TTree file output
  // Set up the ROOT TFile and TTree.
  TFile *outputFile;
  Pythia8::Event *outputEvent;
  TTree *outputTree;

  ClassDef(PythiaEventGenerator,0)
};

#endif /* WAC_PythiaEventGenerator */
