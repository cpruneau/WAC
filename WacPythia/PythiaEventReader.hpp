// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_PythiaEventReader
#define WAC_PythiaEventReader
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"

class PythiaEventReader : public Task
{
public:

  PythiaEventReader(const TString & name,
                    TaskConfiguration * configuration,
                    Event * event,
                    EventFilter * ef,
                    ParticleFilter * pf,
                    LogLevel selectedLevel);
  virtual ~PythiaEventReader();
  virtual void initialize();
  void execute();

  Int_t GetEntry(Long64_t entry);
  Long64_t LoadTree(Long64_t entry);
  void Init(TTree *tree);
  Bool_t Notify();
  void Show(Long64_t entry);
  Int_t Cut(Long64_t entry);

  ////////////////////////////////////////////////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////////////////////////////////////////////////
  TTree  *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t   fCurrent; //!current Tree number in a TChain
  TFile  *inputDataFile;
  Long64_t nentries;
  Long64_t nbytes;
  Long64_t nb;
  long jentry;

  static const int kMaxparticles = 2000;

  // Declaration of leaf types
  Int_t           particles_;
  UInt_t          particles_fUniqueID[kMaxparticles];   //[particles_]
  UInt_t          particles_fBits[kMaxparticles];   //[particles_]
  Short_t         particles_fLineColor[kMaxparticles];   //[particles_]
  Short_t         particles_fLineStyle[kMaxparticles];   //[particles_]
  Short_t         particles_fLineWidth[kMaxparticles];   //[particles_]
  Int_t           particles_fPdgCode[kMaxparticles];   //[particles_]
  Int_t           particles_fStatusCode[kMaxparticles];   //[particles_]
  Int_t           particles_fMother[kMaxparticles][2];   //[particles_]
  Int_t           particles_fDaughter[kMaxparticles][2];   //[particles_]
  Float_t         particles_fWeight[kMaxparticles];   //[particles_]
  Double_t        particles_fCalcMass[kMaxparticles];   //[particles_]
  Double_t        particles_fPx[kMaxparticles];   //[particles_]
  Double_t        particles_fPy[kMaxparticles];   //[particles_]
  Double_t        particles_fPz[kMaxparticles];   //[particles_]
  Double_t        particles_fE[kMaxparticles];   //[particles_]
  Double_t        particles_fVx[kMaxparticles];   //[particles_]
  Double_t        particles_fVy[kMaxparticles];   //[particles_]
  Double_t        particles_fVz[kMaxparticles];   //[particles_]
  Double_t        particles_fVt[kMaxparticles];   //[particles_]
  Double_t        particles_fPolarTheta[kMaxparticles];   //[particles_]
  Double_t        particles_fPolarPhi[kMaxparticles];   //[particles_]

  // List of branches
  TBranch        *b_particles_;   //!
  TBranch        *b_particles_fUniqueID;   //!
  TBranch        *b_particles_fBits;   //!
  TBranch        *b_particles_fLineColor;   //!
  TBranch        *b_particles_fLineStyle;   //!
  TBranch        *b_particles_fLineWidth;   //!
  TBranch        *b_particles_fPdgCode;   //!
  TBranch        *b_particles_fStatusCode;   //!
  TBranch        *b_particles_fMother;   //!
  TBranch        *b_particles_fDaughter;   //!
  TBranch        *b_particles_fWeight;   //!
  TBranch        *b_particles_fCalcMass;   //!
  TBranch        *b_particles_fPx;   //!
  TBranch        *b_particles_fPy;   //!
  TBranch        *b_particles_fPz;   //!
  TBranch        *b_particles_fE;   //!
  TBranch        *b_particles_fVx;   //!
  TBranch        *b_particles_fVy;   //!
  TBranch        *b_particles_fVz;   //!
  TBranch        *b_particles_fVt;   //!
  TBranch        *b_particles_fPolarTheta;   //!
  TBranch        *b_particles_fPolarPhi;   //!

  EventFilter * eventFilter;
  ParticleFilter * particleFilter;

  ClassDef(PythiaEventReader,0)
};

#endif /* WAC_PythiaEventReader */
