// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */
#include "TDatabasePDG.h"
#include "PythiaEventGenerator.hpp"
#include "PythiaConfiguration.hpp"

ClassImp(PythiaEventGenerator);

PythiaEventGenerator::PythiaEventGenerator(const TString & name,
                                           TaskConfiguration * configuration,
                                           Event * event,
                                           EventFilter * ef,
                                           ParticleFilter * pf,
                                           LogLevel selectedLevel)
:
Task(name, configuration, event, selectedLevel),
pythia8(nullptr),
nMax(10000),
particles(nullptr),
eventFilter(ef),
particleFilter(pf),
outputFile(nullptr),
outputEvent(nullptr),
outputTree(nullptr)
{
  if (reportDebug()) cout << "PythiaEventGenerator::PythiaEventGenerator(...) No ops" << endl;
}

PythiaEventGenerator::~PythiaEventGenerator()
{
  if (reportDebug()) cout << "PythiaEventGenerator::~PythiaEventGenerator(...) No ops" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize generator
//pythia8->Initialize(2212 /* p */, 2212 /* p */, 14000. /* GeV */);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PythiaEventGenerator::initialize()
{
  if (reportDebug()) cout << "PythiaEventGenerator::initialize() Started" << endl;


  PythiaConfiguration * pc = (PythiaConfiguration*) getTaskConfiguration();
  pythia8 = new TPythia8();
  for (int iOption=0; iOption<pc->nOptions; iOption++)
  {
  pythia8->ReadString( *pc->options[iOption]);
  }
  pythia8->Initialize(pc->beam,pc->target,pc->energy);

  if (pc->dataOutputUsed)
    {
    TString outputFileName = pc->dataOutputPath;
    outputFileName += "/";
    outputFileName += pc->dataOutputFileName;
    outputFile = TFile::Open(outputFileName,"recreate");
    //outputEvent = &pythia8->Pythia8()->event;
    outputTree  = new TTree(pc->dataOutputTreeName,"PythiaEventTree");
    particles = (TClonesArray*) pythia8->GetListOfParticles();
    outputTree->Branch("particles", &particles);
    //outputTree->Branch("event",&outputEvent);
    }
  if (pc->dataConversionToWac)
    {
    nMax = 10000;
    particles = new TClonesArray("TParticle", nMax);
    }

  if (reportDebug()) cout << "PythiaEventGenerator::initialize() Completed" << endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read an ampt event from file
// Copy the event into Event for convenience...
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PythiaEventGenerator::execute()
{
  if (reportDebug()) cout << "PythiaEventGenerator::execute() Started" << endl;
  PythiaConfiguration * pc = (PythiaConfiguration*) getTaskConfiguration();

  int nparts;
  bool seekingEvent = true;
  while (seekingEvent)
    {
    pythia8->GenerateEvent();
    if (reportDebug()) pythia8->EventListing();
    if (reportDebug()) cout << "PythiaEventGenerator::execute() Calling pythia8->ImportParticles()" << endl;
    pythia8->ImportParticles(particles,"Final");
    if (reportDebug()) cout << "PythiaEventGenerator::execute() pythia8->ImportParticles() completed" << endl;
    nparts = particles->GetEntriesFast();
    if (reportDebug()) cout << "PythiaEventGenerator::execute() with nparts:" << nparts << endl;
    if (nparts>2) seekingEvent = false;
    }
  if (nparts>nMax)
    {
    if (reportError()) cout << " ARRAY TOO SMALL np>nMax. nparts:" << nparts << " nMax:" << nMax << endl;
    postTaskFatal();
    //exit(0);
    }

  if (pc->dataOutputUsed)
    {
    outputTree->Fill();
    }
  if (pc->dataConversionToWac)
    {
    Factory<Particle> * particleFactory = Particle::getFactory();
    particleFactory->reset();
    int thePid;
    double charge, mass, p_x, p_y, p_z, p_e;
    Particle * particle;
    int particleAccepted = 0;
    int particleCounted = 0;

    //------------------- Randomizing the particle phi --------------Starts
    double eventAngle= TMath::TwoPi() * gRandom->Rndm();
    double cosPhi = cos(eventAngle);
    double sinPhi = sin(eventAngle);

    // load particles from TClone storage and copy into event.
    Particle aParticle;
    //if (reportDebug()) cout << "PythiaEventGenerator::execute() starting copy loop into event..." << endl;

    for (int iParticle = 0; iParticle < nparts; iParticle++)
      {
      TParticle & part = * (TParticle*) particles->At(iParticle);
      int ist = part.GetStatusCode();
      //if (reportDebug()) cout << "PythiaEventGenerator::execute() ist: " << ist << endl;
      if (ist <= 0) continue;
      int pdg = part.GetPdgCode();
      mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
      if (mass<0.002) continue;  // no photons, electrons...
      charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      p_x  = cosPhi*part.Px() - sinPhi*part.Py();
      p_y  = sinPhi*part.Px() + cosPhi*part.Py();
      p_z  = part.Pz();
      p_e  = part.Energy();
      aParticle.setPidPxPyPzE(pdg, charge, p_x,p_y,p_z,p_e);
      //aParticle.printProperties(cout);
      //if (reportDebug()) cout << "PythiaEventGenerator::execute() calling filter " << endl;
      particleCounted++;
      if (!particleFilter->accept(aParticle)) continue;
      particle = particleFactory->getNextObject();
      *particle = aParticle;
      particleAccepted++;
      //    if (true)
      //      {
      //      cout << "PythiaEventGenerator::execute() particle: " << iParticle << " / " << particleAccepted << endl;
      //      particle->printProperties(cout);
      //      }
      }

    event->nParticles   = particleCounted;
    event->multiplicity = particleAccepted;
    if (reportDebug()) cout << "PythiaEventGenerator::execute() No of accepted Particles : "<< particleAccepted<<endl;
    if (reportDebug()) cout << "PythiaEventGenerator::execute() No of counted Particles : "<< particleCounted <<endl;
    }
}

void PythiaEventGenerator::finalize()
{
  if (reportDebug()) cout << "PythiaEventGenerator::finalize() started" << endl;
  if (reportInfo()) pythia8->PrintStatistics();
  PythiaConfiguration * pc = (PythiaConfiguration*) getTaskConfiguration();
  if (pc->dataOutputUsed)
    {
    outputTree->Print();
    outputTree->Write();
    delete outputFile;
    }
  if (reportDebug()) cout << "PythiaEventGenerator::finalize() completed" << endl;
}




//pythia8->ReadString("Init:showChangedSettings = on");      // list changed settings
//pythia8->ReadString("Init:showChangedParticleData = off"); // list changed particle data
//pythia8->ReadString("Next:numberCount = 100000");            // print message every n events
//pythia8->ReadString("Next:numberShowInfo = 1");            // print event information n times
//pythia8->ReadString("Next:numberShowProcess = 0");         // print process record n times
//pythia8->ReadString("Next:numberShowEvent = 0");
//
//pythia8->ReadString("SoftQCD:all = off");                   // Allow total sigma = elastic/SD/DD/ND
//                                                           // Optionally only study one or a few processes at a time.
//                                                           //  pythia8->ReadString("SoftQCD:elastic = on");               // Elastic
//                                                           //  pythia8->ReadString("SoftQCD:singleDiffractive = on");     // Single diffractive
//                                                           //  pythia8->ReadString("SoftQCD:doubleDiffractive = on");     // Double diffractive
//                                                           //  pythia8->ReadString("SoftQCD:centralDiffractive = on");    // Central diffractive
//                                                           //  pythia8->ReadString("SoftQCD:nonDiffractive = on");        // Nondiffractive (inelastic)
//                                                           //  pythia8->ReadString("SoftQCD:inelastic = on");             // All inelastic
//                                                           // Optionally switch on hand-set cross section and Coulomb term.
//                                                           // Note: these values are illustrative only, not to be taken seriously.
//                                                           //  pythia8->ReadString("SigmaTotal:setOwn = on");              // set your own total/elastic/diffr. sigma
//                                                           //  pythia8->ReadString("SigmaTotal:sigmaTot = 106.0");         // total cross section (excluding Coulomb)
//                                                           //  pythia8->ReadString("SigmaTotal:sigmaXX = 5.0");        // reduce double diffractive cross section
//                                                           //  pythia8->ReadString("SigmaElastic:setOwn = on");        // suboption for Coulomb term in elastic
//                                                           //  pythia8->ReadString("SigmaElastic:tAbsMin = 4e-5");        // divergence requires lower cutoff
//
//pythia8->ReadString("HardQCD:all = on");
//
//
////  Optionally select diffractive model.
////  pythia8->ReadString("SigmaDiffractive:mode = 2");           // the MBR model for diffraction
//
//// Multiparton interactions and impact parameter picture.
//// Note: these values are illustrative only, not to be taken seriously.
////  pythia8->ReadString("MultipartonInteractions:pT0Ref = 2.");
////  pythia8->ReadString("MultipartonInteractions:ecmRef = 1960.");
////  pythia8->ReadString("MultipartonInteractions:ecmPow = 0.16");
////  pythia8->ReadString("MultipartonInteractions:pTmin = 0.2");
////  pythia8->ReadString("MultipartonInteractions:bProfile = 2");
////  pythia8->ReadString("MultipartonInteractions:coreRadius = 0.4");
////  pythia8->ReadString("MultipartonInteractions:coreFraction = 0.5");
////  pythia8->ReadString("ultipartonInteractions:expPow = 1.");
//
//// use a reproducible seed: always the same results for the tutorial.
////pythia8->ReadString("Random:setSeed = on");
////pythia8->ReadString("Random:seed = 42");
