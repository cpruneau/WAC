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
#include "PythiaEventReader.hpp"
ClassImp(PythiaEventReader);

PythiaEventReader::PythiaEventReader(const TString & name,
                                     TaskConfiguration * configuration,
                                     Event * event,
                                     EventFilter * ef,
                                     ParticleFilter * pf,
                                     LogLevel selectedLevel)
:
Task(name, configuration, event, selectedLevel),
fChain(nullptr),
fCurrent(0),
inputDataFile(nullptr),
nentries(0),
nbytes(0),
nb(0),
jentry(0),
eventFilter(ef),
particleFilter(pf)
{
  if (reportDebug()) cout << "PythiaEventReader::PythiaEventReader(...) No ops" << endl;
}

PythiaEventReader::~PythiaEventReader()
{
   cout << "PythiaEventReader::~PythiaEventReader(...)" << endl;
  if (inputDataFile)
    {
      inputDataFile->Close();
      delete inputDataFile;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize generator
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PythiaEventReader::initialize()
{
  if (reportStart("PythiaEventReader",getTaskName(),"initialize()"))
    ;

  PythiaConfiguration * pc = (PythiaConfiguration*) getTaskConfiguration();

  TString inputFileName = pc->dataInputPath;
  inputFileName += "/";
  inputFileName += pc->dataInputFileName;
  if (reportInfo("PythiaEventReader",getTaskName(),"initialize()")) cout << "Opening file: " << inputFileName << endl;
  inputDataFile = TFile::Open(inputFileName);
  if (!inputDataFile)
    {
    if (reportFatal("PythiaEventReader",getTaskName(),"initialize()")) cout << "Unable to open file: " << inputFileName <<endl;
    postTaskFatal();
    return;
    }
  TTree * inputTree = nullptr;
  inputDataFile->GetObject(pc->dataInputTreeName,inputTree);
  if (!inputTree)
    {
    if (reportFatal("PythiaEventReader",getTaskName(),"initialize()"))
      cout << "No inputTree named: " << pc->dataInputTreeName << " in file: " << inputFileName << endl;
    postTaskFatal();
    return;
    }
  Init(inputTree);
  if (reportDebug("PythiaEventReader",getTaskName(),"initialize()")) cout << "PythiaEventReader::initialize() Completed" << endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read an ampt event from file
// Copy the event into Event for convenience...
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PythiaEventReader::execute()
{
  //if (reportStart("PythiaEventReader",getTaskName(),"execute()"))
  //  ;
  event->reset();
  Factory<Particle> * particleFactory = Particle::getFactory();
  particleFactory->reset();
  int nparts;
  bool seekingEvent = true;
  //if (reportDebug("PythiaEventReader",getTaskName(),"execute()")) cout << "Start seek loop" << endl;

  while (seekingEvent)
    {
    //if (reportDebug("PythiaEventReader",getTaskName(),"execute()")) cout << "jentry:" << jentry << endl;
    // load another event from the root file/TTree
    Long64_t ientry = LoadTree(jentry++);
    //if (reportDebug("PythiaEventReader",getTaskName(),"execute()")) cout << "ientry:" << ientry << endl;

    // returning a null point is an indication that
    // there are no more events in the file or stack of files.
    if (ientry < 0)
      {
      postTaskEod(); // end of data
      return;
      }
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    nparts = particles_;
    //if (reportDebug("PythiaEventReader",getTaskName(),"execute()")) cout << " nb:" << nb << " nparts:" <<  nparts << endl;
    if (nparts>2) seekingEvent = false;
    }

  //if (reportDebug("PythiaEventReader",getTaskName(),"execute()")) cout << "Copy into event" << endl;

  int thePid;
  double charge, baryon, mass, p_x, p_y, p_z, p_e;
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
  //if (reportDebug("PythiaEventReader",getTaskName(),"execute()")) cout << "iParticle: " << iParticle << endl;

  int ist = particles_fStatusCode[iParticle];
  if (ist <= 0) continue;
  int pdg = particles_fPdgCode[iParticle];
  mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  if (mass<0.002) continue;  // no photons, electrons...
  charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/3.0;

  double px = particles_fPx[iParticle];
  double py = particles_fPy[iParticle];
  double pz = particles_fPz[iParticle];
  double e  = particles_fE[iParticle];
  p_x  = cosPhi*px - sinPhi*py;
  p_y  = sinPhi*px + cosPhi*py;
  p_z  = pz;
  p_e  = e;
  aParticle.setPidPxPyPzE(pdg, charge, p_x,p_y,p_z,p_e);

  int work = fabs(pdg);
  float baryon = ((work>=1000) && (work<6000))? 1 : 0;
  if (pdg<0) baryon  = -baryon;
  aParticle.baryon = baryon;
  //aParticle.printProperties(cout);
  //if (reportDebug()) cout << "PythiaEventGenerator::execute() calling filter " << endl;
  particleCounted++;
  if (!particleFilter->accept(aParticle)) continue;
  particle = particleFactory->getNextObject();
  *particle = aParticle;
  particleAccepted++;
  //if (reportDebug("PythiaEventReader",getTaskName(),"execute()")) cout << "particleAccepted: " << particleAccepted << endl;
}
  event->nParticles   = particleAccepted;
  event->multiplicity = particleCounted ;
//  if (reportDebug("PythiaEventReader",getTaskName(),"execute()"))
//    {
//    cout << endl;
//    cout << "No of accepted Particles : "<< particleAccepted<<endl;
//    cout << " No of counted Particles : "<< particleCounted <<endl;
//    }
}

Int_t PythiaEventReader::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t PythiaEventReader::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void PythiaEventReader::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("particles", &particles_, &b_particles_);
  fChain->SetBranchAddress("particles.fUniqueID", particles_fUniqueID, &b_particles_fUniqueID);
  fChain->SetBranchAddress("particles.fBits", particles_fBits, &b_particles_fBits);
  fChain->SetBranchAddress("particles.fLineColor", particles_fLineColor, &b_particles_fLineColor);
  fChain->SetBranchAddress("particles.fLineStyle", particles_fLineStyle, &b_particles_fLineStyle);
  fChain->SetBranchAddress("particles.fLineWidth", particles_fLineWidth, &b_particles_fLineWidth);
  fChain->SetBranchAddress("particles.fPdgCode", particles_fPdgCode, &b_particles_fPdgCode);
  fChain->SetBranchAddress("particles.fStatusCode", particles_fStatusCode, &b_particles_fStatusCode);
  fChain->SetBranchAddress("particles.fMother[2]", particles_fMother, &b_particles_fMother);
  fChain->SetBranchAddress("particles.fDaughter[2]", particles_fDaughter, &b_particles_fDaughter);
  fChain->SetBranchAddress("particles.fWeight", particles_fWeight, &b_particles_fWeight);
  fChain->SetBranchAddress("particles.fCalcMass", particles_fCalcMass, &b_particles_fCalcMass);
  fChain->SetBranchAddress("particles.fPx", particles_fPx, &b_particles_fPx);
  fChain->SetBranchAddress("particles.fPy", particles_fPy, &b_particles_fPy);
  fChain->SetBranchAddress("particles.fPz", particles_fPz, &b_particles_fPz);
  fChain->SetBranchAddress("particles.fE", particles_fE, &b_particles_fE);
  fChain->SetBranchAddress("particles.fVx", particles_fVx, &b_particles_fVx);
  fChain->SetBranchAddress("particles.fVy", particles_fVy, &b_particles_fVy);
  fChain->SetBranchAddress("particles.fVz", particles_fVz, &b_particles_fVz);
  fChain->SetBranchAddress("particles.fVt", particles_fVt, &b_particles_fVt);
  fChain->SetBranchAddress("particles.fPolarTheta", particles_fPolarTheta, &b_particles_fPolarTheta);
  fChain->SetBranchAddress("particles.fPolarPhi", particles_fPolarPhi, &b_particles_fPolarPhi);
  Notify();
  nentries = fChain->GetEntriesFast();
  nbytes = 0;
  nb = 0;
}

Bool_t PythiaEventReader::Notify()
{
  return kTRUE;
}

void PythiaEventReader::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t PythiaEventReader::Cut(Long64_t entry)
{
  return 1;
}


