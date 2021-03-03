#include <ROOT/RDataFrame.hxx>
#include "AACollisionReader.hpp"
#include "TDatabasePDG.h"
#include "Nucleon.hpp"
#include "Particle.hpp"

ClassImp(AACollisionReader);

AACollisionReader::AACollisionReader(const TString &name,
                                     TaskConfiguration *configuration,
                                     Event *event,
                                     EventFilter *ef,
                                     ParticleFilter *pf,
                                     LogLevel selectedLevel,
                                     CollisionGeometry *collisionGeo)
    : PythiaEventReader(name, configuration, event, ef, pf, selectedLevel),
      nCollisions(40000),
      collisionGeometry(collisionGeo),
      ppOnly(true)
{
  if (reportDebug())
    cout << "AACollisionReader::AACollisionReader(...) No ops" << endl;
}

AACollisionReader::~AACollisionReader()
{
  if (reportDebug())
    cout << "AACollisionReader::~AACollisionReader(...) No ops" << endl;
  delete fChain;
}
void AACollisionReader::initialize()
{
  if (reportStart("void AACollisionReader::initialize()", getTaskName(), "initialize()"))
    ;
  PythiaConfiguration *pc = (PythiaConfiguration *)getTaskConfiguration();
  removePhotons = pc->removePhotons;
  TString inputFileName = pc->dataInputPath;
  inputFileName += "/";
  inputFileName += pc->dataInputFileName;
  if (reportInfo("AACollisionReader", getTaskName(), "initialize()"))
    cout << "Opening file: " << inputFileName << endl;
  inputDataFile = TFile::Open(inputFileName);
  if (!inputDataFile)
  {
    if (reportFatal("AACollisionReader", getTaskName(), "initialize()"))
      cout << "Unable to open file: " << inputFileName << endl;
    postTaskFatal();
    return;
  }
  for (int i = 0; i < 3; i++)
  {
    TTree *inputTree = nullptr;
    inputDataFile->GetObject(pc->ppdataInputTreeName[i], inputTree);
    if (!inputTree)
    {
      if (!ppOnly && i > 0)
      {
        if (reportFatal("AACollisionReader", getTaskName(), "initialize()"))
          cout << "No inputTree named: " << pc->ppdataInputTreeName[i] << " in file: " << inputFileName << endl;
        postTaskFatal();
        return;
      }
      continue;
    }
    Init(i, inputTree);
  }

  if (reportDebug("AACollisionReader", getTaskName(), "initialize()"))
    cout << "AACollisionReader::initialize() Completed" << endl;

  Factory<Particle> *particleFactory = Particle::getFactory();
  particleFactory->initialize(Particle::factorySize * 2000);
}

void AACollisionReader::execute()
{
  Factory<Particle> *particleFactory = Particle::getFactory();
  event->reset();
  //particleFactory->reset();
  int thePid;
  double charge, mass, p_x, p_y, p_z, p_e;
  Particle *particle;
  Particle aParticle;
  int nParticles;
  int particleAccepted = 0;
  int particleCounted = 0;
  int nCollisions = collisionGeometry->nBinary; //get the number of binary collisions
  int pp = 0;                                   // corresponds to pp,pn,np,or nn collision
  if (reportDebug("AACollisionPythiaGenerator", getTaskName(), "execute()"))
    cout << "nCollisions:" << nCollisions << endl;
  for (int iCollision = 0; iCollision < nCollisions; iCollision++)
  {
    if (!ppOnly)
    {
      CollisionGeometry::NucleonNucleonType nnType = collisionGeometry->getNNTpye(iCollision);
      switch (nnType)
      {
      case CollisionGeometry::ProtonProton:
        pp = 0;
        break;
      case CollisionGeometry::ProtonNeutron:
        pp = 1;
        break;
      case CollisionGeometry::NeutronProton:
        pp = 2;
        break;
      case CollisionGeometry::NeutronNeutron:
        pp = 3;
        break;
      }
    }

    bool seekingEvent = true;
    while (seekingEvent)
    {
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "jentry:" << jentry << endl;
      // load another event from the root file/TTree
      Long64_t ientry = LoadTree(pp, ppjentry[pp]++);
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "ientry:" << ientry << endl;

      /*if (ientry < 0) // returning a null point is an indication that there are no more events in the file or stack of files.
      {
        postTaskEod(); // end of data
        return;
      }*/
      if (ientry < 0)
      {
        ppjentry[pp] = 0; // move back to the start of the file
      }
      nb = ppfChain[pp]->GetEntry(ppjentry[pp]);
      nbytes += nb;
      nParticles = particles_;
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << " nb:" << nb << " nparts:" <<  nparts << endl;
      if (nParticles > 2)
        seekingEvent = false;
    }
    double eventAngle = TMath::TwoPi() * gRandom->Rndm();
    double cosPhi = cos(eventAngle);
    double sinPhi = sin(eventAngle);
    if (reportDebug("AACollisionReader", getTaskName(), "execute()"))
      cout << "Starting copy loop into event..." << endl;
    for (int iParticle = 0; iParticle < nParticles; iParticle++)
    {
      int ist = particles_fStatusCode[iParticle];
      if (ist <= 0)
        continue;
      int pdg = particles_fPdgCode[iParticle];
      mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
      if (mass < 0.0001)
        continue; // no photons, electrons...
      charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge() / 3.0;
      double px = particles_fPx[iParticle];
      double py = particles_fPy[iParticle];
      double pz = particles_fPz[iParticle];
      double e = particles_fE[iParticle];
      p_x = cosPhi * px - sinPhi * py;
      p_y = sinPhi * px + cosPhi * py;
      p_z = pz;
      p_e = e;
      aParticle.setPidPxPyPzE(pdg, charge, p_x, p_y, p_z, p_e);
      aParticle.setSourceIndex(iCollision);
      int work = fabs(pdg);
      float baryon = ((work >= 1000) && (work < 6000)) ? 1 : 0;
      if (pdg < 0)
        baryon = -baryon;
      aParticle.baryon = baryon;
      particleCounted++;
      if (!particleFilter->accept(aParticle))
        continue;
      particle = particleFactory->getNextObject();
      *particle = aParticle;
      particleAccepted++;
    }

    event->nParticles += particleAccepted;
    event->multiplicity += particleCounted;
    nEventProcessed++;
    nEventAccepted++;

    if (reportDebug("AACollisionReader", getTaskName(), "execute()"))
    {
      cout << "No of accepted Particles : " << particleAccepted << endl;
      cout << " No of counted Particles : " << particleCounted << endl;
    }
  }
}

Long64_t AACollisionReader::LoadTree(int i, Long64_t entry)
{
  // Set the environment to read one entry
  if (!ppfChain[i])
    return -5;
  Long64_t centry = ppfChain[i]->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (ppfChain[i]->GetTreeNumber() != ppfCurrent[i])
  {
    ppfCurrent[i] = ppfChain[i]->GetTreeNumber();
    Notify();
  }
  return centry;
}

void AACollisionReader::Init(int i, TTree *tree)
{
  if (!tree)
    return;
  ppfChain[i] = tree;
  ppfCurrent[i] = -1;
  ppfChain[i]->SetMakeClass(1);

  ppfChain[i]->SetBranchAddress("particles", &particles_, &b_particles_);
  ppfChain[i]->SetBranchAddress("particles.fUniqueID", particles_fUniqueID, &b_particles_fUniqueID);
  ppfChain[i]->SetBranchAddress("particles.fBits", particles_fBits, &b_particles_fBits);
  ppfChain[i]->SetBranchAddress("particles.fLineColor", particles_fLineColor, &b_particles_fLineColor);
  ppfChain[i]->SetBranchAddress("particles.fLineStyle", particles_fLineStyle, &b_particles_fLineStyle);
  ppfChain[i]->SetBranchAddress("particles.fLineWidth", particles_fLineWidth, &b_particles_fLineWidth);
  ppfChain[i]->SetBranchAddress("particles.fPdgCode", particles_fPdgCode, &b_particles_fPdgCode);
  ppfChain[i]->SetBranchAddress("particles.fStatusCode", particles_fStatusCode, &b_particles_fStatusCode);
  ppfChain[i]->SetBranchAddress("particles.fMother[2]", particles_fMother, &b_particles_fMother);
  ppfChain[i]->SetBranchAddress("particles.fDaughter[2]", particles_fDaughter, &b_particles_fDaughter);
  ppfChain[i]->SetBranchAddress("particles.fWeight", particles_fWeight, &b_particles_fWeight);
  ppfChain[i]->SetBranchAddress("particles.fCalcMass", particles_fCalcMass, &b_particles_fCalcMass);
  ppfChain[i]->SetBranchAddress("particles.fPx", particles_fPx, &b_particles_fPx);
  ppfChain[i]->SetBranchAddress("particles.fPy", particles_fPy, &b_particles_fPy);
  ppfChain[i]->SetBranchAddress("particles.fPz", particles_fPz, &b_particles_fPz);
  ppfChain[i]->SetBranchAddress("particles.fE", particles_fE, &b_particles_fE);
  ppfChain[i]->SetBranchAddress("particles.fVx", particles_fVx, &b_particles_fVx);
  ppfChain[i]->SetBranchAddress("particles.fVy", particles_fVy, &b_particles_fVy);
  ppfChain[i]->SetBranchAddress("particles.fVz", particles_fVz, &b_particles_fVz);
  ppfChain[i]->SetBranchAddress("particles.fVt", particles_fVt, &b_particles_fVt);
  ppfChain[i]->SetBranchAddress("particles.fPolarTheta", particles_fPolarTheta, &b_particles_fPolarTheta);
  ppfChain[i]->SetBranchAddress("particles.fPolarPhi", particles_fPolarPhi, &b_particles_fPolarPhi);
  Notify();
  ppnentries[i] = ppfChain[i]->GetEntriesFast();
  ppnbytes[i] = 0;
  ppnb[i] = 0;
}