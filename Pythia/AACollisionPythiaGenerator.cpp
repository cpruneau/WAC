#include "AACollisionPythiaGenerator.hpp"
#include "TDatabasePDG.h"

ClassImp(AACollisionPythiaGenerator);

AACollisionPythiaGenerator::AACollisionPythiaGenerator(const TString & _name,
                                                       PythiaConfiguration * _configuration,
                                                       CollisionGeometry * _collisionGeometry,
                                                       Event * _event,
                                                       EventFilter * _eventFilter,
                                                       ParticleFilter * _particleFilter,
                                                       LogLevel _selectedLevel)
:
Task(_name, _configuration, _event, _selectedLevel),
collisionGeometry(_collisionGeometry),
eventFilter(_eventFilter),
particleFilter(_particleFilter),
ppPythia(nullptr),
pnPythia(nullptr),
npPythia(nullptr),
nnPythia(nullptr),
particles(nullptr),
clonesArraySize(0),
ppOnly(true)
{
  // no ops
}

AACollisionPythiaGenerator::~AACollisionPythiaGenerator()
{
  // no ops
}

void AACollisionPythiaGenerator::initialize()
{
  if (reportStart("AACollisionPythiaGenerator",getTaskName(),"initialize()"))
    ;

  PythiaConfiguration * pc = (PythiaConfiguration*) getTaskConfiguration();
  clonesArraySize  = pc->clonesArraySize;
  removePhotons    = pc->removePhotons;
  ppOnly           = pc->ppOnly;
  int protonId  = 2212;
  int neutronId = 2112;
  particles = new TClonesArray("TParticle", clonesArraySize );
  ppPythia = new TPythia8();
  if (!ppOnly)
    {
    pnPythia = new TPythia8();
    npPythia = new TPythia8();
    nnPythia = new TPythia8();
    }
  for (int iOption=0; iOption<pc->nOptions; iOption++)
  {
  TString option = *pc->options[iOption];
  if (reportInfo()) cout << "AACollisionPythiaGenerator::initialize() Read string:" << option << endl;
  ppPythia->ReadString( option );
  if (!ppOnly)
    {
    pnPythia->ReadString( option );
    npPythia->ReadString( option );
    nnPythia->ReadString( option );
    }
  }
  ppPythia->Initialize(protonId,protonId,pc->energy);
  if (!ppOnly)
    {
    pnPythia->Initialize(protonId,  neutronId, pc->energy);
    npPythia->Initialize(neutronId, protonId,  pc->energy);
    nnPythia->Initialize(neutronId, neutronId, pc->energy);
    }
  Factory<Particle> * particleFactory = Particle::getFactory();
  particleFactory -> initialize(Particle::factorySize * 5000);
  if (reportEnd("AACollisionPythiaGenerator",getTaskName(),"initialize()"))
    ;
}

void AACollisionPythiaGenerator::execute()
{
  Factory<Particle> * particleFactory = Particle::getFactory();
  event->reset();
  particleFactory->reset();
  int thePid;
  double charge, mass, p_x, p_y, p_z, p_e;
  Particle * particle;
  Particle aParticle;
  int nParticles;
  int particleAccepted = 0;
  int particleCounted = 0;
  TPythia8 * pythia = ppPythia;
  int nCollisions = collisionGeometry->nBinary; //get the number of binary collisions
  if (reportDebug("AACollisionPythiaGenerator",getTaskName(),"execute()")) cout << "nCollisions:" << nCollisions << endl;
  for (int iCollision = 0; iCollision < nCollisions; iCollision++)
  {
  if (!ppOnly)
    {
    CollisionGeometry::NucleonNucleonType nnType = collisionGeometry->getNNTpye(iCollision);
    switch (nnType)
      {
        case CollisionGeometry::ProtonProton   : pythia = ppPythia; break;
        case CollisionGeometry::ProtonNeutron  : pythia = pnPythia; break;
        case CollisionGeometry::NeutronProton  : pythia = npPythia; break;
        case CollisionGeometry::NeutronNeutron : pythia = nnPythia; break;
      }
    }
  bool seekingEvent = true;
  while (seekingEvent)
    {
    pythia->GenerateEvent();
    pythia->ImportParticles(particles,"Final");
    nParticles = particles->GetEntriesFast();
    if (reportDebug("AACollisionPythiaGenerator",getTaskName(),"execute()")) cout << "nParticles:" << nParticles << endl;
    if (nParticles>2) seekingEvent = false;
    }
  if (nParticles>clonesArraySize)
    {
    if (reportError("AACollisionPythiaGenerator",getTaskName(),"execute()")) cout << " nParticles:" << nParticles << " clonesArraySize:" << clonesArraySize << endl;
    postTaskFatal();
    }
  double eventAngle= TMath::TwoPi() * gRandom->Rndm();
  double cosPhi = cos(eventAngle);
  double sinPhi = sin(eventAngle);
  if (reportDebug("AACollisionPythiaGenerator",getTaskName(),"execute()")) cout << "Starting copy loop into event..." << endl;
  for (int iParticle = 0; iParticle < nParticles; iParticle++)
    {
    TParticle & part = * (TParticle*) particles->At(iParticle);
    int ist = part.GetStatusCode();
    //if (reportDebug()) cout << "AACollisionPythiaGenerator::execute() ist: " << ist << endl;
    if (ist <= 0) continue;
    int pdg = part.GetPdgCode();
    mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
    if (removePhotons && mass<0.002) continue;
    charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    p_x  = cosPhi*part.Px() - sinPhi*part.Py();
    p_y  = sinPhi*part.Px() + cosPhi*part.Py();
    p_z  = part.Pz();
    p_e  = part.Energy();
    aParticle.setPidPxPyPzE(pdg, charge, p_x,p_y,p_z,p_e);
    aParticle.setSourceIndex(iCollision);
    particleCounted++;
    if (!particleFilter->accept(aParticle)) continue;
    particle = particleFactory->getNextObject();
    *particle = aParticle;
    particleAccepted++;
    }
  event->nParticles   += particleAccepted;
  event->multiplicity += particleCounted;
  //if (reportDebug()) cout << "Generated Event " << eventsProcessed + 1 << ":" << i + 1 << endl;
  }
  incrementEventProcessed();
  //if (!eventFilter->accept(*event)) return;
  incrementEventAccepted(); // count events used to fill histograms and for scaling at the end...

  if (reportDebug("AACollisionPythiaGenerator",getTaskName(),"execute()"))
    {
    cout << "No of accepted Particles : "<< event->nParticles<<endl;
    cout << "No of counted Particles : "<< event->multiplicity <<endl;
    }
}

void AACollisionPythiaGenerator::finalize()
{
  if (reportInfo("AACollisionPythiaGenerator",getTaskName(),"finalize()"))
    {
    ppPythia->PrintStatistics();
    if (!ppOnly)
      {
      pnPythia->PrintStatistics();
      ppPythia->PrintStatistics();
      nnPythia->PrintStatistics();
      }
    }
  if (reportEnd("AACollisionPythiaGenerator",getTaskName(),"finalize()"))
    ;
}
