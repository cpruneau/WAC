#include "AACollisionGenerator.hpp"
#include "TDatabasePDG.h"
#include "AnalysisConfiguration.hpp"

ClassImp(AACollisionGenerator);

AACollisionGenerator::AACollisionGenerator(const TString & name,
	TaskConfiguration * configuration,
	Event * event,
	EventFilter * ef,
	ParticleFilter * pf,
	CollisionGeometry * collisionGeo,
  int* maxCollisions)
:
Task(name, configuration, event),
eventFilter(ef),
particleFilter(pf),
nCollisions(40000),
collisionGeometry(collisionGeo),
nCollisionsMax(maxCollisions)
{
	if (reportDebug()) cout << "AACollisionGenerator::AACollisionGenerator(...) No ops" << endl;
}

AACollisionGenerator::~AACollisionGenerator()
{
	if (reportDebug()) cout << "AACollisionGenerator::~AACollisionGenerator(...) No ops" << endl;
}

void AACollisionGenerator::initialize()
{
	if (reportDebug()) cout << "AACollisionGenerator::initialize() Started" << endl;

  *nCollisionsMax = 0;

  nMax = 10000 ;
  particles = new TClonesArray("TParticle", nMax );
  pythia8 = new TPythia8();


  //pythia8->ReadString("Init:showMultipartonInteractions = off"); // don't list multiparton interaction initialization
  //pythia8->ReadString("Init:showProcesses = off");            // don't list initialization info
  pythia8->ReadString("Init:showChangedSettings = on");      // list changed settings
  pythia8->ReadString("Init:showChangedParticleData = off"); // list changed particle data
  pythia8->ReadString("Next:numberCount = 100000");            // print message every n events
  pythia8->ReadString("Next:numberShowInfo = 1");            // print event information n times
  pythia8->ReadString("Next:numberShowProcess = 0");         // print process record n times
  pythia8->ReadString("Next:numberShowEvent = 0");

  pythia8->ReadString("SoftQCD:all = on");                   // Allow total sigma = elastic/SD/DD/ND
                                                             // Optionally only study one or a few processes at a time.
                                                             //  pythia8->ReadString("SoftQCD:elastic = on");               // Elastic
                                                             //  pythia8->ReadString("SoftQCD:singleDiffractive = on");     // Single diffractive
                                                             //  pythia8->ReadString("SoftQCD:doubleDiffractive = on");     // Double diffractive
                                                             //  pythia8->ReadString("SoftQCD:centralDiffractive = on");    // Central diffractive
                                                             //  pythia8->ReadString("SoftQCD:nonDiffractive = on");        // Nondiffractive (inelastic)
                                                             //  pythia8->ReadString("SoftQCD:inelastic = on");             // All inelastic
                                                             // Optionally switch on hand-set cross section and Coulomb term.
                                                             // Note: these values are illustrative only, not to be taken seriously.
                                                             //  pythia8->ReadString("SigmaTotal:setOwn = on");              // set your own total/elastic/diffr. sigma
                                                             //  pythia8->ReadString("SigmaTotal:sigmaTot = 106.0");         // total cross section (excluding Coulomb)
                                                             //  pythia8->ReadString("SigmaTotal:sigmaXX = 5.0");        // reduce double diffractive cross section
                                                             //  pythia8->ReadString("SigmaElastic:setOwn = on");        // suboption for Coulomb term in elastic
                                                             //  pythia8->ReadString("SigmaElastic:tAbsMin = 4e-5");        // divergence requires lower cutoff

  pythia8->ReadString("HardQCD:all = on");


  //  Optionally select diffractive model.
  //  pythia8->ReadString("SigmaDiffractive:mode = 2");           // the MBR model for diffraction

  // Multiparton interactions and impact parameter picture.
  // Note: these values are illustrative only, not to be taken seriously.
  //  pythia8->ReadString("MultipartonInteractions:pT0Ref = 2.");
  //  pythia8->ReadString("MultipartonInteractions:ecmRef = 1960.");
  //  pythia8->ReadString("MultipartonInteractions:ecmPow = 0.16");
  //  pythia8->ReadString("MultipartonInteractions:pTmin = 0.2");
  //  pythia8->ReadString("MultipartonInteractions:bProfile = 2");
  //  pythia8->ReadString("MultipartonInteractions:coreRadius = 0.4");
  //  pythia8->ReadString("MultipartonInteractions:coreFraction = 0.5");
  //  pythia8->ReadString("MultipartonInteractions:expPow = 1.");

  // use a reproducible seed: always the same results for the tutorial.
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 42");

  // Initialize
  pythia8->Initialize(2212 /* p */, 2212 /* p */, 14000. /* TeV */);

  if (reportDebug()) cout << "AACollisionGenerator::initialize() Ended" << endl;


}

void AACollisionGenerator::execute()
{
	if (reportDebug()) cout << "AACollisionGenerator::execute() Started" << endl;

  nCollisions = collisionGeometry->nBinary;

  *nCollisionsMax = nCollisions> *nCollisionsMax? nCollisions: *nCollisionsMax;

  Factory<Particle> * particleFactory = Particle::getFactory();

  if(nCollisions != 0)
  {
    particleFactory -> initialize(Particle::factorySize * nCollisions);
  }

  if (reportDebug()) cout << "AACollisionGenerator::execute() processing " << nCollisions << " collisions." << endl;

  event->nParticles = 0;
  event->multiplicity = 0;

  for(int i = 0; i < nCollisions; i++)
  {
    int nparts;
    bool seekingEvent = true;
    while (seekingEvent)
    {
     pythia8->GenerateEvent();
     if (reportDebug()) pythia8->EventListing();
     if (reportDebug()) cout << "AACollisionGenerator::execute() Calling pythia8->ImportParticles()" << endl;

     pythia8->ImportParticles(particles,"Final");
     if (reportDebug()) cout << "AACollisionGenerator::execute() pythia8->ImportParticles() completed" << endl;

     nparts = particles->GetEntriesFast();
     if (reportDebug()) cout << "AACollisionGenerator::execute() with nparts:" << nparts << endl;

     if (nparts>2) seekingEvent = false;
   }

   if (nparts>nMax)
   {
     if (reportError()) cout << " ARRAY TOO SMALL np>nMax. nparts:" << nparts << " nMax:" << nMax << endl;
     postTaskFatal();
    //exit(0);
   }


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
  //if (reportDebug()) cout << "AACollisionGenerator::execute() starting copy loop into event..." << endl;

   for (int iParticle = 0; iParticle < nparts; iParticle++)
   {
     TParticle & part = * (TParticle*) particles->At(iParticle);
     int ist = part.GetStatusCode();
    //if (reportDebug()) cout << "AACollisionGenerator::execute() ist: " << ist << endl;
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
    		//if (reportDebug()) cout << "AACollisionGenerator::execute() calling filter " << endl;
    		particleCounted++;
    		if (!particleFilter->accept(aParticle)) continue;
    		particle = particleFactory->getNextObject();
    		*particle = aParticle;
    		particleAccepted++;
    		//    if (true)
    		//      {
    		//      cout << "AACollisionGenerator::execute() particle: " << iParticle << " / " << particleAccepted << endl;
    		//      particle->printProperties(cout);
    		//      }
    	}

    	event->nParticles += particleAccepted;
    	event->multiplicity += particleCounted;
    	if (reportDebug()) cout << "Generated Event " << eventsProcessed + 1 << ":" << i + 1 << endl;
    }
    eventsProcessed++;
    if (reportDebug()) cout << "AACollisionGenerator::execute() No of accepted Particles : "<< event->nParticles<<endl;
    if (reportDebug()) cout << "AACollisionGenerator::execute() No of counted Particles : "<< event->multiplicity <<endl;
    if (reportDebug()) cout << "AACollisionGenerator::execute() event completed!" << endl;
  }

  void AACollisionGenerator::finalize()
  {
   if (reportDebug()) cout << "AACollisionGenerator::finalize() started" << endl;
	if (reportInfo()) //pythia8->PrintStatistics();
	if (reportDebug()) cout << "AACollisionGenerator::finalize() completed" << endl;
}

void AACollisionGenerator::reset()
{
	if (reportDebug()) cout << "AACollisionGenerator::reset() Started" << endl;
	event->reset();
	Particle::getFactory()->reset();
	if (reportDebug()) cout << "AACollisionGenerator::reset() Completed" << endl;
}