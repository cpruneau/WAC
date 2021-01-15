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
      collisionGeometry(collisionGeo)
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

void AACollisionReader::execute()
{
  if (reportDebug())
    cout << "AACollisionReader::execute() Started" << endl;


  nCollisions = collisionGeometry->nBinary; //get the number of binary collisions
  if (reportDebug())
    cout << "AACollisionReader::execute() processing " << nCollisions << " collisions." << endl;

  if (!fChain)
  {
    if (reportFatal())
      cout << " AACollisionReader::execute() no TChain available" << endl;
    postTaskFatal();
    return;
  }
  event->reset();
  Factory<Particle> *particleFactory = Particle::getFactory();
  particleFactory->reset();


  for (int i = 0; i < nCollisions; i++)
  {

    ///////////////////////////////////////////////
    // read events and move particles into the particle factory
    //////////////////////////////////////////////

    //if (reportStart("AACollisionReader",getTaskName(),"execute()"))
    //  ;

    //Following chunk of code comes from the Superclass execute()
    int nparts;
    bool seekingEvent = true;
    //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "Start seek loop" << endl;

    while (seekingEvent)
    {
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "jentry:" << jentry << endl;
      // load another event from the root file/TTree
      Long64_t ientry = LoadTree(jentry++);
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "ientry:" << ientry << endl;

      // returning a null point is an indication that
      // there are no more events in the file or stack of files.
      /*if (ientry < 0)
      {
        postTaskEod(); // end of data
        return;
      }*/
      if(ientry < 0)
      {
        jentry = 0; // move back to the start of the file
      }
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      nparts = particles_;
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << " nb:" << nb << " nparts:" <<  nparts << endl;
      if (nparts > 2)
        seekingEvent = false;
    }

    //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "Copy into event" << endl;

    int thePid;
    double charge, baryon, mass, p_x, p_y, p_z, p_e;
    Particle *particle;
    int particleAccepted = 0;
    int particleCounted = 0;

    //------------------- Randomizing the particle phi --------------Starts
    double eventAngle = TMath::TwoPi() * gRandom->Rndm();
    double cosPhi = cos(eventAngle);
    double sinPhi = sin(eventAngle);

    // load particles from TClone storage and copy into event.
    Particle aParticle;
    //if (reportDebug()) cout << "PythiaEventGenerator::execute() starting copy loop into event..." << endl;

    for (int iParticle = 0; iParticle < nparts; iParticle++)
    {
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "iParticle: " << iParticle << endl;

      int ist = particles_fStatusCode[iParticle];
      if (ist <= 0)
        continue;
      int pdg = particles_fPdgCode[iParticle];
      mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
      if (mass < 0.002)
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

      int work = fabs(pdg);
      float baryon = ((work >= 1000) && (work < 6000)) ? 1 : 0;
      if (pdg < 0)
        baryon = -baryon;
      aParticle.baryon = baryon;
      //aParticle.printProperties(cout);
      //if (reportDebug()) cout << "PythiaEventGenerator::execute() calling filter " << endl;
      particleCounted++;
      if (!particleFilter->accept(aParticle))
        continue;
      particle = particleFactory->getNextObject();
      *particle = aParticle;
      particleAccepted++;
      //if (reportDebug("AACollisionReader",getTaskName(),"execute()")) cout << "particleAccepted: " << particleAccepted << endl;
    }

    event->nParticles += particleAccepted;
    event->multiplicity += particleCounted;

    //  if (reportDebug("AACollisionReader",getTaskName(),"execute()"))
    //    {
    //    cout << endl;
    //    cout << "No of accepted Particles : "<< particleAccepted<<endl;
    //    cout << " No of counted Particles : "<< particleCounted <<endl;
    //    }
    //if (reportDebug())
    //  cout << "Generated Event " << eventsProcessed + 1 << ":" << i + 1 << endl;
  }
  nEventProcessed++;
  if (reportDebug())
    cout << "AACollisionReader::execute() No of accepted Particles : " << event->nParticles << endl;
  if (reportDebug())
    cout << "AACollisionReader::execute() No of counted Particles : " << event->multiplicity << endl;
  if (reportDebug())
    cout << "AACollisionReader::execute() event completed!" << endl;
}
