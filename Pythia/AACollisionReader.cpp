#include <ROOT/RDataFrame.hxx>
#include "AACollisionReader.hpp"
#include "TDatabasePDG.h"
#include "AnalysisConfiguration.hpp"
#include "HeavyIonConfiguration.hpp"
#include "Nucleon.hpp"
#include "Particle.hpp"

ClassImp(AACollisionReader);

AACollisionReader::AACollisionReader(const TString & name,
  TaskConfiguration * configuration,
  Event * event,
  EventFilter * ef,
  ParticleFilter * pf,
  CollisionGeometry * collisionGeo)
:
PythiaEventReader(name, configuration, event, ef, pf),
nCollisions(40000),
collisionGeometry(collisionGeo)
{
  if (reportDebug()) cout << "AACollisionReader::AACollisionReader(...) No ops" << endl;
}

AACollisionReader::~AACollisionReader()
{
  if (reportDebug()) cout << "AACollisionReader::~AACollisionReader(...) No ops" << endl;
  delete fChain;
}

void AACollisionReader::initialize()
{
  if (reportDebug()) cout << "AACollisionReader::initialize() Started" << endl;
  int nthreads = 5;
  ROOT::EnableImplicitMT(nthreads); // to speed up getting the entry
  AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
  HeavyIonConfiguration * hc = (HeavyIonConfiguration *) ac;

  fileLocation = hc->treeFile;
  eventTreeName = hc->eventTreeName;
  branchName_eventNo= hc->branchName_eventNo;
  branchName_mult= hc->branchName_mult;
  branchName_px= hc->branchName_px;
  branchName_py= hc->branchName_py;
  branchName_pz= hc->branchName_pz;
  branchName_ist= hc->branchName_ist;
  branchName_pdg= hc->branchName_pdg;
  branchName_pE= hc->branchName_pE;
  numFiles = hc->numFiles;

  TChain *chain = new TChain(eventTreeName);
  for(int i=0; i<numFiles; i++)
  {
    chain->Add(Form(fileLocation + "%i.root",i));
  }
  Init(chain);
  jentry = 0;
  nentries = fChain->GetEntries();
  if (reportInfo()) cout << "AACollisionReader::initialize() nEntries: " << nentries << endl;
  if (nentries < 1)
  {
    if (reportError()) cout << "AACollisionReader::initialize() no data found. Abort." << endl;
    postTaskFatal();
    return;
  }
  nbytes = 0;
  nb = 0;
Factory<Particle> * particleFactory = Particle::getFactory();
  particleFactory -> initialize(Particle::factorySize * 2000);
  if (reportDebug()) cout << "AACollisionReader::initialize() Ended" << endl;


}

void AACollisionReader::execute()
{
  if (reportDebug()) cout << "AACollisionReader::execute() Started" << endl;



  AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
  HeavyIonConfiguration * hc = (HeavyIonConfiguration *) ac;

  nCollisions = collisionGeometry->nBinary; //get the number of binary collisions
  if (reportDebug()) cout << "AACollisionReader::execute() processing " << nCollisions << " collisions." << endl;

  if (!fChain)
  {
    if (reportFatal()) cout << " AACollisionReader::execute() no TChain available" << endl;
    postTaskFatal();
    return;
  }

  Factory<Particle> * particleFactory = Particle::getFactory();



  hc->nCollisionsMax = nCollisions> hc->nCollisionsMax? nCollisions: hc->nCollisionsMax; //set the max number of binary collisions per event, to set the size of the histos later

  
  event->nParticles = 0;
  event->multiplicity = 0;


  for(int i = 0; i < nCollisions; i++)
  {


    ///////////////////////////////////////////////
    // read events and move particles into the particle factory
    //////////////////////////////////////////////


    bool seekingEvent = true;
    while (seekingEvent)
    {
      // load a random event from the root file/TTree
      jentry = gRandom->Integer(nentries);
      Long64_t ientry = LoadTree(jentry);

    // returning a null point is an indication that
    // there are no more events in the file or stack of files.
      //if (ientry < 0)
      //{
      //postTaskEod(); // end of data
      //return;
    //}
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      seekingEvent = false;
  }

  if (mult > arraySize)
  {
    if (reportError()) cout<< "AACollisionReader::execute() n particles is " << mult << " and exceeds capacity " << arraySize << endl;
    postTaskError();
    return;
  }

   //////////////////////////////////////////////////////////////////////
   // Calculate the boost for the particles
   //////////////////////////////////////////////////////////////////////
  double x_col, y_col, z_col;
  x_col = collisionGeometry->x[i];
  y_col = collisionGeometry->y[i];
  z_col = collisionGeometry->z[i];

  double transverseR = TMath::Sqrt(x_col*x_col + y_col*y_col);
  double phi = TMath::ATan(y_col/x_col);
  if(x_col < 0) phi += TMath::Pi()/2;
   double param_b = hc->param_b; // exponent of order 1
   double param_a = hc->param_a;
   double beta = param_a * TMath::Power(transverseR, param_b);
   if(beta > 1 - TMath::Power(10, -8)) beta = 1 - TMath::Power(10, -8);
   double betax = beta * TMath::Cos(phi);
   double betay = beta * TMath::Sin(phi);

  ///////////////////////////////////////////////////////////////////////////////////////// 
  // load particles from tree storage and copy into event.
  /////////////////////////////////////////////////////////////////////////////////////////


   double charge, mass, p_x, p_y, p_z, p_e;
   int particleAccepted = 0;
   int particleCounted = 0;
   Particle * particle;

  //------------------- Randomizing the particle phi --------------Starts
   double eventAngle= TMath::TwoPi() * gRandom->Rndm();
   double cosPhi = cos(eventAngle);
   double sinPhi = sin(eventAngle);

   int iParticle=0;

   if (reportDebug()) cout << "AACollisionReader::execute() starting copy loop into event..." << endl;
   Particle aParticle;
   for (int iParticle = 0; iParticle < mult; iParticle++)
   {
     int _ist = ist[iParticle];
    //if (reportDebug()) cout << "AACollisionGenerator::execute() ist: " << ist << endl;
     if (_ist <= 0) continue;
     int _pdg = pdg[iParticle];
     mass = TDatabasePDG::Instance()->GetParticle(_pdg)->Mass();
        if (mass<0.002) continue;  // no photons, electrons...
        charge = TDatabasePDG::Instance()->GetParticle(_pdg)->Charge();
        p_x  = cosPhi*px[iParticle] - sinPhi*py[iParticle];
        p_y  = sinPhi*px[iParticle] + cosPhi*py[iParticle];
        p_z  = pz[iParticle];
        p_e  = pE[iParticle];
        aParticle.setPidPxPyPzE(_pdg, charge, p_x,p_y,p_z,p_e);


        aParticle.boost(betax,betay,0.0);

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
    if (reportDebug()) cout << "AACollisionReader::execute() No of accepted Particles : "<< event->nParticles<<endl;
    if (reportDebug()) cout << "AACollisionReader::execute() No of counted Particles : "<< event->multiplicity <<endl;
    if (reportDebug()) cout << "AACollisionReader::execute() event completed!" << endl;
  }

  void AACollisionReader::Init(TTree *tree)
  {
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    fChain->SetBranchAddress(branchName_eventNo, &eventNo);
    fChain->SetBranchAddress(branchName_mult, &mult);
    fChain->SetBranchAddress(branchName_px, px);
    fChain->SetBranchAddress(branchName_py, py);
    fChain->SetBranchAddress(branchName_pz, pz);
    fChain->SetBranchAddress(branchName_ist, ist);
    fChain->SetBranchAddress(branchName_pdg, pdg);
    fChain->SetBranchAddress(branchName_pE, pE);
    Notify();
    nentries = fChain->GetEntriesFast();
    nbytes = 0;
    nb = 0;
  }

